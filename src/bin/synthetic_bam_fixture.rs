use sequenza_utils::errors::{AppError, Result};
use clap::{ArgAction, Parser};
use std::ffi::{OsStr, OsString};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode, Output, Stdio};

#[derive(Debug, Parser)]
#[command(
    name = "synthetic_bam_fixture",
    about = "Generate deterministic normal/tumor BAM fixtures from one input BAM"
)]
struct Cli {
    #[arg(long = "input-bam")]
    input_bam: PathBuf,
    #[arg(long = "output-dir")]
    output_dir: PathBuf,
    #[arg(long = "sample-prefix")]
    sample_prefix: String,
    #[arg(long = "normal-fraction", value_parser = parse_fraction)]
    normal_fraction: f64,
    #[arg(long = "tumor-fraction", value_parser = parse_fraction)]
    tumor_fraction: f64,
    #[arg(long = "seed", default_value_t = 20260218)]
    seed: u64,
    #[arg(long = "threads", default_value_t = 1, value_parser = parse_threads)]
    threads: usize,
    #[arg(long = "samtools", default_value = "samtools")]
    samtools: String,
    #[arg(long = "region", action = ArgAction::Append)]
    region: Vec<String>,
    #[arg(long = "force", action = ArgAction::SetTrue)]
    force: bool,
}

#[derive(Debug)]
struct OutputPaths {
    normal_bam: PathBuf,
    normal_bai: PathBuf,
    tumor_bam: PathBuf,
    tumor_bai: PathBuf,
    normal_flagstat: PathBuf,
    normal_idxstats: PathBuf,
    tumor_flagstat: PathBuf,
    tumor_idxstats: PathBuf,
    manifest: PathBuf,
}

impl OutputPaths {
    fn from_args(args: &Cli) -> Self {
        let normal_bam = args
            .output_dir
            .join(format!("{}.normal.bam", args.sample_prefix));
        let tumor_bam = args
            .output_dir
            .join(format!("{}.tumor.bam", args.sample_prefix));
        Self {
            normal_bai: bai_path_for(&normal_bam),
            tumor_bai: bai_path_for(&tumor_bam),
            normal_flagstat: args
                .output_dir
                .join(format!("{}.normal.flagstat.txt", args.sample_prefix)),
            normal_idxstats: args
                .output_dir
                .join(format!("{}.normal.idxstats.txt", args.sample_prefix)),
            tumor_flagstat: args
                .output_dir
                .join(format!("{}.tumor.flagstat.txt", args.sample_prefix)),
            tumor_idxstats: args
                .output_dir
                .join(format!("{}.tumor.idxstats.txt", args.sample_prefix)),
            manifest: args
                .output_dir
                .join(format!("{}.manifest.tsv", args.sample_prefix)),
            normal_bam,
            tumor_bam,
        }
    }

    fn all_paths(&self) -> [&Path; 9] {
        [
            &self.normal_bam,
            &self.normal_bai,
            &self.tumor_bam,
            &self.tumor_bai,
            &self.normal_flagstat,
            &self.normal_idxstats,
            &self.tumor_flagstat,
            &self.tumor_idxstats,
            &self.manifest,
        ]
    }
}

#[derive(Debug, Clone)]
struct ManifestRow {
    role: &'static str,
    bam_path: PathBuf,
    bai_path: PathBuf,
    fraction: f64,
    seed: u64,
}

fn parse_fraction(raw: &str) -> std::result::Result<f64, String> {
    let parsed = raw
        .parse::<f64>()
        .map_err(|_| "must be a number > 0 and <= 1".to_string())?;
    if !(parsed > 0.0 && parsed <= 1.0) {
        return Err("must be > 0 and <= 1".to_string());
    }
    Ok(parsed)
}

fn parse_threads(raw: &str) -> std::result::Result<usize, String> {
    let parsed = raw
        .parse::<usize>()
        .map_err(|_| "must be a positive integer".to_string())?;
    if parsed == 0 {
        return Err("must be a positive integer".to_string());
    }
    Ok(parsed)
}

fn main() -> ExitCode {
    let args = Cli::parse();
    match run(args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("synthetic_bam_fixture: {error}");
            ExitCode::from(1)
        }
    }
}

fn run(args: Cli) -> Result<()> {
    validate_input_bam(&args.input_bam)?;
    ensure_command_exists(&args.samtools)?;
    ensure_command_exists("sha256sum")?;

    if args.output_dir.exists() && !args.output_dir.is_dir() {
        return Err(AppError::InvalidValue {
            flag: "--output-dir".to_string(),
            value: args.output_dir.display().to_string(),
            reason: "path exists but is not a directory".to_string(),
        });
    }
    fs::create_dir_all(&args.output_dir)?;

    let outputs = OutputPaths::from_args(&args);
    ensure_overwrite_policy(&outputs, args.force)?;

    let tumor_seed = args
        .seed
        .checked_add(1)
        .ok_or_else(|| AppError::InvalidValue {
            flag: "--seed".to_string(),
            value: args.seed.to_string(),
            reason: "seed is too large to derive tumor seed (seed + 1)".to_string(),
        })?;

    process_role(
        &args,
        "normal",
        args.normal_fraction,
        args.seed,
        &outputs.normal_bam,
        &outputs.normal_flagstat,
        &outputs.normal_idxstats,
    )?;
    process_role(
        &args,
        "tumor",
        args.tumor_fraction,
        tumor_seed,
        &outputs.tumor_bam,
        &outputs.tumor_flagstat,
        &outputs.tumor_idxstats,
    )?;

    let normal_manifest_row = ManifestRow {
        role: "normal",
        bam_path: outputs.normal_bam.clone(),
        bai_path: outputs.normal_bai.clone(),
        fraction: args.normal_fraction,
        seed: args.seed,
    };
    let tumor_manifest_row = ManifestRow {
        role: "tumor",
        bam_path: outputs.tumor_bam.clone(),
        bai_path: outputs.tumor_bai.clone(),
        fraction: args.tumor_fraction,
        seed: tumor_seed,
    };
    write_manifest(
        &outputs.manifest,
        &[normal_manifest_row, tumor_manifest_row],
    )?;

    Ok(())
}

fn validate_input_bam(input_bam: &Path) -> Result<()> {
    if input_bam.is_file() {
        return Ok(());
    }
    Err(AppError::InvalidValue {
        flag: "--input-bam".to_string(),
        value: input_bam.display().to_string(),
        reason: "input BAM does not exist or is not a regular file".to_string(),
    })
}

fn ensure_overwrite_policy(outputs: &OutputPaths, force: bool) -> Result<()> {
    if force {
        return Ok(());
    }
    for output in outputs.all_paths() {
        if output.exists() {
            return Err(AppError::InvalidValue {
                flag: "--force".to_string(),
                value: "false".to_string(),
                reason: format!("output already exists: {}", output.display()),
            });
        }
    }
    Ok(())
}

fn process_role(
    args: &Cli,
    role: &str,
    fraction: f64,
    seed: u64,
    output_bam: &Path,
    flagstat_path: &Path,
    idxstats_path: &Path,
) -> Result<()> {
    run_view_sort(
        &args.samtools,
        args.threads,
        &args.input_bam,
        output_bam,
        &args.region,
        fraction,
        seed,
    )?;
    run_samtools_index(&args.samtools, args.threads, output_bam)?;
    run_samtools_quickcheck(&args.samtools, output_bam)?;
    write_samtools_report(&args.samtools, "flagstat", output_bam, flagstat_path)?;
    write_samtools_report(&args.samtools, "idxstats", output_bam, idxstats_path)?;

    if !output_bam.is_file() {
        return Err(AppError::ParseError {
            message: format!(
                "{role} BAM was expected but not found after processing: {}",
                output_bam.display()
            ),
        });
    }
    let output_bai = bai_path_for(output_bam);
    if !output_bai.is_file() {
        return Err(AppError::ParseError {
            message: format!(
                "{role} BAM index was expected but not found after processing: {}",
                output_bai.display()
            ),
        });
    }
    Ok(())
}

fn run_view_sort(
    samtools: &str,
    threads: usize,
    input_bam: &Path,
    output_bam: &Path,
    regions: &[String],
    fraction: f64,
    seed: u64,
) -> Result<()> {
    let threads_arg = threads.to_string();
    let mut view = Command::new(samtools);
    view.arg("view").arg("-@").arg(&threads_arg).arg("-bh");

    if !fraction_is_one(fraction) {
        let sample_arg = samtools_sample_arg(seed, fraction);
        view.arg("-s").arg(sample_arg);
    }

    view.arg(input_bam);
    for region in regions {
        view.arg(region);
    }
    view.stdout(Stdio::piped()).stderr(Stdio::piped());

    let mut view_child = view
        .spawn()
        .map_err(|err| map_command_spawn_error(err, samtools))?;
    let view_stdout = view_child
        .stdout
        .take()
        .ok_or_else(|| AppError::ParseError {
            message: "failed to capture samtools view stdout for piping".to_string(),
        })?;

    let mut sort = Command::new(samtools);
    sort.arg("sort")
        .arg("-@")
        .arg(&threads_arg)
        .arg("-o")
        .arg(output_bam)
        .arg("-")
        .stdin(Stdio::from(view_stdout))
        .stderr(Stdio::piped());

    let sort_output = sort
        .output()
        .map_err(|err| map_command_spawn_error(err, samtools))?;
    let view_output = view_child.wait_with_output()?;

    if view_output.status.success() && sort_output.status.success() {
        return Ok(());
    }

    let view_stderr = String::from_utf8_lossy(&view_output.stderr);
    let sort_stderr = String::from_utf8_lossy(&sort_output.stderr);
    let mut command = format!(
        "{samtools} view -@ {threads_arg} -bh {}",
        input_bam.display()
    );
    if !fraction_is_one(fraction) {
        command.push_str(&format!(" -s {}", samtools_sample_arg(seed, fraction)));
    }
    for region in regions {
        command.push(' ');
        command.push_str(region);
    }
    command.push_str(&format!(
        " | {samtools} sort -@ {threads_arg} -o {} -",
        output_bam.display()
    ));

    Err(AppError::CommandFailed {
        command,
        code: sort_output.status.code().or(view_output.status.code()),
        stderr: format!(
            "samtools view stderr:\n{}\nsamtools sort stderr:\n{}",
            view_stderr.trim(),
            sort_stderr.trim()
        ),
    })
}

fn run_samtools_index(samtools: &str, threads: usize, bam_path: &Path) -> Result<()> {
    let mut command = Command::new(samtools);
    command
        .arg("index")
        .arg("-@")
        .arg(threads.to_string())
        .arg(bam_path);
    run_command_checked(
        command,
        samtools,
        format!("{samtools} index -@ {} {}", threads, bam_path.display()),
    )?;
    Ok(())
}

fn run_samtools_quickcheck(samtools: &str, bam_path: &Path) -> Result<()> {
    let mut command = Command::new(samtools);
    command.arg("quickcheck").arg(bam_path);
    run_command_checked(
        command,
        samtools,
        format!("{samtools} quickcheck {}", bam_path.display()),
    )?;
    Ok(())
}

fn write_samtools_report(
    samtools: &str,
    subcommand: &str,
    bam: &Path,
    output: &Path,
) -> Result<()> {
    let mut command = Command::new(samtools);
    command.arg(subcommand).arg(bam);
    let command_output = run_command_checked(
        command,
        samtools,
        format!("{samtools} {subcommand} {}", bam.display()),
    )?;
    fs::write(output, command_output.stdout)?;
    Ok(())
}

fn write_manifest(manifest_path: &Path, rows: &[ManifestRow]) -> Result<()> {
    let file = File::create(manifest_path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "role\tbam_path\tbai_path\tfraction\tseed\tsha256_bam\tsha256_bai"
    )?;

    for row in rows {
        let bam_hash = sha256sum(&row.bam_path)?;
        let bai_hash = sha256sum(&row.bai_path)?;
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.role,
            row.bam_path.display(),
            row.bai_path.display(),
            format_fraction(row.fraction),
            row.seed,
            bam_hash,
            bai_hash
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn sha256sum(path: &Path) -> Result<String> {
    let mut command = Command::new("sha256sum");
    command.arg(path);
    let output = run_command_checked(
        command,
        "sha256sum",
        format!("sha256sum {}", path.display()),
    )?;
    let text = String::from_utf8(output.stdout).map_err(|_| AppError::ParseError {
        message: format!("sha256sum produced non-utf8 output for {}", path.display()),
    })?;
    let hash = text
        .split_whitespace()
        .next()
        .ok_or_else(|| AppError::ParseError {
            message: format!("failed to parse sha256sum output for {}", path.display()),
        })?;
    Ok(hash.to_string())
}

fn run_command_checked(
    mut command: Command,
    command_name: &str,
    display: String,
) -> Result<Output> {
    let output = command
        .output()
        .map_err(|err| map_command_spawn_error(err, command_name))?;

    if output.status.success() {
        return Ok(output);
    }

    Err(AppError::CommandFailed {
        command: display,
        code: output.status.code(),
        stderr: String::from_utf8_lossy(&output.stderr).to_string(),
    })
}

fn map_command_spawn_error(error: std::io::Error, command_name: &str) -> AppError {
    if error.kind() == std::io::ErrorKind::NotFound {
        return AppError::CommandNotFound {
            command: command_name.to_string(),
        };
    }
    AppError::Io(error)
}

fn ensure_command_exists(command: &str) -> Result<()> {
    if command_exists(command) {
        return Ok(());
    }
    Err(AppError::CommandNotFound {
        command: command.to_string(),
    })
}

fn command_exists(command: &str) -> bool {
    if Path::new(command).is_file() {
        return true;
    }
    std::env::var_os("PATH").is_some_and(|paths| {
        std::env::split_paths(&paths)
            .map(|dir| dir.join(command))
            .any(|full_path| full_path.is_file())
    })
}

fn bai_path_for(bam_path: &Path) -> PathBuf {
    let mut bam_os: OsString = bam_path.as_os_str().to_owned();
    bam_os.push(OsStr::new(".bai"));
    PathBuf::from(bam_os)
}

fn fraction_is_one(fraction: f64) -> bool {
    (fraction - 1.0).abs() <= f64::EPSILON
}

fn samtools_sample_arg(seed: u64, fraction: f64) -> String {
    let suffix = format_fraction(fraction);
    if let Some(stripped) = suffix.strip_prefix("0.") {
        return format!("{seed}.{stripped}");
    }
    format!("{seed}.{suffix}")
}

fn format_fraction(value: f64) -> String {
    let mut text = format!("{value:.12}");
    while text.ends_with('0') {
        text.pop();
    }
    if text.ends_with('.') {
        text.pop();
    }
    text
}

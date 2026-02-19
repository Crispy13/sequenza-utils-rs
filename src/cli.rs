use crate::errors::{AppError, Result};
use clap::error::ErrorKind;
use clap::{ArgAction, Parser};

#[derive(Debug, Clone, Parser)]
#[command(disable_help_flag = true, disable_version_flag = true)]
struct CliArgs {
    #[arg(short = 'p', long = "pileup", action = ArgAction::SetTrue)]
    pileup: bool,
    #[arg(short = 'n', long = "normal")]
    normal: Option<String>,
    #[arg(short = 't', long = "tumor")]
    tumor: Option<String>,
    #[arg(long = "gc")]
    gc: Option<String>,
    #[arg(short = 'F', long = "fasta")]
    fasta: Option<String>,
    #[arg(short = 'o', long = "output", default_value = "-")]
    out: String,
    #[arg(long = "normal2")]
    normal2: Option<String>,
    #[arg(short = 'C', long = "chromosome", num_args = 1..)]
    chr: Vec<String>,
    #[arg(long = "parallel")]
    nproc: Option<String>,
    #[arg(long = "parallel-single-output", action = ArgAction::SetTrue)]
    parallel_single_output: bool,
    #[arg(short = 'S', long = "samtools", default_value = "samtools")]
    samtools: String,
    #[arg(short = 'T', long = "tabix", default_value = "tabix")]
    tabix: String,
    #[arg(short = 'q', long = "qlimit", default_value = "20")]
    qlimit: String,
    #[arg(short = 'f', long = "qformat", default_value = "sanger")]
    qformat: String,
    #[arg(short = 'N', default_value = "20")]
    depth_sum: String,
    #[arg(long = "hom", default_value = "0.9")]
    hom: String,
    #[arg(long = "het", default_value = "0.25")]
    het: String,
    #[arg(long = "het_f", default_value = "-0.2")]
    het_f: String,
    #[arg(long = "progress", action = ArgAction::SetTrue)]
    progress: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QualityFormat {
    Sanger,
    Illumina,
}

impl QualityFormat {
    pub fn parse(value: &str) -> Result<Self> {
        match value {
            "sanger" => Ok(Self::Sanger),
            "illumina" => Ok(Self::Illumina),
            other => Err(AppError::InvalidValue {
                flag: "--qformat".to_string(),
                value: other.to_string(),
                reason: "Supported quality format are only \"illumina\" and \"sanger\"(default)."
                    .to_string(),
            }),
        }
    }

    pub fn offset(&self) -> i32 {
        match self {
            Self::Sanger => 33,
            Self::Illumina => 64,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Bam2SeqzArgs {
    pub pileup: bool,
    pub normal: String,
    pub tumor: String,
    pub gc: String,
    pub fasta: Option<String>,
    pub out: String,
    pub normal2: Option<String>,
    pub chr: Vec<String>,
    pub nproc: usize,
    pub parallel_single_output: bool,
    pub samtools: String,
    pub tabix: String,
    pub qlimit: i32,
    pub qformat: QualityFormat,
    pub depth_sum: i32,
    pub hom: f64,
    pub het: f64,
    pub het_f: f64,
    pub progress: bool,
}

impl Default for Bam2SeqzArgs {
    fn default() -> Self {
        Self {
            pileup: false,
            normal: String::new(),
            tumor: String::new(),
            gc: String::new(),
            fasta: None,
            out: "-".to_string(),
            normal2: None,
            chr: Vec::new(),
            nproc: 1,
            parallel_single_output: false,
            samtools: "samtools".to_string(),
            tabix: "tabix".to_string(),
            qlimit: 20,
            qformat: QualityFormat::Sanger,
            depth_sum: 20,
            hom: 0.9,
            het: 0.25,
            het_f: -0.2,
            progress: false,
        }
    }
}

impl Bam2SeqzArgs {
    pub fn qlimit_ascii(&self) -> i32 {
        self.qlimit + self.qformat.offset()
    }

    pub fn has_explicit_ranged_regions(&self) -> bool {
        self.chr.iter().any(|region| region.contains(':'))
    }

    pub fn validate(&self) -> Result<()> {
        if self.normal.is_empty() {
            return Err(AppError::MissingRequired {
                field: "--normal".to_string(),
            });
        }
        if self.tumor.is_empty() {
            return Err(AppError::MissingRequired {
                field: "--tumor".to_string(),
            });
        }
        if self.gc.is_empty() {
            return Err(AppError::MissingRequired {
                field: "-gc".to_string(),
            });
        }
        if !self.pileup && self.fasta.is_none() {
            return Err(AppError::MissingRequired {
                field: "--fasta (required when input is BAM)".to_string(),
            });
        }
        if self.nproc > 1 && self.chr.len() < 2 {
            let can_auto_bin_without_ranges = !self.pileup && !self.has_explicit_ranged_regions();
            if !can_auto_bin_without_ranges {
                return Err(AppError::InvalidValue {
                    flag: "--parallel".to_string(),
                    value: self.nproc.to_string(),
                    reason: "--chromosome must include at least two regions".to_string(),
                });
            }
        }
        if self.nproc > 1 && self.out == "-" {
            return Err(AppError::InvalidValue {
                flag: "--output".to_string(),
                value: "-".to_string(),
                reason: "ERROR: The --output parameter must be different then STDOUT when used with the --parallel option".to_string(),
            });
        }
        Ok(())
    }
}

pub fn parse_from_env() -> Result<Bam2SeqzArgs> {
    parse_args(std::env::args())
}

pub fn parse_args<I, S>(args: I) -> Result<Bam2SeqzArgs>
where
    I: IntoIterator<Item = S>,
    S: Into<String>,
{
    let mut tokens: Vec<String> = args.into_iter().map(Into::into).collect();
    if tokens.is_empty() {
        tokens.push("bam2seqz_rs".to_string());
    }

    let normalized = normalize_legacy_tokens(tokens);
    let cli = CliArgs::try_parse_from(normalized).map_err(map_clap_error)?;

    let parsed = Bam2SeqzArgs {
        pileup: cli.pileup,
        normal: cli.normal.unwrap_or_default(),
        tumor: cli.tumor.unwrap_or_default(),
        gc: cli.gc.unwrap_or_default(),
        fasta: cli.fasta,
        out: cli.out,
        normal2: cli.normal2,
        chr: cli.chr,
        nproc: cli
            .nproc
            .as_deref()
            .map_or(Ok(1), |value| parse_usize("--parallel", value))?,
        parallel_single_output: cli.parallel_single_output,
        samtools: cli.samtools,
        tabix: cli.tabix,
        qlimit: parse_i32("-q", &cli.qlimit)?,
        qformat: QualityFormat::parse(&cli.qformat)?,
        depth_sum: parse_i32("-N", &cli.depth_sum)?,
        hom: parse_f64("--hom", &cli.hom)?,
        het: parse_f64("--het", &cli.het)?,
        het_f: parse_f64("--het_f", &cli.het_f)?,
        progress: cli.progress,
    };

    parsed.validate()?;
    Ok(parsed)
}

fn normalize_legacy_tokens(tokens: Vec<String>) -> Vec<String> {
    tokens
        .into_iter()
        .map(|token| {
            if token == "-gc" {
                "--gc".to_string()
            } else if token == "-n2" {
                "--normal2".to_string()
            } else if let Some(value) = token.strip_prefix("-gc=") {
                format!("--gc={value}")
            } else if let Some(value) = token.strip_prefix("-n2=") {
                format!("--normal2={value}")
            } else {
                token
            }
        })
        .collect()
}

fn map_clap_error(error: clap::Error) -> AppError {
    let kind = error.kind();
    let rendered = error.to_string();
    match kind {
        ErrorKind::UnknownArgument => AppError::UnsupportedArgument {
            arg: first_quoted_token(&rendered).unwrap_or(rendered),
        },
        ErrorKind::TooFewValues | ErrorKind::WrongNumberOfValues => AppError::MissingValue {
            flag: first_quoted_token(&rendered).unwrap_or_else(|| "argument".to_string()),
        },
        _ => AppError::ParseError {
            message: clap_error_message(&rendered),
        },
    }
}

fn first_quoted_token(message: &str) -> Option<String> {
    let start = message.find('\'')?;
    let end = message[start + 1..].find('\'')?;
    Some(message[start + 1..start + 1 + end].to_string())
}

fn clap_error_message(message: &str) -> String {
    message
        .lines()
        .find_map(|line| line.strip_prefix("error: "))
        .or_else(|| message.lines().next())
        .unwrap_or("failed to parse arguments")
        .to_string()
}

fn parse_usize(flag: &str, value: &str) -> Result<usize> {
    value.parse::<usize>().map_err(|_| AppError::InvalidValue {
        flag: flag.to_string(),
        value: value.to_string(),
        reason: "must be a positive integer".to_string(),
    })
}

fn parse_i32(flag: &str, value: &str) -> Result<i32> {
    value.parse::<i32>().map_err(|_| AppError::InvalidValue {
        flag: flag.to_string(),
        value: value.to_string(),
        reason: "must be an integer".to_string(),
    })
}

fn parse_f64(flag: &str, value: &str) -> Result<f64> {
    value.parse::<f64>().map_err(|_| AppError::InvalidValue {
        flag: flag.to_string(),
        value: value.to_string(),
        reason: "must be a floating-point number".to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::{QualityFormat, parse_args};

    #[test]
    fn parses_minimal_bam_arguments() {
        let args = parse_args([
            "bam2seqz_rs",
            "-n",
            "normal.bam",
            "-t",
            "tumor.bam",
            "-gc",
            "gc.wig.gz",
            "-F",
            "ref.fa",
        ])
        .expect("expected parse success");

        assert!(!args.pileup);
        assert_eq!(args.normal, "normal.bam");
        assert_eq!(args.tumor, "tumor.bam");
        assert_eq!(args.qformat, QualityFormat::Sanger);
        assert_eq!(args.qlimit_ascii(), 53);
    }

    #[test]
    fn parses_pileup_without_fasta() {
        let args = parse_args([
            "bam2seqz_rs",
            "--pileup",
            "-n",
            "normal.pileup.gz",
            "-t",
            "tumor.pileup.gz",
            "-gc",
            "gc.wig.gz",
        ])
        .expect("expected parse success");

        assert!(args.pileup);
        assert!(args.fasta.is_none());
    }

    #[test]
    fn rejects_missing_fasta_for_bam_input() {
        let result = parse_args(["bam2seqz_rs", "-n", "n.bam", "-t", "t.bam", "-gc", "gc.wig"]);
        assert!(result.is_err());
    }

    #[test]
    fn rejects_parallel_stdout() {
        let result = parse_args([
            "bam2seqz_rs",
            "-n",
            "n.bam",
            "-t",
            "t.bam",
            "-gc",
            "gc.wig",
            "-F",
            "ref.fa",
            "-C",
            "7",
            "12",
            "--parallel",
            "2",
            "-o",
            "-",
        ]);
        assert!(result.is_err());
    }

    #[test]
    fn parses_quality_format_illumina() {
        let args = parse_args([
            "bam2seqz_rs",
            "-n",
            "normal.bam",
            "-t",
            "tumor.bam",
            "-gc",
            "gc.wig.gz",
            "-F",
            "ref.fa",
            "-f",
            "illumina",
        ])
        .expect("expected parse success");

        assert_eq!(args.qformat, QualityFormat::Illumina);
        assert_eq!(args.qlimit_ascii(), 84);
    }

    #[test]
    fn parses_legacy_normal2_flag() {
        let args = parse_args([
            "bam2seqz_rs",
            "-n",
            "normal.bam",
            "-t",
            "tumor.bam",
            "-gc",
            "gc.wig.gz",
            "-F",
            "ref.fa",
            "-n2",
            "normal2.bam",
        ])
        .expect("expected parse success");

        assert_eq!(args.normal2.as_deref(), Some("normal2.bam"));
    }

    #[test]
    fn parses_progress_flag() {
        let args = parse_args([
            "bam2seqz_rs",
            "-n",
            "normal.bam",
            "-t",
            "tumor.bam",
            "-gc",
            "gc.wig.gz",
            "-F",
            "ref.fa",
            "--progress",
        ])
        .expect("expected parse success");

        assert!(args.progress);
    }

    #[test]
    fn parses_parallel_single_output_flag() {
        let args = parse_args([
            "bam2seqz_rs",
            "-n",
            "normal.bam",
            "-t",
            "tumor.bam",
            "-gc",
            "gc.wig.gz",
            "-F",
            "ref.fa",
            "-C",
            "chr20:1-100",
            "chr20:101-200",
            "--parallel",
            "2",
            "--parallel-single-output",
            "-o",
            "out.seqz.gz",
        ])
        .expect("expected parse success");

        assert!(args.parallel_single_output);
    }

    #[test]
    fn accepts_parallel_without_chromosome_for_bam_input() {
        let args = parse_args([
            "bam2seqz_rs",
            "-n",
            "normal.bam",
            "-t",
            "tumor.bam",
            "-gc",
            "gc.wig.gz",
            "-F",
            "ref.fa",
            "--parallel",
            "4",
            "-o",
            "out.seqz.gz",
        ])
        .expect("expected parse success");

        assert_eq!(args.nproc, 4);
        assert!(args.chr.is_empty());
    }

    #[test]
    fn rejects_single_ranged_region_with_parallel() {
        let result = parse_args([
            "bam2seqz_rs",
            "-n",
            "normal.bam",
            "-t",
            "tumor.bam",
            "-gc",
            "gc.wig.gz",
            "-F",
            "ref.fa",
            "-C",
            "chr20:1-100",
            "--parallel",
            "2",
            "-o",
            "out.seqz.gz",
        ]);

        assert!(result.is_err());
    }
}

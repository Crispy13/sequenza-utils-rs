use std::fs;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;

fn workspace_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn golden_dir() -> PathBuf {
    workspace_dir()
        .join("copilot-desk")
        .join("golden")
        .join("python")
}

fn run_binary(args: &[&str]) -> std::process::Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_bam2seqz_rs"));
    command.args(args).current_dir(workspace_dir());
    command
        .output()
        .expect("expected bam2seqz_rs binary to execute")
}

const NORMAL_BAM: &str = "testdata/NA12878.chrom20.synthetic.seed20260218.normal.chr.bam";
const TUMOR_BAM: &str = "testdata/NA12878.chrom20.synthetic.seed20260218.tumor.chr.bam";
const NORMAL_PILEUP: &str =
    "testdata/NA12878.chrom20.synthetic.seed20260218.normal.region59900_64000.pileup.gz";
const TUMOR_PILEUP: &str =
    "testdata/NA12878.chrom20.synthetic.seed20260218.tumor.region59900_64000.pileup.gz";
const GC_WIG: &str = "testdata/hg19_gc50.wig.gz";
const FASTA: &str = "testdata/hg19.fa";
const REGION_A: &str = "chr20:59900-62000";
const REGION_B: &str = "chr20:62001-64000";

fn tabix_index_path(output_gz: &Path) -> PathBuf {
    let mut path = output_gz.as_os_str().to_os_string();
    path.push(".tbi");
    PathBuf::from(path)
}

fn assert_output_and_index_exist(output_gz: &Path) {
    assert!(
        output_gz.is_file(),
        "expected output artifact to exist: {}",
        output_gz.display()
    );
    let index = tabix_index_path(output_gz);
    assert!(
        index.is_file(),
        "expected tabix index artifact to exist: {}",
        index.display()
    );
}

fn assert_stderr_contains_golden_final_semantic_line(actual_stderr: &str, golden_stderr: &str) {
    let expected_final = golden_stderr
        .lines()
        .rev()
        .find(|line| !line.trim().is_empty())
        .expect("expected non-empty final stderr line in golden artifact");
    let semantic_tail = expected_final
        .strip_prefix("Exception: ")
        .unwrap_or(expected_final);
    assert!(
        actual_stderr.contains(expected_final) || actual_stderr.contains(semantic_tail),
        "stderr did not contain expected final semantic line: {expected_final}"
    );
}

fn gunzip_to_text(gz_path: &Path) -> String {
    let output = Command::new("gzip")
        .arg("-cd")
        .arg(gz_path)
        .output()
        .expect("expected gzip command to run");
    assert!(
        output.status.success(),
        "gzip failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).expect("expected utf8 output")
}

#[test]
fn golden_directory_exists() {
    assert!(golden_dir().is_dir());
}

#[test]
fn golden_manifest_and_seqz_files_present() {
    let dir = golden_dir();
    let manifest = dir.join("SHA256SUMS.txt");
    assert!(manifest.is_file());

    let entries = fs::read_dir(dir)
        .expect("expected golden directory to be readable")
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.path())
        .collect::<Vec<_>>();

    let has_seqz = entries.iter().any(|path| {
        path.extension()
            .and_then(|ext| ext.to_str())
            .map(|ext| ext == "seqz")
            .unwrap_or(false)
    });

    assert!(has_seqz, "expected at least one .seqz golden file");
}

#[test]
fn rust_matches_golden_bam_default() {
    let out = workspace_dir()
        .join("target")
        .join("it_bam_default.seqz.gz");
    let _ = fs::remove_file(&out);
    let _ = fs::remove_file(tabix_index_path(&out));
    let output = run_binary(&[
        "-n",
        NORMAL_BAM,
        "-t",
        TUMOR_BAM,
        "-gc",
        GC_WIG,
        "-F",
        FASTA,
        "-o",
        "target/it_bam_default.seqz.gz",
    ]);
    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_output_and_index_exist(&out);

    let actual = gunzip_to_text(&out);
    let expected =
        fs::read_to_string(golden_dir().join("bam_default.seqz")).expect("expected golden text");
    assert_eq!(actual, expected);
}

#[test]
fn rust_matches_golden_pileup_default() {
    let out = workspace_dir()
        .join("target")
        .join("it_pileup_default.seqz.gz");
    let _ = fs::remove_file(&out);
    let _ = fs::remove_file(tabix_index_path(&out));
    let output = run_binary(&[
        "--pileup",
        "-n",
        NORMAL_PILEUP,
        "-t",
        TUMOR_PILEUP,
        "-gc",
        GC_WIG,
        "-o",
        "target/it_pileup_default.seqz.gz",
    ]);
    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_output_and_index_exist(&out);

    let actual = gunzip_to_text(&out);
    let expected =
        fs::read_to_string(golden_dir().join("pileup_default.seqz")).expect("expected golden text");
    assert_eq!(actual, expected);
}

#[test]
fn rust_matches_golden_region_12() {
    let region_out = workspace_dir().join("target").join("it_region_12.seqz.gz");
    let _ = fs::remove_file(&region_out);
    let _ = fs::remove_file(tabix_index_path(&region_out));
    let region = run_binary(&[
        "-n",
        NORMAL_BAM,
        "-t",
        TUMOR_BAM,
        "-gc",
        GC_WIG,
        "-F",
        FASTA,
        "-C",
        REGION_B,
        "-o",
        "target/it_region_12.seqz.gz",
    ]);
    assert!(region.status.success(), "region run failed");
    assert_output_and_index_exist(&region_out);
    assert_eq!(
        gunzip_to_text(&region_out),
        fs::read_to_string(golden_dir().join("bam_region_12.seqz")).expect("golden region")
    );
}

#[test]
fn rust_matches_golden_parallel_regions() {
    let parallel_prefix = workspace_dir().join("target").join("it_parallel.seqz.gz");
    let _ = fs::remove_file(&parallel_prefix);
    let _ = fs::remove_file(
        workspace_dir()
            .join("target")
            .join("it_parallel_chr20_59900_62000.seqz.gz"),
    );
    let _ = fs::remove_file(
        workspace_dir()
            .join("target")
            .join("it_parallel_chr20_62001_64000.seqz.gz"),
    );
    let _ = fs::remove_file(
        workspace_dir()
            .join("target")
            .join("it_parallel_chr20_59900_62000.seqz.gz.tbi"),
    );
    let _ = fs::remove_file(
        workspace_dir()
            .join("target")
            .join("it_parallel_chr20_62001_64000.seqz.gz.tbi"),
    );
    let parallel = run_binary(&[
        "-n",
        NORMAL_BAM,
        "-t",
        TUMOR_BAM,
        "-gc",
        GC_WIG,
        "-F",
        FASTA,
        "-C",
        REGION_A,
        REGION_B,
        "--parallel",
        "2",
        "-o",
        "target/it_parallel.seqz.gz",
    ]);
    assert!(parallel.status.success(), "parallel run failed");

    let parallel_7 = workspace_dir()
        .join("target")
        .join("it_parallel_chr20_59900_62000.seqz.gz");
    let parallel_12 = workspace_dir()
        .join("target")
        .join("it_parallel_chr20_62001_64000.seqz.gz");
    assert_output_and_index_exist(&parallel_7);
    assert_output_and_index_exist(&parallel_12);

    let actual7 = gunzip_to_text(&parallel_7);
    let actual12 = gunzip_to_text(&parallel_12);
    let expected7 =
        fs::read_to_string(golden_dir().join("bam_parallel_7.seqz")).expect("golden parallel 7");
    let expected12 =
        fs::read_to_string(golden_dir().join("bam_parallel_12.seqz")).expect("golden parallel 12");
    assert_eq!(actual7, expected7);
    assert_eq!(actual12, expected12);
}

#[test]
fn rust_parallel_single_output_preserves_region_order() {
    let out = workspace_dir().join("target").join("it_parallel_single.seqz.gz");
    let _ = fs::remove_file(&out);
    let _ = fs::remove_file(tabix_index_path(&out));

    let run = run_binary(&[
        "-n",
        NORMAL_BAM,
        "-t",
        TUMOR_BAM,
        "-gc",
        GC_WIG,
        "-F",
        FASTA,
        "-C",
        REGION_A,
        REGION_B,
        "--parallel",
        "2",
        "--parallel-single-output",
        "-o",
        "target/it_parallel_single.seqz.gz",
    ]);
    assert!(
        run.status.success(),
        "parallel single-output run failed: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    assert_output_and_index_exist(&out);

    let actual = gunzip_to_text(&out);
    let expected_a =
        fs::read_to_string(golden_dir().join("bam_parallel_7.seqz")).expect("golden parallel 7");
    let expected_b = fs::read_to_string(golden_dir().join("bam_parallel_12.seqz"))
        .expect("golden parallel 12");
    let expected_b_body = expected_b
        .split_once('\n')
        .map(|(_, body)| body)
        .unwrap_or("");
    let expected = format!("{expected_a}{expected_b_body}");

    assert_eq!(actual, expected);
}

#[test]
fn rust_matches_golden_normal2() {
    let normal2_out = workspace_dir().join("target").join("it_normal2.seqz.gz");
    let _ = fs::remove_file(&normal2_out);
    let _ = fs::remove_file(tabix_index_path(&normal2_out));
    let normal2 = run_binary(&[
        "-n",
        NORMAL_BAM,
        "-t",
        TUMOR_BAM,
        "-gc",
        GC_WIG,
        "-F",
        FASTA,
        "--normal2",
        NORMAL_BAM,
        "-o",
        "target/it_normal2.seqz.gz",
    ]);
    assert!(normal2.status.success(), "normal2 run failed");
    assert_output_and_index_exist(&normal2_out);
    assert_eq!(
        gunzip_to_text(&normal2_out),
        fs::read_to_string(golden_dir().join("bam_normal2.seqz")).expect("golden normal2")
    );
}

#[test]
fn rust_matches_error_bad_qformat_semantics() {
    let bad_q = run_binary(&[
        "-n",
        NORMAL_BAM,
        "-t",
        TUMOR_BAM,
        "-gc",
        GC_WIG,
        "-F",
        FASTA,
        "-f",
        "invalid",
        "-o",
        "target/it_bad_q.seqz.gz",
    ]);
    assert!(!bad_q.status.success());

    let expected_bad_exit = fs::read_to_string(golden_dir().join("error_bad_qformat.exit"))
        .expect("expected golden exit")
        .trim()
        .to_string();
    assert_eq!(
        bad_q.status.code().unwrap_or(1).to_string(),
        expected_bad_exit
    );

    let expected_bad_stdout =
        fs::read_to_string(golden_dir().join("error_bad_qformat.stdout")).expect("golden stdout");
    let actual_bad_stdout = String::from_utf8_lossy(&bad_q.stdout);
    assert_eq!(actual_bad_stdout, expected_bad_stdout);

    let expected_bad_stderr =
        fs::read_to_string(golden_dir().join("error_bad_qformat.stderr")).expect("golden stderr");
    let actual_bad_stderr = String::from_utf8_lossy(&bad_q.stderr);
    assert_stderr_contains_golden_final_semantic_line(&actual_bad_stderr, &expected_bad_stderr);
}

#[test]
fn rust_matches_error_parallel_stdout_semantics() {
    let parallel_stdout = run_binary(&[
        "-n",
        NORMAL_BAM,
        "-t",
        TUMOR_BAM,
        "-gc",
        GC_WIG,
        "-F",
        FASTA,
        "-C",
        REGION_A,
        REGION_B,
        "--parallel",
        "2",
        "-o",
        "-",
    ]);
    assert!(!parallel_stdout.status.success());

    let expected_parallel_exit =
        fs::read_to_string(golden_dir().join("error_parallel_stdout.exit"))
            .expect("expected golden exit")
            .trim()
            .to_string();
    assert_eq!(
        parallel_stdout.status.code().unwrap_or(1).to_string(),
        expected_parallel_exit
    );

    let expected_parallel_stdout =
        fs::read_to_string(golden_dir().join("error_parallel_stdout.stdout"))
            .expect("golden stdout");
    let actual_parallel_stdout = String::from_utf8_lossy(&parallel_stdout.stdout);
    assert_eq!(actual_parallel_stdout, expected_parallel_stdout);

    let expected_parallel_stderr =
        fs::read_to_string(golden_dir().join("error_parallel_stdout.stderr"))
            .expect("golden stderr");
    let actual_parallel_stderr = String::from_utf8_lossy(&parallel_stdout.stderr);
    assert_stderr_contains_golden_final_semantic_line(
        &actual_parallel_stderr,
        &expected_parallel_stderr,
    );
}

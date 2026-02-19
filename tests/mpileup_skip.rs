use sequenza_utils::seqz_core::seqz_header;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const NORMAL_BAM: &str = "testdata/NA12878.chrom20.synthetic.seed20260218.normal.chr.bam";
const TUMOR_BAM: &str = "testdata/NA12878.chrom20.synthetic.seed20260218.tumor.chr.bam";
const GC_WIG: &str = "testdata/hg19_gc50.wig.gz";
const FASTA: &str = "testdata/hg19.fa";

const MISSING_SAMTOOLS: &str = "missing_samtools_for_mpileup_skip_test";
const MISSING_TABIX: &str = "missing_tabix_for_mpileup_skip_test";

fn workspace_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn run_binary(args: &[&str]) -> std::process::Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_bam2seqz"));
    command.args(args).current_dir(workspace_dir());
    command
        .output()
        .expect("expected bam2seqz binary to execute")
}

fn expected_header() -> String {
    format!("{}\n", seqz_header().join("\t"))
}

fn assert_header_only(path: &Path) {
    assert!(path.is_file(), "expected output file to exist: {}", path.display());
    let content = fs::read_to_string(path).expect("expected output content");
    assert_eq!(content, expected_header());
}

#[test]
fn single_thread_region_without_gc_overlap_skips_mpileup() {
    let out = workspace_dir().join("target").join("it_skip_single.seqz");
    let _ = fs::remove_file(&out);

    let output = run_binary(&[
        "-n",
        NORMAL_BAM,
        "-t",
        TUMOR_BAM,
        "-gc",
        GC_WIG,
        "-F",
        FASTA,
        "-C",
        "chrNoSuch:1-10",
        "-S",
        MISSING_SAMTOOLS,
        "-T",
        MISSING_TABIX,
        "-o",
        "target/it_skip_single.seqz",
    ]);

    assert!(
        output.status.success(),
        "expected skip path success without samtools: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_header_only(&out);
}

#[test]
fn parallel_multi_output_regions_without_gc_overlap_skip_mpileup() {
    let out_a = workspace_dir()
        .join("target")
        .join("it_skip_parallel_chrNoSuch_1_10.seqz");
    let out_b = workspace_dir()
        .join("target")
        .join("it_skip_parallel_chrNoSuch_11_20.seqz");
    let _ = fs::remove_file(&out_a);
    let _ = fs::remove_file(&out_b);

    let output = run_binary(&[
        "-n",
        NORMAL_BAM,
        "-t",
        TUMOR_BAM,
        "-gc",
        GC_WIG,
        "-F",
        FASTA,
        "-C",
        "chrNoSuch:1-10",
        "chrNoSuch:11-20",
        "--parallel",
        "2",
        "-S",
        MISSING_SAMTOOLS,
        "-T",
        MISSING_TABIX,
        "-o",
        "target/it_skip_parallel.seqz",
    ]);

    assert!(
        output.status.success(),
        "expected skip path success without samtools: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_header_only(&out_a);
    assert_header_only(&out_b);
}

#[test]
fn parallel_single_output_regions_without_gc_overlap_skip_mpileup() {
    let out = workspace_dir()
        .join("target")
        .join("it_skip_parallel_single.seqz");
    let _ = fs::remove_file(&out);

    let output = run_binary(&[
        "-n",
        NORMAL_BAM,
        "-t",
        TUMOR_BAM,
        "-gc",
        GC_WIG,
        "-F",
        FASTA,
        "-C",
        "chrNoSuch:1-10",
        "chrNoSuch:11-20",
        "--parallel",
        "2",
        "--parallel-single-output",
        "-S",
        MISSING_SAMTOOLS,
        "-T",
        MISSING_TABIX,
        "-o",
        "target/it_skip_parallel_single.seqz",
    ]);

    assert!(
        output.status.success(),
        "expected skip path success without samtools: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_header_only(&out);
}

#[test]
fn normal2_region_without_gc_overlap_skips_all_mpileup_streams() {
    let out = workspace_dir().join("target").join("it_skip_normal2.seqz");
    let _ = fs::remove_file(&out);

    let output = run_binary(&[
        "-n",
        NORMAL_BAM,
        "-t",
        TUMOR_BAM,
        "--normal2",
        NORMAL_BAM,
        "-gc",
        GC_WIG,
        "-F",
        FASTA,
        "-C",
        "chrNoSuch:1-10",
        "-S",
        MISSING_SAMTOOLS,
        "-T",
        MISSING_TABIX,
        "-o",
        "target/it_skip_normal2.seqz",
    ]);

    assert!(
        output.status.success(),
        "expected skip path success without samtools: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_header_only(&out);
}

#[test]
fn parallel_auto_binning_contig_without_gc_overlap_skips_mpileup() {
    let out = workspace_dir().join("target").join("it_skip_auto_bin.seqz");
    let _ = fs::remove_file(&out);

    let output = run_binary(&[
        "-n",
        NORMAL_BAM,
        "-t",
        TUMOR_BAM,
        "-gc",
        GC_WIG,
        "-F",
        FASTA,
        "-C",
        "chrNoSuch",
        "--parallel",
        "2",
        "-S",
        MISSING_SAMTOOLS,
        "-T",
        MISSING_TABIX,
        "-o",
        "target/it_skip_auto_bin.seqz",
    ]);

    assert!(
        output.status.success(),
        "expected auto-bin skip path success without samtools: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_header_only(&out);
}

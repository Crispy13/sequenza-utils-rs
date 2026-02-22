use std::fs;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;

const NORMAL_BAM: &str = "testdata/NA12878.chrom20.synthetic.seed20260218.normal.chr.bam";
const TUMOR_BAM: &str = "testdata/NA12878.chrom20.synthetic.seed20260218.tumor.chr.bam";
const GC_WIG: &str = "testdata/hg19_gc50.wig.gz";
const FASTA: &str = "testdata/hg19.fa";
const REGIONS: [&str; 8] = [
    "chr20:59900-60412",
    "chr20:60413-60925",
    "chr20:60926-61438",
    "chr20:61439-61951",
    "chr20:61952-62464",
    "chr20:62465-62977",
    "chr20:62978-63490",
    "chr20:63491-64000",
];

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

fn region_output_path(base_output: &str, region: &str) -> PathBuf {
    let region_token = region.replace([':', '-'], "_");
    let stem = base_output
        .strip_suffix(".seqz.gz")
        .expect("expected .seqz.gz suffix");
    workspace_dir().join(format!("{stem}_{region_token}.seqz.gz"))
}

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

fn remove_outputs(base_output: &str) {
    for region in REGIONS {
        let output = region_output_path(base_output, region);
        let _ = fs::remove_file(&output);
        let _ = fs::remove_file(tabix_index_path(&output));
    }
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
#[ignore = "long-running; run manually with --ignored"]
fn long_running_parallel_regions_p8() {
    let output = "target/it_long_parallel_p8.seqz.gz";
    remove_outputs(output);

    let result = run_binary(&[
        "-n",
        NORMAL_BAM,
        "-t",
        TUMOR_BAM,
        "-gc",
        GC_WIG,
        "-F",
        FASTA,
        "-C",
        REGIONS[0],
        REGIONS[1],
        REGIONS[2],
        REGIONS[3],
        REGIONS[4],
        REGIONS[5],
        REGIONS[6],
        REGIONS[7],
        "--parallel",
        "8",
        "-o",
        output,
    ]);
    assert!(
        result.status.success(),
        "parallel p8 run failed: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    for region in REGIONS {
        let region_output = region_output_path(output, region);
        assert_output_and_index_exist(&region_output);
    }
}

#[test]
#[ignore = "long-running; run manually with --ignored"]
fn long_running_parallel_normal2_p8() {
    let output = "target/it_long_parallel_normal2_p8.seqz.gz";
    remove_outputs(output);

    let result = run_binary(&[
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
        "-C",
        REGIONS[0],
        REGIONS[1],
        REGIONS[2],
        REGIONS[3],
        REGIONS[4],
        REGIONS[5],
        REGIONS[6],
        REGIONS[7],
        "--parallel",
        "8",
        "-o",
        output,
    ]);
    assert!(
        result.status.success(),
        "parallel normal2 p8 run failed: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    for region in REGIONS {
        let region_output = region_output_path(output, region);
        assert_output_and_index_exist(&region_output);
    }
}

#[test]
#[ignore = "long-running; run manually with --ignored"]
fn long_running_parallel_p8_repeatability() {
    let output_a = "target/it_long_parallel_p8_a.seqz.gz";
    let output_b = "target/it_long_parallel_p8_b.seqz.gz";
    remove_outputs(output_a);
    remove_outputs(output_b);

    let args_prefix = [
        "-n",
        NORMAL_BAM,
        "-t",
        TUMOR_BAM,
        "-gc",
        GC_WIG,
        "-F",
        FASTA,
        "-C",
        REGIONS[0],
        REGIONS[1],
        REGIONS[2],
        REGIONS[3],
        REGIONS[4],
        REGIONS[5],
        REGIONS[6],
        REGIONS[7],
        "--parallel",
        "8",
    ];

    let mut args_a = args_prefix.to_vec();
    args_a.extend(["-o", output_a]);
    let run_a = run_binary(&args_a);
    assert!(
        run_a.status.success(),
        "repeatability run A failed: {}",
        String::from_utf8_lossy(&run_a.stderr)
    );

    let mut args_b = args_prefix.to_vec();
    args_b.extend(["-o", output_b]);
    let run_b = run_binary(&args_b);
    assert!(
        run_b.status.success(),
        "repeatability run B failed: {}",
        String::from_utf8_lossy(&run_b.stderr)
    );

    for region in REGIONS {
        let file_a = region_output_path(output_a, region);
        let file_b = region_output_path(output_b, region);
        assert_output_and_index_exist(&file_a);
        assert_output_and_index_exist(&file_b);
        assert_eq!(gunzip_to_text(&file_a), gunzip_to_text(&file_b));
    }
}

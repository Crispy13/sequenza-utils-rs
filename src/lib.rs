pub mod cli;
pub mod errors;
pub mod external_tools;
pub mod pipeline;
pub mod seqz_core;
pub mod writer;

use cli::Bam2SeqzArgs;
use errors::Result;
use std::sync::Once;
use tracing_subscriber::EnvFilter;

static TRACING_INIT: Once = Once::new();

pub fn init_tracing() {
    TRACING_INIT.call_once(|| {
        let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
        let _ = tracing_subscriber::fmt()
            .with_env_filter(filter)
            .with_target(false)
            .try_init();
    });
}

pub fn run_from_args(args: Bam2SeqzArgs) -> Result<()> {
    pipeline::run(&args)
}

#[cfg(test)]
mod tests {
    use crate::cli::{QualityFormat, parse_args};
    use std::path::PathBuf;

    #[test]
    fn wiring_parses_and_runs_stub_pipeline() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let normal = root
            .join("testdata/NA12878.chrom20.synthetic.seed20260218.normal.region59900_64000.pileup.gz")
            .to_string_lossy()
            .to_string();
        let tumor = root
            .join(
                "testdata/NA12878.chrom20.synthetic.seed20260218.tumor.region59900_64000.pileup.gz",
            )
            .to_string_lossy()
            .to_string();
        let gc = root
            .join("testdata/hg19_gc50.wig.gz")
            .to_string_lossy()
            .to_string();
        let out = root
            .join("target/test_wiring.seqz")
            .to_string_lossy()
            .to_string();

        let args = parse_args([
            "bam2seqz_rs",
            "--pileup",
            "-n",
            &normal,
            "-t",
            &tumor,
            "-gc",
            &gc,
            "-o",
            &out,
            "-f",
            "sanger",
        ])
        .expect("expected valid args");
        assert_eq!(args.qformat, QualityFormat::Sanger);
        assert!(crate::run_from_args(args).is_ok());
    }
}

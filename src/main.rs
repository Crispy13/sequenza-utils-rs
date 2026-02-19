use std::process::ExitCode;

fn main() -> ExitCode {
    match bam2seqz_rs::cli::parse_from_env().and_then(bam2seqz_rs::run_from_args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("bam2seqz_rs: {error}");
            ExitCode::from(1)
        }
    }
}

use std::process::ExitCode;
use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() -> ExitCode {
    sequenza_utils::init_tracing();
    match sequenza_utils::cli::parse_from_env().and_then(sequenza_utils::run_from_args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("bam2seqz: {error}");
            ExitCode::from(1)
        }
    }
}

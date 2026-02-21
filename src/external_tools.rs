use crate::cli::Bam2SeqzArgs;
use crate::errors::{AppError, Result};
use flate2::read::GzDecoder;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::process::{Child, ChildStdout, Command, Stdio};
use tempfile::{Builder, NamedTempFile};
use tracing::info;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExternalTools {
    pub samtools: String,
    pub tabix: String,
}

#[derive(Debug)]
pub struct CommandStream {
    pub child: Child,
    pub stdout: ChildStdout,
    pub stderr_capture: NamedTempFile,
    pub command: String,
}

impl CommandStream {
    pub fn read_stderr(&self) -> String {
        let mut bytes = Vec::new();
        if let Ok(mut file) = self.stderr_capture.reopen() {
            let _ = file.read_to_end(&mut bytes);
        }
        String::from_utf8_lossy(&bytes).to_string()
    }
}

impl ExternalTools {
    pub fn from_args(args: &Bam2SeqzArgs) -> Self {
        Self {
            samtools: args.samtools.clone(),
            tabix: args.tabix.clone(),
        }
    }

    pub fn bam_mpileup_command(&self, bam: &str, fasta: &str, regions: &[String]) -> Vec<String> {
        let mut command = vec![
            self.samtools.clone(),
            "view".to_string(),
            "-u".to_string(),
            bam.to_string(),
        ];
        command.extend(regions.iter().cloned());
        command.push("|".to_string());
        command.push(self.samtools.clone());
        command.push("mpileup".to_string());
        command.push("-f".to_string());
        command.push(fasta.to_string());
        command.push("-q".to_string());
        command.push("20".to_string());
        command.push("-Q".to_string());
        command.push("20".to_string());
        command.push("-".to_string());
        command
    }

    pub fn indexed_pileup_command(&self, pileup: &str, regions: &[String]) -> Vec<String> {
        let mut command = vec![self.tabix.clone(), pileup.to_string()];
        command.extend(regions.iter().cloned());
        command
    }

    pub fn bgzip(&self) -> String {
        if self.tabix.ends_with("tabix") {
            self.tabix.trim_end_matches("tabix").to_string() + "bgzip"
        } else {
            "bgzip".to_string()
        }
    }

    pub fn list_bam_chromosomes(&self, bam: &str) -> Result<Vec<String>> {
        let mut command = Command::new(&self.samtools);
        command.arg("idxstats").arg(bam);
        let lines = self.capture_stdout_lines(command, &self.samtools)?;

        let mut chromosomes = Vec::new();
        for line in lines {
            let parts = line.split('\t').collect::<Vec<_>>();
            if parts.len() < 4 {
                continue;
            }
            let chrom = parts[0];
            if chrom == "*" {
                continue;
            }
            let mapped = parts[2].parse::<u64>().unwrap_or(0);
            let unmapped = parts[3].parse::<u64>().unwrap_or(0);
            if (mapped > 0 || unmapped > 0) && !chromosomes.iter().any(|c| c == chrom) {
                chromosomes.push(chrom.to_string());
            }
        }
        Ok(chromosomes)
    }

    pub fn read_text_lines(&self, path: &str) -> Result<Vec<String>> {
        if path.ends_with(".gz") {
            let file = File::open(path)?;
            let decoder = GzDecoder::new(file);
            Self::collect_non_empty_lines(BufReader::new(decoder))
        } else {
            let file = File::open(path)?;
            Self::collect_non_empty_lines(BufReader::new(file))
        }
    }

    pub fn run_tabix_pileup(&self, pileup: &str, regions: &[String]) -> Result<Vec<String>> {
        if regions.is_empty() {
            return self.read_text_lines(pileup);
        }
        let mut command = Command::new(&self.tabix);
        command.arg(pileup);
        for region in regions {
            command.arg(region);
        }
        self.capture_stdout_lines(command, &self.tabix)
    }

    pub fn run_samtools_mpileup(
        &self,
        bam: &str,
        fasta: &str,
        regions: &[String],
    ) -> Result<Vec<String>> {
        let mut collected = Vec::new();
        if regions.is_empty() {
            let mut command = Command::new(&self.samtools);
            command
                .arg("mpileup")
                .arg("-f")
                .arg(fasta)
                .arg("-q")
                .arg("20")
                .arg("-Q")
                .arg("20")
                .arg(bam);
            collected.extend(self.capture_stdout_lines(command, &self.samtools)?);
            return Ok(collected);
        }

        for region in regions {
            let mut command = Command::new(&self.samtools);
            command
                .arg("mpileup")
                .arg("-f")
                .arg(fasta)
                .arg("-q")
                .arg("20")
                .arg("-Q")
                .arg("20")
                .arg("-r")
                .arg(region)
                .arg(bam);
            collected.extend(self.capture_stdout_lines(command, &self.samtools)?);
        }
        Ok(collected)
    }

    pub fn run_tabix_pileup_to_tempfile(
        &self,
        pileup: &str,
        regions: &[String],
    ) -> Result<NamedTempFile> {
        let mut command = Command::new(&self.tabix);
        command.arg(pileup);
        for region in regions {
            command.arg(region);
        }
        self.capture_stdout_to_tempfile(command, &self.tabix)
    }

    pub fn run_samtools_mpileup_to_tempfile(
        &self,
        bam: &str,
        fasta: &str,
        regions: &[String],
    ) -> Result<NamedTempFile> {
        if regions.is_empty() {
            let mut command = Command::new(&self.samtools);
            command
                .arg("mpileup")
                .arg("-f")
                .arg(fasta)
                .arg("-q")
                .arg("20")
                .arg("-Q")
                .arg("20")
                .arg(bam);
            return self.capture_stdout_to_tempfile(command, &self.samtools);
        }

        let temp = Builder::new()
            .prefix("bam2seqz_cmd_")
            .suffix(".txt")
            .tempfile_in("tmp")?;

        for region in regions {
            let append_file = OpenOptions::new().append(true).open(temp.path())?;
            let mut command = Command::new(&self.samtools);
            command
                .arg("mpileup")
                .arg("-f")
                .arg(fasta)
                .arg("-q")
                .arg("20")
                .arg("-Q")
                .arg("20")
                .arg("-r")
                .arg(region)
                .arg(bam)
                .stdout(Stdio::from(append_file));

            let output = command.output().map_err(|err| {
                if err.kind() == std::io::ErrorKind::NotFound {
                    AppError::CommandNotFound {
                        command: self.samtools.clone(),
                    }
                } else {
                    AppError::Io(err)
                }
            })?;

            if !output.status.success() {
                return Err(AppError::CommandFailed {
                    command: format!("{} mpileup -r {}", self.samtools, region),
                    code: output.status.code(),
                    stderr: String::from_utf8_lossy(&output.stderr).to_string(),
                });
            }
        }

        Ok(temp)
    }

    pub fn spawn_samtools_mpileup_stream(
        &self,
        bam: &str,
        fasta: &str,
        region: Option<&str>,
        target_bed: Option<&str>,
    ) -> Result<CommandStream> {
        let stderr_capture = Builder::new()
            .prefix("bam2seqz_cmd_stderr_")
            .suffix(".log")
            .tempfile_in("tmp")?;
        let stderr_file = stderr_capture.reopen()?;

        let mut command = Command::new(&self.samtools);
        command
            .arg("mpileup")
            .arg("-f")
            .arg(fasta)
            .arg("-q")
            .arg("20")
            .arg("-Q")
            .arg("20");

        if let Some(path) = target_bed {
            command.arg("-l").arg(path);
        }

        let command_label = if let Some(value) = region {
            command.arg("-r").arg(value).arg(bam);
            if let Some(path) = target_bed {
                format!(
                    "{} mpileup -l {} -r {} {}",
                    self.samtools, path, value, bam
                )
            } else {
                format!("{} mpileup -r {} {}", self.samtools, value, bam)
            }
        } else {
            command.arg(bam);
            if let Some(path) = target_bed {
                format!("{} mpileup -l {} {}", self.samtools, path, bam)
            } else {
                format!("{} mpileup {}", self.samtools, bam)
            }
        };
        info!(command = %command_label, "spawning mpileup stream");

        command
            .stdout(Stdio::piped())
            .stderr(Stdio::from(stderr_file));

        let mut child = command.spawn().map_err(|err| {
            if err.kind() == std::io::ErrorKind::NotFound {
                AppError::CommandNotFound {
                    command: self.samtools.clone(),
                }
            } else {
                AppError::Io(err)
            }
        })?;

        let stdout = child.stdout.take().ok_or_else(|| AppError::ParseError {
            message: format!("failed to capture stdout for command: {command_label}"),
        })?;

        Ok(CommandStream {
            child,
            stdout,
            stderr_capture,
            command: command_label,
        })
    }

    pub fn bgzip_compress_to(&self, input_path: &str, output_path: &str) -> Result<()> {
        let output_file = File::create(output_path)?;
        let mut command = Command::new(self.bgzip());
        command
            .arg("-c")
            .arg(input_path)
            .stdout(Stdio::from(output_file));
        let status = command.status().map_err(|err| {
            if err.kind() == std::io::ErrorKind::NotFound {
                AppError::CommandNotFound {
                    command: self.bgzip(),
                }
            } else {
                AppError::Io(err)
            }
        })?;

        if !status.success() {
            return Err(AppError::CommandFailed {
                command: format!("{} -c {}", self.bgzip(), input_path),
                code: status.code(),
                stderr: String::new(),
            });
        }
        Ok(())
    }

    pub fn tabix_index_seqz(&self, seqz_gz_path: &str) -> Result<()> {
        let status = Command::new(&self.tabix)
            .arg("-f")
            .arg("-s")
            .arg("1")
            .arg("-b")
            .arg("2")
            .arg("-e")
            .arg("2")
            .arg("-S")
            .arg("1")
            .arg(seqz_gz_path)
            .status()
            .map_err(|err| {
                if err.kind() == std::io::ErrorKind::NotFound {
                    AppError::CommandNotFound {
                        command: self.tabix.clone(),
                    }
                } else {
                    AppError::Io(err)
                }
            })?;

        if !status.success() {
            return Err(AppError::CommandFailed {
                command: format!("{} -f -s 1 -b 2 -e 2 -S 1 {}", self.tabix, seqz_gz_path),
                code: status.code(),
                stderr: String::new(),
            });
        }
        Ok(())
    }

    fn capture_stdout_lines(&self, command: Command, command_name: &str) -> Result<Vec<String>> {
        let stdout_temp = self.capture_stdout_to_tempfile(command, command_name)?;

        self.read_text_lines(stdout_temp.path().to_string_lossy().as_ref())
    }

    fn capture_stdout_to_tempfile(
        &self,
        mut command: Command,
        command_name: &str,
    ) -> Result<NamedTempFile> {
        let stdout_temp = Builder::new()
            .prefix("bam2seqz_cmd_")
            .suffix(".txt")
            .tempfile_in("tmp")?;
        let stdout_file = stdout_temp.reopen()?;
        command.stdout(Stdio::from(stdout_file));

        let output = command.output().map_err(|err| {
            if err.kind() == std::io::ErrorKind::NotFound {
                AppError::CommandNotFound {
                    command: command_name.to_string(),
                }
            } else {
                AppError::Io(err)
            }
        })?;

        if !output.status.success() {
            return Err(AppError::CommandFailed {
                command: command_name.to_string(),
                code: output.status.code(),
                stderr: String::from_utf8_lossy(&output.stderr).to_string(),
            });
        }

        Ok(stdout_temp)
    }

    fn collect_non_empty_lines<R: BufRead>(mut reader: R) -> Result<Vec<String>> {
        let mut lines = Vec::new();
        let mut buf = String::new();
        loop {
            buf.clear();
            let read = reader.read_line(&mut buf)?;
            if read == 0 {
                break;
            }
            let trimmed = buf.trim_end_matches(['\n', '\r']);
            if !trimmed.trim().is_empty() {
                lines.push(trimmed.to_string());
            }
        }
        Ok(lines)
    }

    pub fn exists_command(&self, cmd: &str) -> bool {
        if Path::new(cmd).exists() {
            return true;
        }
        std::env::var_os("PATH").is_some_and(|paths| {
            std::env::split_paths(&paths)
                .map(|dir| dir.join(cmd))
                .any(|full| full.exists())
        })
    }
}

#[cfg(test)]
mod tests {
    use super::ExternalTools;
    use crate::cli::parse_args;

    #[test]
    fn builds_mpileup_command_with_regions() {
        let args = parse_args([
            "bam2seqz",
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
        ])
        .expect("expected parse success");
        let tools = ExternalTools::from_args(&args);
        let command = tools.bam_mpileup_command("n.bam", "ref.fa", &args.chr);

        assert!(command.iter().any(|value| value == "mpileup"));
        assert!(command.iter().any(|value| value == "7"));
        assert!(command.iter().any(|value| value == "12"));
    }
}

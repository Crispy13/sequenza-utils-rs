use std::fs::File;
use std::io::{BufWriter, Write, stdout};
use std::process::{Command, Stdio};

use crate::errors::Result;
use crate::external_tools::ExternalTools;
use crate::seqz_core::seqz_header;

pub fn write_seqz_header<W: Write + ?Sized>(writer: &mut W) -> Result<()> {
    writer.write_all(seqz_header().join("\t").as_bytes())?;
    writer.write_all(b"\n")?;
    Ok(())
}

pub fn with_text_output_writer<F>(path: &str, tools: &ExternalTools, write_fn: F) -> Result<()>
where
    F: FnOnce(&mut dyn Write) -> Result<()>,
{
    with_text_output_writer_with_threads(path, tools, 1, write_fn)
}

pub fn with_text_output_writer_with_threads<F>(
    path: &str,
    tools: &ExternalTools,
    compression_threads: usize,
    write_fn: F,
) -> Result<()>
where
    F: FnOnce(&mut dyn Write) -> Result<()>,
{
    let mut write_fn = Some(write_fn);

    if path == "-" {
        let mut out = stdout().lock();
        let run = write_fn
            .take()
            .ok_or_else(|| std::io::Error::other("writer callback already consumed"))?;
        return run(&mut out);
    }

    if path.ends_with(".gz") {
        let output_file = File::create(path)?;
        let mut command = Command::new(tools.bgzip());
        command.arg("-c");
        if compression_threads > 1 {
            command.arg("-@").arg(compression_threads.to_string());
        }
        command
            .stdin(Stdio::piped())
            .stdout(Stdio::from(output_file));

        let mut child = command.spawn().map_err(|err| {
            if err.kind() == std::io::ErrorKind::NotFound {
                std::io::Error::new(
                    std::io::ErrorKind::NotFound,
                    format!("required command not found in PATH: {}", tools.bgzip()),
                )
            } else {
                err
            }
        })?;

        {
            let mut stdin = child
                .stdin
                .take()
                .ok_or_else(|| std::io::Error::other("failed to open bgzip stdin"))?;
            let run = write_fn
                .take()
                .ok_or_else(|| std::io::Error::other("writer callback already consumed"))?;
            run(&mut stdin)?;
            stdin.flush()?;
        }

        let status = child.wait()?;
        if !status.success() {
            return Err(std::io::Error::other(format!(
                "bgzip compression failed with status: {status}"
            ))
            .into());
        }
        return Ok(());
    }

    let mut file = BufWriter::new(File::create(path)?);
    let run = write_fn
        .take()
        .ok_or_else(|| std::io::Error::other("writer callback already consumed"))?;
    run(&mut file)?;
    file.flush()?;
    Ok(())
}

pub fn write_text_output(path: &str, content: &str, tools: &ExternalTools) -> Result<()> {
    with_text_output_writer(path, tools, |out| {
        out.write_all(content.as_bytes())?;
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use super::write_seqz_header;

    #[test]
    fn writes_tab_delimited_header() {
        let mut output = Vec::new();
        write_seqz_header(&mut output).expect("expected header write success");
        let line = String::from_utf8(output).expect("expected utf8 output");
        assert!(line.starts_with("chromosome\tposition\tbase.ref"));
        assert!(line.ends_with("tumor.strand\n"));
    }
}

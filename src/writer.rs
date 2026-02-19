use std::fs::File;
use std::io::{BufWriter, Write, stdout};
use std::path::Path;

use crate::errors::Result;
use crate::external_tools::ExternalTools;
use crate::seqz_core::seqz_header;
use tempfile::Builder;

pub fn write_seqz_header<W: Write + ?Sized>(writer: &mut W) -> Result<()> {
    writer.write_all(seqz_header().join("\t").as_bytes())?;
    writer.write_all(b"\n")?;
    Ok(())
}

pub fn with_text_output_writer<F>(path: &str, tools: &ExternalTools, write_fn: F) -> Result<()>
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
        let output_path = Path::new(path);
        let parent_dir = output_path.parent().unwrap_or_else(|| Path::new("."));
        let mut plain_file = Builder::new()
            .prefix("bam2seqz_")
            .suffix(".tmp_plain")
            .tempfile_in(parent_dir)?;

        {
            let mut buf = BufWriter::new(plain_file.as_file_mut());
            let run = write_fn
                .take()
                .ok_or_else(|| std::io::Error::other("writer callback already consumed"))?;
            run(&mut buf)?;
            buf.flush()?;
        }

        let plain_path = plain_file.path().to_string_lossy().into_owned();
        tools.bgzip_compress_to(&plain_path, path)?;
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

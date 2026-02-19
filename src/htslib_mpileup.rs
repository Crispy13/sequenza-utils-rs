use crate::cli::Bam2SeqzArgs;
use crate::errors::{AppError, Result};
use bio::io::fasta;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::fmt::Write as _;
use std::fs::File;

const MIN_MAPQ: u8 = 20;
const MIN_BASEQ: u8 = 20;

#[derive(Debug, Clone)]
pub(crate) struct HtslibPileupRecord {
    pub chromosome: String,
    pub position: i32,
    pub reference: String,
    pub depth: i32,
    pub pileup: String,
    pub quality: String,
}

impl HtslibPileupRecord {
    pub fn write_data_line(&self, buffer: &mut String) {
        buffer.clear();
        buffer.push_str(&self.reference);
        buffer.push('\t');
        let _ = write!(buffer, "{}", self.depth);
        buffer.push('\t');
        buffer.push_str(&self.pileup);
        buffer.push('\t');
        buffer.push_str(&self.quality);
    }
}

#[derive(Debug, Clone)]
pub(crate) struct RegionSpec {
    pub chrom: String,
    pub start0: u64,
    pub end0: u64,
}

pub(crate) struct HtslibReadersMut<'a> {
    pub normal: &'a mut bam::IndexedReader,
    pub tumor: &'a mut bam::IndexedReader,
    pub normal2: Option<&'a mut bam::IndexedReader>,
    pub fasta: &'a mut fasta::IndexedReader<File>,
}

#[derive(Debug)]
pub(crate) struct HtslibWorkerContext {
    normal_path: String,
    tumor_path: String,
    normal2_path: Option<String>,
    fasta_path: String,
    normal_reader: Option<bam::IndexedReader>,
    tumor_reader: Option<bam::IndexedReader>,
    normal2_reader: Option<bam::IndexedReader>,
    fasta_reader: Option<fasta::IndexedReader<File>>,
}

impl HtslibWorkerContext {
    pub fn from_args(args: &Bam2SeqzArgs) -> Self {
        Self {
            normal_path: args.normal.clone(),
            tumor_path: args.tumor.clone(),
            normal2_path: args.normal2.clone(),
            fasta_path: args.fasta.clone().unwrap_or_default(),
            normal_reader: None,
            tumor_reader: None,
            normal2_reader: None,
            fasta_reader: None,
        }
    }

    pub fn resolve_regions(&mut self, requested: &[String]) -> Result<Vec<String>> {
        if !requested.is_empty() {
            return Ok(requested.to_vec());
        }

        self.ensure_open()?;
        let tumor_reader = self.tumor_reader.as_ref().ok_or_else(|| AppError::ParseError {
            message: "missing tumor BAM reader in htslib worker context".to_string(),
        })?;

        let header = tumor_reader.header();
        let mut regions = Vec::new();
        for tid in 0..header.target_count() {
            if header.target_len(tid).is_some_and(|value| value > 0) {
                regions.push(String::from_utf8_lossy(header.tid2name(tid)).into_owned());
            }
        }
        Ok(regions)
    }

    pub fn readers_mut(
        &mut self,
    ) -> Result<HtslibReadersMut<'_>> {
        self.ensure_open()?;
        let normal = self
            .normal_reader
            .as_mut()
            .ok_or_else(|| AppError::ParseError {
                message: "missing normal BAM reader in htslib worker context".to_string(),
            })?;
        let tumor = self
            .tumor_reader
            .as_mut()
            .ok_or_else(|| AppError::ParseError {
                message: "missing tumor BAM reader in htslib worker context".to_string(),
            })?;
        let normal2 = self.normal2_reader.as_mut();
        let fasta = self
            .fasta_reader
            .as_mut()
            .ok_or_else(|| AppError::ParseError {
                message: "missing FASTA reader in htslib worker context".to_string(),
            })?;

        Ok(HtslibReadersMut {
            normal,
            tumor,
            normal2,
            fasta,
        })
    }

    fn ensure_open(&mut self) -> Result<()> {
        if self.normal_reader.is_none() {
            self.normal_reader = Some(bam::IndexedReader::from_path(&self.normal_path).map_err(
                |err| AppError::ParseError {
                    message: format!(
                        "failed to open normal BAM index reader for {}: {err}",
                        self.normal_path
                    ),
                },
            )?);
        }

        if self.tumor_reader.is_none() {
            self.tumor_reader = Some(bam::IndexedReader::from_path(&self.tumor_path).map_err(
                |err| AppError::ParseError {
                    message: format!(
                        "failed to open tumor BAM index reader for {}: {err}",
                        self.tumor_path
                    ),
                },
            )?);
        }

        if self.fasta_reader.is_none() {
            self.fasta_reader = Some(
                fasta::IndexedReader::from_file(&self.fasta_path).map_err(|err| {
                    AppError::ParseError {
                        message: format!(
                            "failed to open indexed FASTA reader for {}: {err}",
                            self.fasta_path
                        ),
                    }
                })?,
            );
        }

        if self.normal2_reader.is_none()
            && let Some(path) = self.normal2_path.as_deref()
        {
            self.normal2_reader = Some(bam::IndexedReader::from_path(path).map_err(|err| {
                AppError::ParseError {
                    message: format!("failed to open normal2 BAM index reader for {}: {err}", path),
                }
            })?);
        }

        Ok(())
    }
}

pub(crate) fn parse_region_spec(reader: &bam::IndexedReader, region: &str) -> Result<RegionSpec> {
    if let Some((chrom, bounds)) = region.split_once(':') {
        let (start_raw, end_raw) = bounds.split_once('-').ok_or_else(|| AppError::ParseError {
            message: format!("invalid region bounds in {region}"),
        })?;

        let start1 = start_raw
            .replace(',', "")
            .parse::<u64>()
            .map_err(|_| AppError::ParseError {
                message: format!("invalid region start in {region}"),
            })?;
        let end1 = end_raw
            .replace(',', "")
            .parse::<u64>()
            .map_err(|_| AppError::ParseError {
                message: format!("invalid region end in {region}"),
            })?;

        if start1 == 0 || end1 < start1 {
            return Err(AppError::ParseError {
                message: format!("invalid region coordinates in {region}"),
            });
        }

        let target_len = target_length(reader, chrom)?;
        let bounded_end = end1.min(target_len);
        if start1 > bounded_end {
            return Err(AppError::ParseError {
                message: format!("region start exceeds contig length in {region}"),
            });
        }

        return Ok(RegionSpec {
            chrom: chrom.to_string(),
            start0: start1 - 1,
            end0: bounded_end,
        });
    }

    let target_len = target_length(reader, region)?;
    Ok(RegionSpec {
        chrom: region.to_string(),
        start0: 0,
        end0: target_len,
    })
}

pub(crate) fn fetch_reference_bases(
    reader: &mut fasta::IndexedReader<File>,
    spec: &RegionSpec,
) -> Result<Vec<u8>> {
    reader
        .fetch(&spec.chrom, spec.start0, spec.end0)
        .map_err(|err| AppError::ParseError {
            message: format!(
                "failed FASTA fetch for {}:{}-{}: {err}",
                spec.chrom,
                spec.start0 + 1,
                spec.end0
            ),
        })?;

    let mut reference_bases = Vec::new();
    reader
        .read(&mut reference_bases)
        .map_err(|err| AppError::ParseError {
            message: format!(
                "failed FASTA read for {}:{}-{}: {err}",
                spec.chrom,
                spec.start0 + 1,
                spec.end0
            ),
        })?;

    Ok(reference_bases)
}

pub(crate) fn fetch_region(reader: &mut bam::IndexedReader, spec: &RegionSpec) -> Result<()> {
    let tid = reader
        .header()
        .tid(spec.chrom.as_bytes())
        .ok_or_else(|| AppError::ParseError {
            message: format!("chromosome not found in BAM header: {}", spec.chrom),
        })?;

    reader
        .fetch((tid, spec.start0, spec.end0))
        .map_err(|err| AppError::ParseError {
            message: format!(
                "failed BAM fetch for {}:{}-{}: {err}",
                spec.chrom,
                spec.start0 + 1,
                spec.end0
            ),
        })?;
    Ok(())
}

pub(crate) fn next_record_from_pileups(
    pileups: &mut bam::pileup::Pileups<'_, bam::IndexedReader>,
    spec: &RegionSpec,
    reference_bases: &[u8],
) -> Result<Option<HtslibPileupRecord>> {
    let mut pileup_line = String::with_capacity(64);
    let mut quality_line = String::with_capacity(64);

    for pileup_result in pileups {
        let pileup = pileup_result.map_err(|err| AppError::ParseError {
            message: format!("failed pileup iteration for {}: {err}", spec.chrom),
        })?;

        let position0 = u64::from(pileup.pos());
        if position0 < spec.start0 || position0 >= spec.end0 {
            continue;
        }

        let Some(reference_base) = reference_bases.get((position0 - spec.start0) as usize)
        else {
            continue;
        };

        pileup_line.clear();
        quality_line.clear();
        let mut depth = 0_i32;
        let normalized_ref = reference_base.to_ascii_uppercase();

        for alignment in pileup.alignments() {
            if alignment.is_refskip() {
                continue;
            }

            let record = alignment.record();
            if record.mapq() < MIN_MAPQ
                || record.is_unmapped()
                || record.is_duplicate()
                || record.is_secondary()
                || record.is_supplementary()
                || record.is_quality_check_failed()
            {
                continue;
            }

            if alignment.is_del() {
                depth += 1;
                pileup_line.push('*');
                quality_line.push('!');
                continue;
            }

            let Some(query_position) = alignment.qpos() else {
                continue;
            };
            let Some(base) = record.seq().as_bytes().get(query_position).copied() else {
                continue;
            };
            let Some(quality) = record.qual().get(query_position).copied() else {
                continue;
            };

            if quality < MIN_BASEQ {
                continue;
            }

            depth += 1;
            pileup_line.push(pileup_base_char(base, normalized_ref, record.is_reverse()));
            quality_line.push(char::from(quality.saturating_add(33)));
        }

        if depth == 0 {
            continue;
        }

        let position = i32::try_from(position0 + 1).map_err(|_| AppError::ParseError {
            message: format!("position exceeds supported i32 range: {}", position0 + 1),
        })?;

        return Ok(Some(HtslibPileupRecord {
            chromosome: spec.chrom.clone(),
            position,
            reference: char::from(normalized_ref).to_string(),
            depth,
            pileup: pileup_line,
            quality: quality_line,
        }));
    }

    Ok(None)
}

fn target_length(reader: &bam::IndexedReader, chrom: &str) -> Result<u64> {
    let header = reader.header();
    let tid = header.tid(chrom.as_bytes()).ok_or_else(|| AppError::ParseError {
        message: format!("chromosome not found in BAM header: {chrom}"),
    })?;

    header.target_len(tid).ok_or_else(|| AppError::ParseError {
        message: format!("missing target length for chromosome: {chrom}"),
    })
}

fn pileup_base_char(base: u8, reference_base: u8, reverse: bool) -> char {
    let normalized_base = normalize_base(base);
    if normalized_base == reference_base {
        return if reverse { ',' } else { '.' };
    }

    if reverse {
        char::from(normalized_base.to_ascii_lowercase())
    } else {
        char::from(normalized_base)
    }
}

fn normalize_base(base: u8) -> u8 {
    let upper = base.to_ascii_uppercase();
    match upper {
        b'A' | b'C' | b'G' | b'T' | b'N' => upper,
        _ => b'N',
    }
}
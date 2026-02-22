pub const SEQZ_HEADER: [&str; 14] = [
    "chromosome",
    "position",
    "base.ref",
    "depth.normal",
    "depth.tumor",
    "depth.ratio",
    "Af",
    "Bf",
    "zygosity.normal",
    "GC.percent",
    "good.reads",
    "AB.normal",
    "AB.tumor",
    "tumor.strand",
];

pub fn seqz_header() -> &'static [&'static str; 14] {
    &SEQZ_HEADER
}

#[derive(Debug, Clone, PartialEq)]
pub struct SeqzParams {
    pub depth_sum: i32,
    pub qlimit: u8,
    pub hom_t: f64,
    pub het_t: f64,
    pub het_f: f64,
    pub het_only: bool,
}

impl Default for SeqzParams {
    fn default() -> Self {
        Self {
            depth_sum: 20,
            qlimit: 53,
            hom_t: 0.85,
            het_t: 0.35,
            het_f: -0.1,
            het_only: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ACGTCounts {
    counts: [i32; 4],
    strand: [i32; 4],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeqzInput<'a> {
    pub reference: &'a str,
    pub normal_depth: i32,
    pub normal_pileup: &'a str,
    pub normal_quality: &'a str,
    pub tumor_depth: i32,
    pub tumor_pileup: &'a str,
    pub tumor_quality: &'a str,
    pub gc: &'a str,
    pub normal2_depth_override: Option<i32>,
}

pub fn do_seqz(data: &[&str], params: &SeqzParams) -> Option<Vec<String>> {
    let unpacked = unpack_data(data)?;
    do_seqz_typed(&unpacked, params)
}

pub fn do_seqz_typed(input: &SeqzInput<'_>, params: &SeqzParams) -> Option<Vec<String>> {
    let reference = parse_reference_base(input.reference)?;
    let normal_depth = input
        .normal2_depth_override
        .filter(|depth| *depth > 0)
        .unwrap_or(input.normal_depth);
    if normal_depth + input.tumor_depth < params.depth_sum {
        return None;
    }

    let normal = acgt(
        input.normal_pileup,
        input.normal_quality,
        reference,
        params.qlimit,
    );
    let tumor = acgt(
        input.tumor_pileup,
        input.tumor_quality,
        reference,
        params.qlimit,
    );

    acgt_to_seqz(
        normal,
        tumor,
        params,
        normal_depth,
        input.tumor_depth,
        reference,
        input.gc,
    )
}

fn parse_reference_base(reference: &str) -> Option<char> {
    reference
        .chars()
        .next()
        .map(|base| base.to_ascii_uppercase())
}

fn unpack_data<'a>(data: &'a [&'a str]) -> Option<SeqzInput<'a>> {
    if data.len() != 3 && data.len() != 4 {
        return None;
    }

    let normal_line = data[0];
    let tumor_line = data[1];
    let gc = data[2];

    let mut normal_parts = normal_line.split('\t');
    let reference = normal_parts.next()?;
    let normal_depth = normal_parts.next()?.parse::<i32>().ok()?;
    let normal_pileup = normal_parts.next()?;
    let normal_quality = normal_parts.next()?;

    let mut tumor_parts = tumor_line.split('\t');
    let _tumor_ref = tumor_parts.next()?;
    let tumor_depth = tumor_parts.next()?.parse::<i32>().ok()?;
    let tumor_pileup = tumor_parts.next()?;
    let tumor_quality = tumor_parts.next()?;

    let normal2_depth_override = if data.len() == 4 {
        data[3]
            .split('\t')
            .nth(1)
            .and_then(|value| value.parse::<i32>().ok())
            .filter(|depth| *depth > 0)
    } else {
        None
    };

    Some(SeqzInput {
        reference,
        normal_depth,
        normal_pileup,
        normal_quality,
        tumor_depth,
        tumor_pileup,
        tumor_quality,
        gc,
        normal2_depth_override,
    })
}

fn acgt(pileup: &str, quality: &str, reference: char, qlimit: u8) -> ACGTCounts {
    let mut counts = [0_i32; 4];
    let mut strand = [0_i32; 4];

    let pileup_bytes = pileup.as_bytes();
    let quality_bytes = quality.as_bytes();
    let mut n = 0usize;
    let mut q = 0usize;
    let mut last_base: Option<usize> = None;

    while n < pileup_bytes.len() {
        let mut base = pileup_bytes[n] as char;

        if base == '^' {
            n += 2;
            if n >= pileup_bytes.len() {
                break;
            }
            base = pileup_bytes[n] as char;
        }

        if base == '.' {
            base = reference;
        } else if base == ',' {
            base = reference.to_ascii_lowercase();
        }

        match base {
            'A' | 'C' | 'G' | 'T' => {
                let idx = base_to_index(base).expect("valid base index");
                if q < quality_bytes.len() && quality_bytes[q] >= qlimit {
                    counts[idx] += 1;
                    strand[idx] += 1;
                    last_base = Some(idx);
                } else {
                    last_base = None;
                }
                n += 1;
                q += 1;
            }
            'a' | 'c' | 'g' | 't' => {
                let idx = base_to_index(base.to_ascii_uppercase()).expect("valid base index");
                if q < quality_bytes.len() && quality_bytes[q] >= qlimit {
                    counts[idx] += 1;
                    last_base = Some(idx + 4);
                } else {
                    last_base = None;
                }
                n += 1;
                q += 1;
            }
            '$' => {
                last_base = None;
                n += 1;
            }
            '*' => {
                last_base = None;
                n += 1;
                q += 1;
            }
            '+' | '-' => {
                let mut offset = n + 1;
                while offset < pileup_bytes.len() && (pileup_bytes[offset] as char).is_ascii_digit()
                {
                    offset += 1;
                }
                if offset <= n + 1 {
                    n += 1;
                    continue;
                }
                let step = std::str::from_utf8(&pileup_bytes[(n + 1)..offset])
                    .ok()
                    .and_then(|value| value.parse::<usize>().ok())
                    .unwrap_or(0);
                n = offset + step;
            }
            _ => {
                if last_base.is_some() {
                    last_base = None;
                }
                n += 1;
            }
        }
    }

    ACGTCounts { counts, strand }
}

fn acgt_to_seqz(
    normal: ACGTCounts,
    tumor: ACGTCounts,
    params: &SeqzParams,
    normal_depth: i32,
    tumor_depth: i32,
    reference: char,
    gc: &str,
) -> Option<Vec<String>> {
    let sum_normal: i32 = normal.counts.iter().sum();
    if sum_normal <= 0 {
        return None;
    }
    let sum_tumor: i32 = tumor.counts.iter().sum();
    if sum_tumor <= 0 {
        return None;
    }

    let normal_freq = frequencies(&normal.counts, sum_normal);
    let tumor_freq = frequencies(&tumor.counts, sum_tumor);
    let alleles = acgt_genotype(
        &normal.counts,
        &normal_freq,
        &normal.strand,
        params.hom_t,
        params.het_t,
        params.het_f,
    );

    if alleles.len() == 1 && !params.het_only {
        Some(parse_homoz(
            &tumor,
            &tumor_freq,
            sum_tumor,
            &normal,
            normal_depth,
            tumor_depth,
            alleles[0],
            reference,
            gc,
        ))
    } else if alleles.len() == 2 {
        let mut sorted_alleles = alleles;
        if tumor_freq[sorted_alleles[0]] < tumor_freq[sorted_alleles[1]] {
            sorted_alleles.swap(0, 1);
        }
        Some(parse_heteroz(
            &tumor_freq,
            sum_tumor,
            normal_depth,
            tumor_depth,
            sorted_alleles,
            reference,
            gc,
        ))
    } else {
        None
    }
}

fn frequencies(values: &[i32; 4], sum: i32) -> [f64; 4] {
    [
        values[0] as f64 / sum as f64,
        values[1] as f64 / sum as f64,
        values[2] as f64 / sum as f64,
        values[3] as f64 / sum as f64,
    ]
}

fn acgt_genotype(
    counts: &[i32; 4],
    freq_list: &[f64; 4],
    strand_list: &[i32; 4],
    hom_t: f64,
    het_t: f64,
    het_f: f64,
) -> Vec<usize> {
    let mut indices = [0usize, 1, 2, 3];
    indices.sort_by(|lhs, rhs| {
        freq_list[*rhs]
            .partial_cmp(&freq_list[*lhs])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let first = indices[0];
    let second = indices[1];
    let mut alleles = vec![first];

    if freq_list[first] < hom_t && freq_list[second] >= het_t {
        let alt_count = counts[second];
        if alt_count > 0 {
            let alt_fw_freq = strand_list[second] as f64 / alt_count as f64;
            if het_f < alt_fw_freq && alt_fw_freq < (1.0 - het_f) {
                alleles.push(second);
            }
        }
    }
    alleles
}

#[allow(clippy::too_many_arguments)]
fn parse_homoz(
    tumor: &ACGTCounts,
    tumor_freq: &[f64; 4],
    sum_tumor: i32,
    normal: &ACGTCounts,
    normal_depth: i32,
    tumor_depth: i32,
    allele_idx: usize,
    reference: char,
    gc: &str,
) -> Vec<String> {
    let mut no_zero_bases = Vec::new();
    let mut strand_bases = Vec::new();

    for (idx, freq) in tumor_freq.iter().enumerate().take(4) {
        if *freq > 0.0 && normal.counts[idx] == 0 && idx != allele_idx {
            no_zero_bases.push(format!("{}{}", index_to_base(idx), py_str_round3(*freq)));
            if tumor.counts[idx] > 0 {
                let strand_ratio = tumor.strand[idx] as f64 / tumor.counts[idx] as f64;
                strand_bases.push(format!(
                    "{}{}",
                    index_to_base(idx),
                    py_str_round3(strand_ratio)
                ));
            }
        }
    }

    let no_zero = if no_zero_bases.is_empty() {
        ".".to_string()
    } else {
        no_zero_bases.join(":")
    };
    let strand_repr = if strand_bases.is_empty() {
        "0".to_string()
    } else {
        strand_bases.join(":")
    };

    vec![
        reference.to_string(),
        normal_depth.to_string(),
        tumor_depth.to_string(),
        depth_ratio_str(normal_depth, tumor_depth, true),
        py_str_round3(tumor_freq[allele_idx]),
        "0".to_string(),
        "hom".to_string(),
        gc.to_string(),
        sum_tumor.to_string(),
        index_to_base(allele_idx).to_string(),
        no_zero,
        strand_repr,
    ]
}

fn parse_heteroz(
    tumor_freq: &[f64; 4],
    sum_tumor: i32,
    normal_depth: i32,
    tumor_depth: i32,
    alleles: Vec<usize>,
    reference: char,
    gc: &str,
) -> Vec<String> {
    let genotype = format!("{}{}", index_to_base(alleles[0]), index_to_base(alleles[1]));
    let mut selected = [tumor_freq[alleles[0]], tumor_freq[alleles[1]]];
    selected.sort_by(|lhs, rhs| rhs.partial_cmp(lhs).unwrap_or(std::cmp::Ordering::Equal));

    vec![
        reference.to_string(),
        normal_depth.to_string(),
        tumor_depth.to_string(),
        depth_ratio_str(normal_depth, tumor_depth, false),
        py_str_round3(selected[0]),
        py_str_round3(selected[1]),
        "het".to_string(),
        gc.to_string(),
        sum_tumor.to_string(),
        genotype,
        ".".to_string(),
        "0".to_string(),
    ]
}

fn base_to_index(base: char) -> Option<usize> {
    match base {
        'A' => Some(0),
        'C' => Some(1),
        'G' => Some(2),
        'T' => Some(3),
        _ => None,
    }
}

fn index_to_base(index: usize) -> char {
    match index {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        _ => 'T',
    }
}

fn depth_ratio_str(normal_depth: i32, tumor_depth: i32, strip_trailing: bool) -> String {
    let n = normal_depth as i64;
    let t = tumor_depth as i64;
    if n == 0 {
        return "0.000".to_string();
    }

    let (q, rem) = ((t * 1000) / n, (t * 1000) % n);
    let rounded = if strip_trailing {
        if rem > n / 2 {
            q + 1
        } else if rem < n / 2 {
            q
        } else if n % 2 == 0 {
            if q % 2 == 0 {
                q
            } else {
                q + 1
            }
        } else {
            q
        }
    } else {
        (t * 1000 + n / 2) / n
    };

    let int_part = rounded / 1000;
    let frac = rounded % 1000;
    let text = format!("{int_part}.{frac:03}");
    if strip_trailing {
        let mut trimmed = text.trim_end_matches('0').to_string();
        if trimmed.ends_with('.') {
            trimmed.push('0');
        }
        trimmed
    } else {
        text
    }
}

fn py_str_round3(value: f64) -> String {
    let rounded = (value * 1000.0).round_ties_even() / 1000.0;
    if (rounded.fract()).abs() < f64::EPSILON {
        return format!("{rounded:.1}");
    }
    let mut text = format!("{rounded:.3}");
    while text.ends_with('0') {
        text.pop();
    }
    if text.ends_with('.') {
        text.push('0');
    }
    text
}

#[cfg(test)]
mod tests {
    use super::{SeqzInput, SeqzParams, do_seqz, do_seqz_typed, py_str_round3, seqz_header};

    fn typed_input_from_lines<'a>(
        normal_line: &'a str,
        tumor_line: &'a str,
        gc: &'a str,
        normal2_line: Option<&'a str>,
    ) -> SeqzInput<'a> {
        let mut normal_parts = normal_line.split('\t');
        let reference = normal_parts.next().expect("normal reference");
        let normal_depth = normal_parts
            .next()
            .expect("normal depth")
            .parse::<i32>()
            .expect("valid normal depth");
        let normal_pileup = normal_parts.next().expect("normal pileup");
        let normal_quality = normal_parts.next().expect("normal quality");

        let mut tumor_parts = tumor_line.split('\t');
        let _ = tumor_parts.next().expect("tumor reference");
        let tumor_depth = tumor_parts
            .next()
            .expect("tumor depth")
            .parse::<i32>()
            .expect("valid tumor depth");
        let tumor_pileup = tumor_parts.next().expect("tumor pileup");
        let tumor_quality = tumor_parts.next().expect("tumor quality");

        let normal2_depth_override = normal2_line
            .and_then(|line| line.split('\t').nth(1))
            .and_then(|depth| depth.parse::<i32>().ok())
            .filter(|depth| *depth > 0);

        SeqzInput {
            reference,
            normal_depth,
            normal_pileup,
            normal_quality,
            tumor_depth,
            tumor_pileup,
            tumor_quality,
            gc,
            normal2_depth_override,
        }
    }

    #[test]
    fn seqz_header_shape_matches_python() {
        let header = seqz_header();
        assert_eq!(header.len(), 14);
        assert_eq!(header[0], "chromosome");
        assert_eq!(header[1], "position");
        assert_eq!(header[13], "tumor.strand");
    }

    #[test]
    fn do_seqz_matches_python_reference_cases() {
        let normal = "T\t29\t,C.C,c,C,c,,,c,cCccC,c,,c,c,,\tBB/<FFFBFFFFFFBFFFF/7/7FBFFFF";
        let tumor = "T\t46\tc$ccc,cCcc.cc,cGcC.Ccc,c,C.CC,.CcC.ccc,Cc,ccccg\t/FFFFFgBF/FF/F/F<//FF/FFk</BF/<F/BF/FFB/FBFF</";
        let line = do_seqz(&[normal, tumor, "50"], &SeqzParams::default()).expect("line expected");
        assert_eq!(line[9], "CT");

        let normal_hom = "T\t29\t,...,,,.,,,,,,,,.,,.,,,,,,,,,\tBB/<FFFBFFFFFFBFFFF/7/7FBFFFF";
        let line_hom =
            do_seqz(&[normal_hom, tumor, "50"], &SeqzParams::default()).expect("line expected");
        assert_eq!(line_hom[9], "T");
        assert_eq!(line_hom[10], "C0.788");
        assert_eq!(line_hom[11], "C0.231");

        let line_alt = do_seqz(&[normal_hom, tumor, "50", tumor], &SeqzParams::default())
            .expect("line expected");
        assert_eq!(line_alt[1], line_alt[2]);
        assert_eq!(line_alt[9], "T");
        assert_eq!(line_alt[10], "C0.788");
        assert_eq!(line_alt[11], "C0.231");
    }

    #[test]
    fn typed_api_matches_do_seqz_for_representative_cases() {
        let params = SeqzParams::default();
        let normal_het = "T\t29\t,C.C,c,C,c,,,c,cCccC,c,,c,c,,\tBB/<FFFBFFFFFFBFFFF/7/7FBFFFF";
        let normal_hom = "T\t29\t,...,,,.,,,,,,,,.,,.,,,,,,,,,\tBB/<FFFBFFFFFFBFFFF/7/7FBFFFF";
        let tumor = "T\t46\tc$ccc,cCcc.cc,cGcC.Ccc,c,C.CC,.CcC.ccc,Cc,ccccg\t/FFFFFgBF/FF/F/F<//FF/FFk</BF/<F/BF/FFB/FBFF</";

        let cases = [
            (normal_het, None),
            (normal_hom, None),
            (normal_hom, Some(tumor)),
        ];

        for (normal, normal2) in cases {
            let mut data = vec![normal, tumor, "50"];
            if let Some(alt) = normal2 {
                data.push(alt);
            }

            let expected = do_seqz(&data, &params);
            let typed = typed_input_from_lines(normal, tumor, "50", normal2);
            let actual = do_seqz_typed(&typed, &params);
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn py_round3_matches_python_ties_to_even() {
        assert_eq!(py_str_round3(77.0 / 16.0), "4.812");
        assert_eq!(py_str_round3(81.0 / 16.0), "5.062");
        assert_eq!(py_str_round3(89.0 / 16.0), "5.562");
        assert_eq!(py_str_round3(1.0), "1.0");
    }
}

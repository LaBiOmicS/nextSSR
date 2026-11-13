//! Per-sequence analysis pipeline.

use serde::Serialize;

use crate::fasta::Record;
use crate::imperfect::chain_seeds;
use crate::locus::{build_loci, Locus, Ssr};
use crate::params::{Config, Mode};
use crate::primer::design;
use crate::scan::{eq_runs, select_maximal, select_misa, PerfectHit};

/// Result of analysing one FASTA record.
#[derive(Debug, Clone, Serialize)]
pub struct RecordResult {
    pub index: usize,
    /// First token of the header.
    pub id: String,
    /// Header as written by misa.pl in `.misa` files.
    pub misa_id: String,
    pub length: usize,
    pub loci: Vec<Locus>,
}

fn to_ssr(seq: &[u8], h: &PerfectHit) -> Ssr {
    Ssr {
        start: h.start + 1,
        end: h.end,
        unit: h.unit,
        motif: seq[h.start..h.start + h.unit]
            .iter()
            .map(|b| b.to_ascii_uppercase() as char)
            .collect(),
        repeats: h.copies,
        imperfect: None,
    }
}

/// Detect SSRs of one unit size. Returns (perfect, imperfect).
fn scan_unit(seq: &[u8], unit: usize, min_rep: usize, cfg: &Config) -> (Vec<Ssr>, Vec<Ssr>) {
    let seed = cfg
        .imperfect
        .as_ref()
        .map(|p| p.seed_for(unit).min(min_rep));
    let min_copies = seed.unwrap_or(min_rep).min(min_rep);
    let runs = eq_runs(seq, unit, (min_copies - 1) * unit);
    let hits = match cfg.mode {
        Mode::Misa => select_misa(seq, &runs, unit, min_rep),
        Mode::Native => select_maximal(seq, &runs, unit, min_rep),
    };
    let perfect: Vec<Ssr> = hits.iter().map(|h| to_ssr(seq, h)).collect();
    let imperfect = match (&cfg.imperfect, seed) {
        (Some(ip), Some(seed)) => {
            let seeds = select_maximal(seq, &runs, unit, seed);
            chain_seeds(seq, &seeds, unit, min_rep * unit, ip)
        }
        _ => Vec::new(),
    };
    (perfect, imperfect)
}

/// Remove perfect SSRs fully contained in an imperfect SSR of the same unit.
fn merge_imperfect(perfect: Vec<Ssr>, imperfect: Vec<Ssr>) -> Vec<Ssr> {
    if imperfect.is_empty() {
        return perfect;
    }
    let mut out = Vec::with_capacity(perfect.len() + imperfect.len());
    let mut j = 0;
    for p in perfect {
        while j < imperfect.len() && imperfect[j].end < p.start {
            j += 1;
        }
        let contained = imperfect[j..]
            .iter()
            .take_while(|i| i.start <= p.start)
            .any(|i| i.end >= p.end);
        if !contained {
            out.push(p);
        }
    }
    out.extend(imperfect);
    out
}

/// Analyse one record: SSR detection, compound grouping, flanks and primers.
pub fn analyze_record(rec: &Record, cfg: &Config) -> RecordResult {
    let seq = rec.seq.as_slice();
    let units: Vec<(usize, usize)> = cfg
        .definition
        .min_repeats
        .iter()
        .map(|(&u, &r)| (u, r))
        .collect();

    #[cfg(feature = "parallel")]
    let per_unit: Vec<(Vec<Ssr>, Vec<Ssr>)> = {
        use rayon::prelude::*;
        if seq.len() >= 1 << 16 {
            units
                .par_iter()
                .map(|&(u, r)| scan_unit(seq, u, r, cfg))
                .collect()
        } else {
            units
                .iter()
                .map(|&(u, r)| scan_unit(seq, u, r, cfg))
                .collect()
        }
    };
    #[cfg(not(feature = "parallel"))]
    let per_unit: Vec<(Vec<Ssr>, Vec<Ssr>)> = units
        .iter()
        .map(|&(u, r)| scan_unit(seq, u, r, cfg))
        .collect();

    let mut ssrs: Vec<Ssr> = Vec::new();
    for (perfect, imperfect) in per_unit {
        ssrs.extend(merge_imperfect(perfect, imperfect));
    }
    ssrs.sort_by_key(|s| (s.start, s.end, s.unit));

    let mut loci = build_loci(ssrs, seq, cfg.max_interruption, cfg.mode);

    if cfg.flank > 0 {
        for l in &mut loci {
            let a = l.start - 1;
            let f5 = &seq[a.saturating_sub(cfg.flank)..a];
            let f3 = &seq[l.end..(l.end + cfg.flank).min(seq.len())];
            l.flank5 = Some(String::from_utf8_lossy(f5).into_owned());
            l.flank3 = Some(String::from_utf8_lossy(f3).into_owned());
        }
    }

    if let Some(pp) = &cfg.primers {
        let regions: Vec<(usize, usize)> = loci.iter().map(|l| (l.start - 1, l.end)).collect();
        let design_one = |idx: usize| {
            let (a, b) = regions[idx];
            let lo = a.saturating_sub(pp.product_max);
            let hi = b + pp.product_max;
            let excluded: Vec<(usize, usize)> = regions
                .iter()
                .enumerate()
                .filter(|&(j, &(s, e))| j != idx && e > lo && s < hi)
                .map(|(_, &r)| r)
                .collect();
            design(seq, a, b, &excluded, pp)
        };
        #[cfg(feature = "parallel")]
        let outcomes: Vec<_> = {
            use rayon::prelude::*;
            (0..loci.len()).into_par_iter().map(design_one).collect()
        };
        #[cfg(not(feature = "parallel"))]
        let outcomes: Vec<_> = (0..loci.len()).map(design_one).collect();
        for (l, o) in loci.iter_mut().zip(outcomes) {
            l.primers = Some(o);
        }
    }

    RecordResult {
        index: rec.index,
        id: rec.id().to_string(),
        misa_id: rec.misa_id(),
        length: seq.len(),
        loci,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::parse_all;
    use crate::params::ImperfectParams;

    #[test]
    fn native_vs_misa_masking() {
        let fa = b">s1 test\nCCCCGAAAAAAAAAAAAAAAGAAGAAGAAGAAGCCTCG\n";
        let rec = &parse_all(fa)[0];
        let native = analyze_record(rec, &Config::default());
        let misa = analyze_record(
            rec,
            &Config {
                mode: Mode::Misa,
                ..Config::default()
            },
        );
        let pats = |r: &RecordResult| r.loci.iter().map(|l| l.pattern.clone()).collect::<Vec<_>>();
        assert_eq!(pats(&misa), vec!["(A)15"]);
        assert_eq!(pats(&native), vec!["(A)15{-2}(AAG)5"]);
        assert_eq!(native.loci[0].ssr_type, "c*");
    }

    #[test]
    fn imperfect_replaces_contained_perfect() {
        let fa = b">s\nCCCCAGAGAGAGAGAGTAGAGAGAGCCCC\n";
        let rec = &parse_all(fa)[0];
        let cfg = Config {
            imperfect: Some(ImperfectParams::default()),
            ..Config::default()
        };
        let r = analyze_record(rec, &cfg);
        assert_eq!(r.loci.len(), 1);
        assert_eq!(r.loci[0].ssr_type, "i2");
        assert_eq!(r.loci[0].pattern, "(AG)6t(AG)4");
    }
}

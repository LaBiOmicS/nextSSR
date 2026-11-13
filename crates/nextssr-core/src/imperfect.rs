//! Imperfect (interrupted) SSR detection.
//!
//! Perfect "seeds" of the same motif family (same unit size and same rotation
//! class on the same strand, e.g. `AG`/`GA`) separated by at most
//! [`ImperfectParams::max_gap`] bp are chained greedily from left to right.
//! A chain is accepted as an imperfect SSR when it contains at least two seeds,
//! spans at least `min_repeats × unit` bp and its purity (fraction of the span
//! covered by seeds) is at least [`ImperfectParams::min_purity`]. Chaining
//! seeds that differ in phase naturally absorbs substitutions, insertions and
//! deletions inside the repeat.

use std::collections::HashMap;

use crate::locus::{ImperfectInfo, Ssr};
use crate::motif::min_rotation;
use crate::params::ImperfectParams;
use crate::scan::PerfectHit;

#[derive(Debug, Clone)]
struct Chain {
    seeds: Vec<PerfectHit>,
    start: usize,
    end: usize,
    covered: usize,
}

impl Chain {
    fn new(h: PerfectHit) -> Self {
        Self {
            seeds: vec![h],
            start: h.start,
            end: h.end,
            covered: h.end - h.start,
        }
    }

    /// Coverage gained and new end if `h` is appended.
    fn extension(&self, h: &PerfectHit) -> (usize, usize) {
        let new_end = self.end.max(h.end);
        let add = if h.start >= self.end {
            h.end - h.start
        } else {
            h.end.saturating_sub(self.end)
        };
        (self.covered + add, new_end)
    }
}

fn upper(seq: &[u8], a: usize, b: usize) -> String {
    seq[a..b]
        .iter()
        .map(|c| c.to_ascii_uppercase() as char)
        .collect()
}

fn finish(
    chain: Chain,
    seq: &[u8],
    unit: usize,
    min_len: usize,
    p: &ImperfectParams,
) -> Option<Ssr> {
    if chain.seeds.len() < 2 {
        return None;
    }
    let span = chain.end - chain.start;
    let purity = chain.covered as f64 / span as f64;
    if span < min_len || purity + 1e-12 < p.min_purity {
        return None;
    }
    let mut pattern = String::new();
    let mut prev_end: Option<usize> = None;
    for s in &chain.seeds {
        if let Some(pe) = prev_end {
            if s.start > pe {
                pattern.push_str(&String::from_utf8_lossy(&seq[pe..s.start]).to_ascii_lowercase());
            } else if s.start < pe {
                pattern.push_str(&format!("{{-{}}}", pe - s.start));
            }
        }
        pattern.push_str(&format!(
            "({}){}",
            upper(seq, s.start, s.start + unit),
            s.copies
        ));
        prev_end = Some(s.end);
    }
    let first = chain.seeds[0];
    Some(Ssr {
        start: chain.start + 1,
        end: chain.end,
        unit,
        motif: upper(seq, first.start, first.start + unit),
        repeats: span / unit,
        imperfect: Some(ImperfectInfo {
            purity,
            interruptions: chain.seeds.len() - 1,
            pattern,
        }),
    })
}

/// Chain seeds (sorted by start, all of unit `unit`) into imperfect SSRs.
/// `min_len` is the minimum span (`min_repeats × unit`).
pub fn chain_seeds(
    seq: &[u8],
    seeds: &[PerfectHit],
    unit: usize,
    min_len: usize,
    p: &ImperfectParams,
) -> Vec<Ssr> {
    let mut open: HashMap<String, Chain> = HashMap::new();
    let mut out = Vec::new();
    for &h in seeds {
        let class = min_rotation(&upper(seq, h.start, h.start + unit));
        match open.get_mut(&class) {
            Some(chain) => {
                let gap = h.start as i64 - chain.end as i64;
                let (covered, new_end) = chain.extension(&h);
                let purity = covered as f64 / (new_end - chain.start) as f64;
                if gap <= p.max_gap as i64 && purity + 1e-12 >= p.min_purity {
                    chain.seeds.push(h);
                    chain.covered = covered;
                    chain.end = new_end;
                } else {
                    let done = std::mem::replace(chain, Chain::new(h));
                    out.extend(finish(done, seq, unit, min_len, p));
                }
            }
            None => {
                open.insert(class, Chain::new(h));
            }
        }
    }
    for (_, chain) in open {
        out.extend(finish(chain, seq, unit, min_len, p));
    }
    out.sort_by_key(|s| (s.start, s.end));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scan::{eq_runs, select_maximal};

    fn imperfect(seq: &[u8], unit: usize, min_rep: usize) -> Vec<Ssr> {
        let p = ImperfectParams::default();
        let seed = p.seed_for(unit);
        let runs = eq_runs(seq, unit, (seed - 1) * unit);
        let seeds = select_maximal(seq, &runs, unit, seed);
        chain_seeds(seq, &seeds, unit, min_rep * unit, &p)
    }

    #[test]
    fn substitution_and_indel() {
        // (AG)4 t (AG)4 : substitution-like interruption
        let s = b"CCCCAGAGAGAGTAGAGAGAGCCCC";
        let r = imperfect(s, 2, 6);
        assert_eq!(r.len(), 1);
        let x = &r[0];
        assert_eq!((x.start, x.end), (5, 21));
        let info = x.imperfect.as_ref().unwrap();
        assert_eq!(info.pattern, "(AG)4t(AG)4");
        assert_eq!(info.interruptions, 1);
        // single-base insertion: (AG)4 a (AG)4
        let s = b"CCCCAGAGAGAGAAGAGAGAGCCC";
        let r = imperfect(s, 2, 6);
        assert_eq!(r.len(), 1);
        let info = r[0].imperfect.as_ref().unwrap();
        assert_eq!(info.pattern, "(AG)4a(AG)4");
        assert!((info.purity - 16.0 / 17.0).abs() < 1e-9);
    }

    #[test]
    fn low_purity_rejected() {
        // two short seeds separated by a long gap
        let s = b"CCAGAGAGTTTTTAGAGAGCC";
        assert!(imperfect(s, 2, 6).is_empty());
    }
}

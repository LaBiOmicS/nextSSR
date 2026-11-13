//! Native PCR primer design flanking SSR loci.
//!
//! Defaults follow Primer3 (size 18–27/opt 20, Tm 57–63/opt 60 °C, SantaLucia
//! thermodynamics with 50 mM K⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM oligo,
//! max self/pair complementarity 8 (any) and 3 (3' end), max poly-X 5).
//! Complementarity is scored with ungapped antiparallel alignments
//! (+1 per Watson–Crick pair, −1 otherwise), close to Primer3's legacy
//! `PRIMER_MAX_SELF_ANY` / `PRIMER_MAX_SELF_END` scores.
//!
//! Primers never overlap the target locus, never contain non-ACGT bases and,
//! by default, avoid other SSR loci (so that a marker is not anchored inside a
//! second repeat).

use std::cmp::Ordering;

use serde::{Deserialize, Serialize};

use crate::motif::{revcomp, CODE};
use crate::params::ParamError;
use crate::thermo::{NnPrefix, ThermoParams};

/// Primer design constraints.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PrimerParams {
    pub min_size: usize,
    pub opt_size: usize,
    pub max_size: usize,
    pub min_tm: f64,
    pub opt_tm: f64,
    pub max_tm: f64,
    pub min_gc: f64,
    pub max_gc: f64,
    /// Longest allowed mononucleotide run inside a primer.
    pub max_poly_x: usize,
    /// Number of G/C required among the last 5 bases at the 3' end (0 = off).
    pub gc_clamp: usize,
    /// Maximum Tm difference between the two primers (°C).
    pub max_tm_diff: f64,
    pub product_min: usize,
    pub product_max: usize,
    pub max_self_any: i32,
    pub max_self_end: i32,
    pub max_pair_any: i32,
    pub max_pair_end: i32,
    /// Do not place primers over other SSR loci.
    pub avoid_other_ssrs: bool,
    /// Best candidates kept per side before pairing.
    pub max_candidates: usize,
    pub thermo: ThermoParams,
}

impl Default for PrimerParams {
    fn default() -> Self {
        Self {
            min_size: 18,
            opt_size: 20,
            max_size: 27,
            min_tm: 57.0,
            opt_tm: 60.0,
            max_tm: 63.0,
            min_gc: 30.0,
            max_gc: 70.0,
            max_poly_x: 5,
            gc_clamp: 0,
            max_tm_diff: 3.0,
            product_min: 100,
            product_max: 300,
            max_self_any: 8,
            max_self_end: 3,
            max_pair_any: 8,
            max_pair_end: 3,
            avoid_other_ssrs: true,
            max_candidates: 150,
            thermo: ThermoParams::default(),
        }
    }
}

impl PrimerParams {
    // Negated comparisons are deliberate: they also reject NaN values.
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    pub fn validate(&self) -> Result<(), ParamError> {
        let e = |m: &str| Err(ParamError(m.to_string()));
        if !(self.min_size >= 8
            && self.min_size <= self.opt_size
            && self.opt_size <= self.max_size
            && self.max_size <= 36)
        {
            return e("primer sizes must satisfy 8 <= min <= opt <= max <= 36");
        }
        if !(self.min_tm <= self.opt_tm && self.opt_tm <= self.max_tm) {
            return e("primer Tm must satisfy min <= opt <= max");
        }
        if !(self.min_gc <= self.max_gc) {
            return e("primer GC must satisfy min <= max");
        }
        if self.product_min > self.product_max {
            return e("product size must satisfy min <= max");
        }
        if self.gc_clamp > 5 {
            return e("GC clamp must be <= 5");
        }
        Ok(())
    }
}

/// One designed primer. Coordinates are 1-based, inclusive, on the + strand
/// (for the reverse primer: the binding site whose reverse complement is the primer).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Primer {
    pub seq: String,
    pub start: usize,
    pub end: usize,
    pub tm: f64,
    pub gc: f64,
    pub penalty: f64,
}

/// A primer pair amplifying a locus.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PrimerPair {
    pub left: Primer,
    pub right: Primer,
    pub product_size: usize,
    pub penalty: f64,
}

/// Result of primer design for one locus.
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(tag = "status", rename_all = "lowercase")]
pub enum PrimerOutcome {
    Designed(PrimerPair),
    Failed { reason: String },
}

impl PrimerOutcome {
    fn fail(reason: &str) -> Self {
        PrimerOutcome::Failed {
            reason: reason.to_string(),
        }
    }

    pub fn pair(&self) -> Option<&PrimerPair> {
        match self {
            PrimerOutcome::Designed(p) => Some(p),
            PrimerOutcome::Failed { .. } => None,
        }
    }
}

#[inline]
fn pairs(a: u8, b: u8) -> bool {
    let (x, y) = (CODE[a as usize], CODE[b as usize]);
    x < 4 && y < 4 && x + y == 3
}

/// Max ungapped local complementarity between `x` and `y` (both 5'->3').
pub fn compl_any(x: &[u8], y: &[u8]) -> i32 {
    let (n, m) = (x.len() as i64, y.len() as i64);
    let mut best = 0;
    // x[i] pairs with y[m-1-(i-s)]
    for s in -(m - 1)..n {
        let mut cur = 0;
        for i in s.max(0)..n.min(s + m) {
            let j = (m - 1 - (i - s)) as usize;
            cur += if pairs(x[i as usize], y[j]) { 1 } else { -1 };
            if cur < 0 {
                cur = 0;
            }
            best = best.max(cur);
        }
    }
    best
}

/// Max complementarity of alignments anchored at the 3' end of `x`.
pub fn compl_end(x: &[u8], y: &[u8]) -> i32 {
    let (n, m) = (x.len() as i64, y.len() as i64);
    let mut best = 0;
    let last = n - 1;
    for s in (last - (m - 1))..=last {
        let mut cur = 0;
        let mut i = last;
        while i >= 0 && i - s >= 0 && i - s < m {
            let j = (m - 1 - (i - s)) as usize;
            cur += if pairs(x[i as usize], y[j]) { 1 } else { -1 };
            best = best.max(cur);
            i -= 1;
        }
    }
    best
}

fn max_homopolymer(s: &[u8]) -> usize {
    let mut best = 0;
    let mut run = 0;
    let mut prev = 255u8;
    for &b in s {
        let c = CODE[b as usize];
        if c == prev {
            run += 1;
        } else {
            run = 1;
            prev = c;
        }
        best = best.max(run);
    }
    best
}

#[derive(Debug, Clone)]
struct Cand {
    pos: usize, // 0-based start of the binding site on +
    len: usize,
    seq: Vec<u8>, // primer sequence 5'->3' (uppercase)
    tm: f64,
    gc: f64,
    penalty: f64,
}

fn overlaps(a: usize, b: usize, regions: &[(usize, usize)]) -> bool {
    regions.iter().any(|&(s, e)| a < e && s < b)
}

fn by_penalty(a: &Cand, b: &Cand) -> Ordering {
    a.penalty
        .partial_cmp(&b.penalty)
        .unwrap_or(Ordering::Equal)
        .then(a.pos.cmp(&b.pos))
        .then(a.len.cmp(&b.len))
}

/// Per-locus template window with O(1) composition and Tm queries.
struct Window {
    up: Vec<u8>,
    rc: Vec<u8>,
    bad: Vec<u32>,
    gc: Vec<u32>,
    nn: NnPrefix,
}

impl Window {
    fn new(seq: &[u8]) -> Self {
        let up: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
        let rc = revcomp(&up);
        let mut bad = vec![0u32; up.len() + 1];
        let mut gc = vec![0u32; up.len() + 1];
        for (i, &b) in up.iter().enumerate() {
            let c = CODE[b as usize];
            bad[i + 1] = bad[i] + u32::from(c >= 4);
            gc[i + 1] = gc[i] + u32::from(c == 1 || c == 2);
        }
        let nn = NnPrefix::new(&up);
        Self {
            up,
            rc,
            bad,
            gc,
            nn,
        }
    }

    /// Evaluate the primer whose binding site is `[i, i+len)` of the window;
    /// `primer` is the oligo 5'->3' (the site itself or its reverse complement).
    fn evaluate(
        &self,
        i: usize,
        len: usize,
        primer: &[u8],
        p: &PrimerParams,
    ) -> Option<(f64, f64, f64)> {
        if self.bad[i + len] - self.bad[i] > 0 {
            return None;
        }
        let gc = 100.0 * f64::from(self.gc[i + len] - self.gc[i]) / len as f64;
        if gc < p.min_gc || gc > p.max_gc {
            return None;
        }
        if max_homopolymer(primer) > p.max_poly_x {
            return None;
        }
        if p.gc_clamp > 0 {
            let tail = &primer[len.saturating_sub(5)..];
            let n = tail
                .iter()
                .filter(|&&b| matches!(CODE[b as usize], 1 | 2))
                .count();
            if n < p.gc_clamp {
                return None;
            }
        }
        let t = self.nn.tm(i, len, &p.thermo);
        if t < p.min_tm || t > p.max_tm {
            return None;
        }
        let penalty = (t - p.opt_tm).abs() + (len as f64 - p.opt_size as f64).abs();
        Some((t, gc, penalty))
    }
}

/// Sort by penalty and keep the best `max_candidates` that pass the (costlier)
/// self-complementarity filters.
fn best_candidates(mut c: Vec<Cand>, p: &PrimerParams) -> Vec<Cand> {
    c.sort_by(by_penalty);
    let mut out = Vec::with_capacity(p.max_candidates.min(c.len()));
    for x in c {
        if out.len() >= p.max_candidates {
            break;
        }
        if compl_any(&x.seq, &x.seq) <= p.max_self_any
            && compl_end(&x.seq, &x.seq) <= p.max_self_end
        {
            out.push(x);
        }
    }
    out
}

/// Design the best primer pair around the target `[tstart, tend)` (0-based,
/// half-open). `excluded` lists other loci (0-based, half-open) to avoid.
pub fn design(
    seq: &[u8],
    tstart: usize,
    tend: usize,
    excluded: &[(usize, usize)],
    p: &PrimerParams,
) -> PrimerOutcome {
    let n = seq.len();
    if tend <= tstart || tend > n {
        return PrimerOutcome::fail("invalid_target");
    }
    if (tend - tstart) + 2 * p.min_size > p.product_max {
        return PrimerOutcome::fail("target_too_long");
    }
    let excl: &[(usize, usize)] = if p.avoid_other_ssrs { excluded } else { &[] };
    let wlo = tend.saturating_sub(p.product_max);
    let whi = (tstart + p.product_max).min(n);
    let w = Window::new(&seq[wlo..whi]);
    let wl = w.up.len();
    let (ts, te) = (tstart - wlo, tend - wlo);

    let mut left: Vec<Cand> = Vec::new();
    for len in p.min_size..=p.max_size {
        if ts < len {
            break;
        }
        for i in 0..=(ts - len) {
            if overlaps(wlo + i, wlo + i + len, excl) {
                continue;
            }
            let primer = &w.up[i..i + len];
            if let Some((t, gc, pen)) = w.evaluate(i, len, primer, p) {
                left.push(Cand {
                    pos: wlo + i,
                    len,
                    seq: primer.to_vec(),
                    tm: t,
                    gc,
                    penalty: pen,
                });
            }
        }
    }
    if left.is_empty() {
        return PrimerOutcome::fail("no_left_primer");
    }

    let mut right: Vec<Cand> = Vec::new();
    for len in p.min_size..=p.max_size {
        if wl < te + len {
            break;
        }
        for j in te..=(wl - len) {
            if overlaps(wlo + j, wlo + j + len, excl) {
                continue;
            }
            let primer = &w.rc[wl - j - len..wl - j];
            if let Some((t, gc, pen)) = w.evaluate(j, len, primer, p) {
                right.push(Cand {
                    pos: wlo + j,
                    len,
                    seq: primer.to_vec(),
                    tm: t,
                    gc,
                    penalty: pen,
                });
            }
        }
    }
    if right.is_empty() {
        return PrimerOutcome::fail("no_right_primer");
    }

    let left = best_candidates(left, p);
    let right = best_candidates(right, p);
    if left.is_empty() {
        return PrimerOutcome::fail("no_left_primer");
    }
    if right.is_empty() {
        return PrimerOutcome::fail("no_right_primer");
    }

    let mut best: Option<(f64, &Cand, &Cand, usize)> = None;
    for l in &left {
        if let Some((bp, ..)) = best {
            if l.penalty >= bp {
                break; // sorted: no better pair possible
            }
        }
        for r in &right {
            if let Some((bp, ..)) = best {
                if l.penalty + r.penalty >= bp {
                    break;
                }
            }
            let product = r.pos + r.len - l.pos;
            if product < p.product_min || product > p.product_max {
                continue;
            }
            let dtm = (l.tm - r.tm).abs();
            if dtm > p.max_tm_diff {
                continue;
            }
            let pen = l.penalty + r.penalty + dtm;
            if let Some((bp, ..)) = best {
                if pen >= bp {
                    continue;
                }
            }
            if compl_any(&l.seq, &r.seq) > p.max_pair_any
                || compl_end(&l.seq, &r.seq).max(compl_end(&r.seq, &l.seq)) > p.max_pair_end
            {
                continue;
            }
            best = Some((pen, l, r, product));
        }
    }
    let Some((pen, l, r, product)) = best else {
        return PrimerOutcome::fail("no_compatible_pair");
    };
    let mk = |c: &Cand| Primer {
        seq: String::from_utf8(c.seq.clone()).expect("ASCII"),
        start: c.pos + 1,
        end: c.pos + c.len,
        tm: round2(c.tm),
        gc: round2(c.gc),
        penalty: round3(c.penalty),
    };
    PrimerOutcome::Designed(PrimerPair {
        left: mk(l),
        right: mk(r),
        product_size: product,
        penalty: round3(pen),
    })
}

fn round2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}

fn round3(x: f64) -> f64 {
    (x * 1000.0).round() / 1000.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complementarity() {
        // perfectly self-complementary palindrome
        assert_eq!(compl_any(b"GAATTC", b"GAATTC"), 6);
        assert_eq!(compl_any(b"AAAA", b"AAAA"), 0);
        assert_eq!(compl_any(b"AAAA", b"TTTT"), 4);
        // fully complementary duplex
        assert_eq!(compl_end(b"AAAAGG", b"CCTTTT"), 6);
        // only the 3'-terminal GG of x pairs with the 3'-terminal CC of y
        assert_eq!(compl_end(b"AAAAGG", b"AAAACC"), 2);
        // complementarity away from the 3' end does not count for END
        assert_eq!(compl_end(b"GGAAAA", b"AAAACC"), 0);
        assert_eq!(compl_any(b"GGAAAA", b"AAAACC"), 2);
    }

    fn rnd_seq(len: usize, mut x: u64) -> Vec<u8> {
        (0..len)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x >> 11) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn designs_valid_pair() {
        let mut s = rnd_seq(400, 7);
        let ssr = b"AGAGAGAGAGAGAGAGAGAG";
        s.extend_from_slice(ssr);
        s.extend(rnd_seq(400, 99));
        let p = PrimerParams::default();
        let out = design(&s, 400, 420, &[], &p);
        let pair = out.pair().expect("designed");
        assert!(pair.left.end <= 400);
        assert!(pair.right.start > 420);
        assert!((p.product_min..=p.product_max).contains(&pair.product_size));
        assert_eq!(pair.product_size, pair.right.end - pair.left.start + 1);
        for pr in [&pair.left, &pair.right] {
            assert!(pr.tm >= p.min_tm && pr.tm <= p.max_tm);
            assert!(pr.seq.len() >= p.min_size && pr.seq.len() <= p.max_size);
        }
        // the reverse primer is the reverse complement of its binding site
        let site = &s[pair.right.start - 1..pair.right.end];
        assert_eq!(revcomp(site), pair.right.seq.as_bytes());
        assert_eq!(
            &s[pair.left.start - 1..pair.left.end],
            pair.left.seq.as_bytes()
        );
    }

    #[test]
    fn fails_at_sequence_edge() {
        let s = b"AGAGAGAGAGAGAGAGAGAGCCGATCGTAGCTAGCTAGCTGATCG".to_vec();
        let out = design(&s, 0, 20, &[], &PrimerParams::default());
        assert_eq!(out, PrimerOutcome::fail("no_left_primer"));
    }
}

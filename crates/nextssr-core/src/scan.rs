//! Linear-time perfect repeat scanner.
//!
//! For a unit size `k`, let `eq(i)` be true when `s[i]` and `s[i+k]` are the same
//! nucleotide (A/C/G/T, case-insensitive). A maximal block of consecutive true
//! `eq` positions `[a, b)` ("run") describes the periodic region `[a, b+k)` with
//! period `k`, which holds `1 + (b-a)/k` complete copies of the motif `s[a..a+k]`.
//!
//! * [`select_maximal`] reports every run long enough (nextSSR semantics).
//! * [`select_misa`] replays the left-to-right, non-overlapping, greedy scan of
//!   the `misa.pl` regular expression `(([acgt]{k})\2{m,})`, including the way a
//!   rejected non-primitive match (e.g. `(AAA)5` at k=3) consumes sequence and
//!   masks a following SSR.
//!
//! Both are exact reformulations of the respective algorithms, verified against
//! `misa.pl` in the integration tests.

use crate::motif::{is_primitive, CODE};

/// Sequences at least this long are scanned in parallel chunks.
#[cfg(feature = "parallel")]
const PAR_MIN_LEN: usize = 1 << 22;
#[cfg(feature = "parallel")]
const CHUNK: usize = 1 << 21;

/// A perfect repeat (0-based, half-open coordinates).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PerfectHit {
    pub start: usize,
    pub end: usize,
    pub unit: usize,
    pub copies: usize,
}

/// Runs of `eq` restricted to positions `[lo, hi)`. Runs shorter than `min_len`
/// are dropped unless they touch the range limits (they may continue in a
/// neighbouring chunk).
fn runs_range(seq: &[u8], k: usize, lo: usize, hi: usize, min_len: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut i = lo;
    while i < hi {
        let c = CODE[seq[i] as usize];
        if c >= 4 || c != CODE[seq[i + k] as usize] {
            i += 1;
            continue;
        }
        let a = i;
        i += 1;
        while i < hi {
            let c = CODE[seq[i] as usize];
            if c >= 4 || c != CODE[seq[i + k] as usize] {
                break;
            }
            i += 1;
        }
        let b = i;
        if b - a >= min_len || a == lo || b == hi {
            out.push((a, b));
        }
    }
    out
}

/// All `eq` runs of length `>= min_len` for unit `k` (sorted, disjoint).
pub fn eq_runs(seq: &[u8], k: usize, min_len: usize) -> Vec<(usize, usize)> {
    if k == 0 || seq.len() <= k {
        return Vec::new();
    }
    let n = seq.len() - k; // eq(i) defined for i in [0, n)

    #[cfg(feature = "parallel")]
    if n >= PAR_MIN_LEN {
        use rayon::prelude::*;
        let chunks: Vec<(usize, usize)> = (0..n)
            .step_by(CHUNK)
            .map(|lo| (lo, (lo + CHUNK).min(n)))
            .collect();
        let parts: Vec<Vec<(usize, usize)>> = chunks
            .par_iter()
            .map(|&(lo, hi)| runs_range(seq, k, lo, hi, min_len))
            .collect();
        let mut merged: Vec<(usize, usize)> = Vec::new();
        for part in parts {
            for r in part {
                match merged.last_mut() {
                    Some(last) if last.1 == r.0 => last.1 = r.1,
                    _ => merged.push(r),
                }
            }
        }
        merged.retain(|&(a, b)| b - a >= min_len);
        return merged;
    }

    let mut runs = runs_range(seq, k, 0, n, min_len);
    runs.retain(|&(a, b)| b - a >= min_len);
    runs
}

/// nextSSR semantics: every maximal perfect repeat with at least `min_repeats`
/// copies of a primitive motif.
pub fn select_maximal(
    seq: &[u8],
    runs: &[(usize, usize)],
    k: usize,
    min_repeats: usize,
) -> Vec<PerfectHit> {
    let need = (min_repeats - 1) * k;
    runs.iter()
        .filter(|&&(a, b)| b - a >= need)
        .filter(|&&(a, _)| is_primitive(&seq[a..a + k]))
        .map(|&(a, b)| {
            let copies = 1 + (b - a) / k;
            PerfectHit {
                start: a,
                end: a + copies * k,
                unit: k,
                copies,
            }
        })
        .collect()
}

/// `misa.pl` semantics (greedy, non-overlapping regex scan with masking).
pub fn select_misa(
    seq: &[u8],
    runs: &[(usize, usize)],
    k: usize,
    min_repeats: usize,
) -> Vec<PerfectHit> {
    let need = (min_repeats - 1) * k;
    let mut pos = 0usize; // regex resume position (`pos($seq)`)
    let mut out = Vec::new();
    for &(a, b) in runs {
        let p = a.max(pos);
        if p >= b || b - p < need {
            continue;
        }
        let copies = 1 + (b - p) / k;
        let end = p + copies * k;
        pos = end;
        if is_primitive(&seq[p..p + k]) {
            out.push(PerfectHit {
                start: p,
                end,
                unit: k,
                copies,
            });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn maximal(seq: &[u8], k: usize, r: usize) -> Vec<(usize, usize, usize)> {
        let runs = eq_runs(seq, k, (r - 1) * k);
        select_maximal(seq, &runs, k, r)
            .into_iter()
            .map(|h| (h.start, h.end, h.copies))
            .collect()
    }

    fn misa(seq: &[u8], k: usize, r: usize) -> Vec<(usize, usize, usize)> {
        let runs = eq_runs(seq, k, (r - 1) * k);
        select_misa(seq, &runs, k, r)
            .into_iter()
            .map(|h| (h.start, h.end, h.copies))
            .collect()
    }

    #[test]
    fn simple_dinucleotide() {
        let s = b"GGGACACACACACACTTT";
        assert_eq!(maximal(s, 2, 6), vec![(3, 15, 6)]);
        assert_eq!(misa(s, 2, 6), vec![(3, 15, 6)]);
    }

    #[test]
    fn rejects_non_primitive_and_partial_copies() {
        // AAAAAAAAAAAA at k=2 is (AA)6: non-primitive -> rejected.
        assert!(maximal(b"AAAAAAAAAAAA", 2, 6).is_empty());
        // partial trailing copy is not counted
        assert_eq!(maximal(b"ACACACACACACA", 2, 6), vec![(0, 12, 6)]);
    }

    #[test]
    fn misa_masking_vs_maximal() {
        // (AAG)5 preceded by 15 A: the regex consumes AAA x5 first and misses it.
        let mut s = b"CCCCG".to_vec();
        s.extend(std::iter::repeat(b'A').take(15));
        s.extend_from_slice(b"GAAGAAGAAGAAGCCTCG");
        assert!(misa(&s, 3, 5).is_empty());
        assert_eq!(maximal(&s, 3, 5), vec![(18, 33, 5)]);
    }

    #[test]
    fn case_insensitive_and_n_breaks() {
        assert_eq!(maximal(b"acACacACacAC", 2, 6), vec![(0, 12, 6)]);
        assert!(maximal(b"ACACACNCACACAC", 2, 6).is_empty());
    }

    #[test]
    fn chunked_equals_serial() {
        // Deterministic pseudo-random sequence with planted repeats.
        let mut x: u64 = 42;
        let mut s = Vec::new();
        while s.len() < 300_000 {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            if x % 97 == 0 {
                let m: &[u8] = [b"AG".as_slice(), b"AAT", b"C", b"ACGT"][(x % 4) as usize];
                for _ in 0..(5 + x % 20) {
                    s.extend_from_slice(m);
                }
            } else {
                s.push(b"ACGT"[(x >> 3) as usize % 4]);
            }
        }
        for k in 1..=6 {
            let serial = runs_range(&s, k, 0, s.len() - k, k);
            let serial: Vec<_> = serial.into_iter().filter(|&(a, b)| b - a >= k).collect();
            let mut chunked: Vec<(usize, usize)> = Vec::new();
            let n = s.len() - k;
            let mut lo = 0;
            while lo < n {
                let hi = (lo + 1000).min(n);
                for r in runs_range(&s, k, lo, hi, k) {
                    match chunked.last_mut() {
                        Some(last) if last.1 == r.0 => last.1 = r.1,
                        _ => chunked.push(r),
                    }
                }
                lo = hi;
            }
            chunked.retain(|&(a, b)| b - a >= k);
            assert_eq!(serial, chunked, "k={k}");
        }
    }
}

//! Oligonucleotide melting temperature (nearest-neighbour thermodynamics).
//!
//! Implements the unified parameters of SantaLucia (1998, PNAS 95:1460) with the
//! SantaLucia salt correction and the divalent-cation conversion of
//! von Ahsen et al. (2001), i.e. the default model of Primer3 (≥ 2.0):
//! `PRIMER_TM_FORMULA=1`, `PRIMER_SALT_CORRECTIONS=1`, 50 mM monovalent,
//! 1.5 mM Mg²⁺, 0.6 mM dNTP and 50 nM oligo.

use serde::{Deserialize, Serialize};

use crate::motif::CODE;

/// Reaction conditions.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ThermoParams {
    /// Monovalent cation concentration (mM).
    pub mv_conc: f64,
    /// Divalent cation (Mg²⁺) concentration (mM).
    pub dv_conc: f64,
    /// dNTP concentration (mM).
    pub dntp_conc: f64,
    /// Oligonucleotide concentration (nM).
    pub dna_conc: f64,
}

impl Default for ThermoParams {
    fn default() -> Self {
        Self {
            mv_conc: 50.0,
            dv_conc: 1.5,
            dntp_conc: 0.6,
            dna_conc: 50.0,
        }
    }
}

impl ThermoParams {
    /// Monovalent-equivalent salt concentration (mM): `mv + 120·√(Mg − dNTP)`.
    pub fn salt_equivalent(&self) -> f64 {
        let free_mg = (self.dv_conc - self.dntp_conc).max(0.0);
        self.mv_conc + 120.0 * free_mg.sqrt()
    }
}

/// SantaLucia (1998) unified NN parameters: (ΔH kcal/mol, ΔS cal/K·mol),
/// indexed by `4 * code(x) + code(y)` for the 5'->3' dinucleotide `xy`.
const NN: [(f64, f64); 16] = [
    (-7.9, -22.2),  // AA
    (-8.4, -22.4),  // AC
    (-7.8, -21.0),  // AG
    (-7.2, -20.4),  // AT
    (-8.5, -22.7),  // CA
    (-8.0, -19.9),  // CC
    (-10.6, -27.2), // CG
    (-7.8, -21.0),  // CT
    (-8.2, -22.2),  // GA
    (-9.8, -24.4),  // GC
    (-8.0, -19.9),  // GG
    (-8.4, -22.4),  // GT
    (-7.2, -21.3),  // TA
    (-8.2, -22.2),  // TC
    (-8.5, -22.7),  // TG
    (-7.9, -22.2),  // TT
];

const R: f64 = 1.987;
const T_KELVIN: f64 = 273.15;

fn is_self_complementary(codes: &[u8]) -> bool {
    let n = codes.len();
    (0..n).all(|i| codes[i] == 3 - codes[n - 1 - i])
}

/// Melting temperature (°C), or `None` for sequences shorter than 2 nt or
/// containing non-ACGT characters.
pub fn tm(seq: &[u8], p: &ThermoParams) -> Option<f64> {
    if seq.len() < 2 {
        return None;
    }
    let codes: Vec<u8> = seq.iter().map(|&b| CODE[b as usize]).collect();
    if codes.iter().any(|&c| c >= 4) {
        return None;
    }
    let mut dh = 0.0;
    let mut ds = 0.0;
    for w in codes.windows(2) {
        let (h, s) = NN[(4 * w[0] + w[1]) as usize];
        dh += h;
        ds += s;
    }
    // Initiation with terminal G·C or A·T pairs.
    for &end in [codes[0], codes[codes.len() - 1]].iter() {
        if end == 1 || end == 2 {
            dh += 0.1;
            ds += -2.8;
        } else {
            dh += 2.3;
            ds += 4.1;
        }
    }
    let ct = p.dna_conc * 1e-9;
    let conc_term = if is_self_complementary(&codes) {
        ds += -1.4;
        ct
    } else {
        ct / 4.0
    };
    // SantaLucia (1998) entropy salt correction.
    ds += 0.368 * (codes.len() as f64 - 1.0) * (p.salt_equivalent() / 1000.0).ln();
    Some(dh * 1000.0 / (ds + R * conc_term.ln()) - T_KELVIN)
}

/// Prefix sums of nearest-neighbour ΔH/ΔS over a template, giving the Tm of
/// any window in O(1). Because the NN table is strand-symmetric, the Tm of a
/// window equals the Tm of its reverse complement, so the same prefix serves
/// forward and reverse primers.
pub struct NnPrefix {
    codes: Vec<u8>,
    h: Vec<f64>,
    s: Vec<f64>,
}

impl NnPrefix {
    pub fn new(seq: &[u8]) -> Self {
        let codes: Vec<u8> = seq.iter().map(|&b| CODE[b as usize]).collect();
        let mut h = vec![0.0; codes.len().max(1)];
        let mut s = vec![0.0; codes.len().max(1)];
        // h[i] = sum of stacks (k, k+1) for k < i
        for i in 1..codes.len() {
            let (x, y) = (codes[i - 1], codes[i]);
            let (dh, ds) = if x < 4 && y < 4 {
                NN[(4 * x + y) as usize]
            } else {
                (0.0, 0.0)
            };
            h[i] = h[i - 1] + dh;
            s[i] = s[i - 1] + ds;
        }
        Self { codes, h, s }
    }

    /// Tm of `seq[i..i+len]` (caller guarantees ACGT-only content and len >= 2).
    pub fn tm(&self, i: usize, len: usize, p: &ThermoParams) -> f64 {
        let j = i + len - 1;
        let mut dh = self.h[j] - self.h[i];
        let mut ds = self.s[j] - self.s[i];
        for &end in [self.codes[i], self.codes[j]].iter() {
            if end == 1 || end == 2 {
                dh += 0.1;
                ds += -2.8;
            } else {
                dh += 2.3;
                ds += 4.1;
            }
        }
        let ct = p.dna_conc * 1e-9;
        let conc_term = if is_self_complementary(&self.codes[i..=j]) {
            ds += -1.4;
            ct
        } else {
            ct / 4.0
        };
        ds += 0.368 * (len as f64 - 1.0) * (p.salt_equivalent() / 1000.0).ln();
        dh * 1000.0 / (ds + R * conc_term.ln()) - T_KELVIN
    }
}

/// GC content (%) of an ACGT sequence.
pub fn gc_percent(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|&&b| matches!(CODE[b as usize], 1 | 2))
        .count();
    100.0 * gc as f64 / seq.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn salt_equivalent_default() {
        let p = ThermoParams::default();
        assert!((p.salt_equivalent() - 163.8).abs() < 0.1);
    }

    #[test]
    fn known_values() {
        // Reference values from an independent implementation of the same model
        // (SantaLucia 1998 + entropy salt correction), 50 mM Na+, 250 nM oligo.
        let p = ThermoParams {
            mv_conc: 50.0,
            dv_conc: 0.0,
            dntp_conc: 0.0,
            dna_conc: 250.0,
        };
        for (s, expected) in [
            (&b"AGCTGCATGCAGTCAGTCAG"[..], 57.10),
            (b"GATTACAAGCTACGGATC", 46.47),
            (b"ATATATTTAATATAATAAAT", 30.29),
        ] {
            let t = tm(s, &p).unwrap();
            assert!((t - expected).abs() < 0.01, "{t} vs {expected}");
        }
        // default (Primer3-like) conditions put a typical 20-mer near 60 °C
        let t = tm(b"AGCTGCATGCAGTCAGTCAG", &ThermoParams::default()).unwrap();
        assert!((55.0..65.0).contains(&t), "{t}");
        assert!(tm(b"ACGN", &p).is_none());
    }

    #[test]
    fn prefix_matches_direct_and_is_strand_symmetric() {
        let seq = b"GATTACAAGCTACGGATCAGCTGCATGCAGTCAGTCAGGAATTCTT";
        let pre = NnPrefix::new(seq);
        let p = ThermoParams::default();
        for i in 0..seq.len() - 18 {
            for len in 18..=(seq.len() - i).min(27) {
                let w = &seq[i..i + len];
                let direct = tm(w, &p).unwrap();
                assert!((pre.tm(i, len, &p) - direct).abs() < 1e-9);
                let rc = crate::motif::revcomp(w);
                assert!((tm(&rc, &p).unwrap() - direct).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn gc() {
        assert_eq!(gc_percent(b"GGCC"), 100.0);
        assert_eq!(gc_percent(b"ATGC"), 50.0);
    }
}

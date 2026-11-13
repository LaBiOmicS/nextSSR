//! SSRs, loci (single or compound SSRs) and their textual notation.

use serde::Serialize;

use crate::params::Mode;
use crate::primer::PrimerOutcome;

/// Extra information for imperfect SSRs.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ImperfectInfo {
    /// Fraction of the span covered by perfect seeds.
    pub purity: f64,
    /// Number of interruptions (seeds - 1).
    pub interruptions: usize,
    /// Notation such as `(AG)4t(AG)5`.
    pub pattern: String,
}

/// A single (perfect or imperfect) SSR. Coordinates are 1-based, inclusive.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Ssr {
    pub start: usize,
    pub end: usize,
    pub unit: usize,
    /// Motif in the phase found in the sequence (uppercase).
    pub motif: String,
    /// Number of complete motif copies (imperfect: span / unit).
    pub repeats: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub imperfect: Option<ImperfectInfo>,
}

impl Ssr {
    pub fn len(&self) -> usize {
        self.end + 1 - self.start
    }

    pub fn is_empty(&self) -> bool {
        false
    }

    /// `(AG)8`, or the imperfect pattern.
    pub fn notation(&self) -> String {
        match &self.imperfect {
            Some(i) => i.pattern.clone(),
            None => format!("({}){}", self.motif, self.repeats),
        }
    }
}

/// A reported locus: one SSR, or several SSRs forming a compound SSR.
#[derive(Debug, Clone, Serialize)]
pub struct Locus {
    /// 1-based number of the locus within its sequence.
    pub nr: usize,
    /// `p1`..`p6`, `i1`..`i6` (imperfect), `c` or `c*` (compound).
    pub ssr_type: String,
    /// Notation, e.g. `(AG)8`, `(A)10ggtc(AG)7`.
    pub pattern: String,
    pub start: usize,
    pub end: usize,
    pub members: Vec<Ssr>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub primers: Option<PrimerOutcome>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub flank5: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub flank3: Option<String>,
}

impl Locus {
    pub fn is_compound(&self) -> bool {
        self.members.len() > 1
    }

    /// Size as reported by MISA: `end - start + 1`.
    pub fn size(&self) -> i64 {
        self.end as i64 - self.start as i64 + 1
    }
}

/// Perl `lc substr($seq, $off, $len)` including the negative-length semantics
/// (a negative length leaves that many characters off the end of the string).
fn perl_substr_lc(seq: &[u8], off: usize, len: i64) -> String {
    let n = seq.len();
    if off > n {
        return String::new();
    }
    let end = if len >= 0 {
        (off + len as usize).min(n)
    } else {
        let cut = (-len) as usize;
        if cut > n {
            return String::new();
        }
        n - cut
    };
    if end <= off {
        return String::new();
    }
    String::from_utf8_lossy(&seq[off..end]).to_ascii_lowercase()
}

fn single_type(s: &Ssr) -> String {
    match s.imperfect {
        Some(_) => format!("i{}", s.unit),
        None => format!("p{}", s.unit),
    }
}

/// Group sorted SSRs into loci.
///
/// Two consecutive SSRs belong to the same compound when the number of bases
/// between them (`start2 - end1 - 1`, negative for overlaps) is at most
/// `max_interruption`, exactly as in `misa.pl` (the comparison is always made
/// with the previous member).
pub fn build_loci(
    ssrs: Vec<Ssr>,
    seq: &[u8],
    max_interruption: Option<i64>,
    mode: Mode,
) -> Vec<Locus> {
    let mut loci = Vec::new();
    let mut groups: Vec<Vec<Ssr>> = Vec::new();
    for s in ssrs {
        let join = match (max_interruption, groups.last().and_then(|g| g.last())) {
            (Some(d), Some(prev)) => (s.start as i64 - prev.end as i64 - 1) <= d,
            _ => false,
        };
        if join {
            groups.last_mut().expect("non-empty").push(s);
        } else {
            groups.push(vec![s]);
        }
    }
    for (i, members) in groups.into_iter().enumerate() {
        let (ssr_type, pattern) = match mode {
            Mode::Misa => misa_notation(&members, seq),
            Mode::Native => native_notation(&members, seq),
        };
        let start = members[0].start;
        let end = members.last().expect("non-empty").end;
        loci.push(Locus {
            nr: i + 1,
            ssr_type,
            pattern,
            start,
            end,
            members,
            primers: None,
            flank5: None,
            flank3: None,
        });
    }
    loci
}

/// Notation exactly as written by `misa.pl` v2.0 (non-GFF mode).
fn misa_notation(m: &[Ssr], seq: &[u8]) -> (String, String) {
    let one = |s: &Ssr| format!("({}){}", s.motif, s.repeats);
    if m.len() == 1 {
        return (format!("p{}", m[0].unit), one(&m[0]));
    }
    let diff = m[1].start as i64 - m[0].end as i64;
    let (mut ty, mut pat) = if diff == 1 {
        ("c*".to_string(), format!("{}{}*", one(&m[0]), one(&m[1])))
    } else {
        let inter = perl_substr_lc(seq, m[0].end, diff - 1);
        (
            "c".to_string(),
            format!("{}{}{}", one(&m[0]), inter, one(&m[1])),
        )
    };
    for w in m[1..].windows(2) {
        let diff = w[1].start as i64 - w[0].end as i64;
        if diff == 1 {
            pat.push_str(&format!("{}*", one(&w[1])));
            ty = "c*".to_string();
        } else {
            pat.push_str(&perl_substr_lc(seq, w[0].end, diff - 1));
            pat.push_str(&one(&w[1]));
        }
    }
    (ty, pat)
}

/// nextSSR notation: interruptions in lowercase, nothing for adjacent SSRs,
/// `{-n}` for SSRs overlapping by `n` bp; imperfect members in brackets.
fn native_notation(m: &[Ssr], seq: &[u8]) -> (String, String) {
    if m.len() == 1 {
        return (single_type(&m[0]), m[0].notation());
    }
    let member = |s: &Ssr| match s.imperfect {
        Some(_) => format!("[{}]", s.notation()),
        None => s.notation(),
    };
    let mut pat = member(&m[0]);
    let mut contiguous = false;
    for w in m.windows(2) {
        let gap = w[1].start as i64 - w[0].end as i64 - 1;
        if gap > 0 {
            pat.push_str(
                &String::from_utf8_lossy(&seq[w[0].end..w[1].start - 1]).to_ascii_lowercase(),
            );
        } else {
            contiguous = true;
            if gap < 0 {
                pat.push_str(&format!("{{{gap}}}"));
            }
        }
        pat.push_str(&member(&w[1]));
    }
    (if contiguous { "c*" } else { "c" }.to_string(), pat)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ssr(start: usize, end: usize, motif: &str) -> Ssr {
        Ssr {
            start,
            end,
            unit: motif.len(),
            motif: motif.into(),
            repeats: (end + 1 - start) / motif.len(),
            imperfect: None,
        }
    }

    #[test]
    fn perl_substr() {
        let s = b"ABCDEFGH";
        assert_eq!(perl_substr_lc(s, 2, 3), "cde");
        assert_eq!(perl_substr_lc(s, 2, -1), "cdefg");
        assert_eq!(perl_substr_lc(s, 6, -3), "");
        assert_eq!(perl_substr_lc(s, 2, 0), "");
    }

    #[test]
    fn compound_notation() {
        //            1234567890123456789012345678
        let seq = b"ACACACACACACggtcTGTGTGTGTGTG";
        let a = ssr(1, 12, "AC");
        let b = ssr(17, 28, "TG");
        let l = build_loci(vec![a.clone(), b.clone()], seq, Some(100), Mode::Misa);
        assert_eq!(l.len(), 1);
        assert_eq!(l[0].ssr_type, "c");
        assert_eq!(l[0].pattern, "(AC)6ggtc(TG)6");
        let l = build_loci(vec![a.clone(), b.clone()], seq, Some(3), Mode::Misa);
        assert_eq!(l.len(), 2);
        let n = build_loci(vec![a, b], seq, Some(100), Mode::Native);
        assert_eq!(n[0].pattern, "(AC)6ggtc(TG)6");
    }

    #[test]
    fn adjacent_and_overlap() {
        let seq = b"AAAAAAAAAAAAAGAGAGAGAGAGAGAG";
        let a = ssr(1, 13, "A");
        let b = ssr(13, 28, "AG");
        // MISA: overlap -> 'c' with Perl negative-length substr garbage (empty here)
        let m = build_loci(vec![a.clone(), b.clone()], seq, Some(100), Mode::Misa);
        assert_eq!(m[0].ssr_type, "c");
        assert_eq!(
            m[0].pattern,
            "(A)13agagagagagagaga(AG)8".replace("agagagagagagaga", &perl_substr_lc(seq, 13, -1))
        );
        let n = build_loci(vec![a, b], seq, Some(100), Mode::Native);
        assert_eq!(n[0].ssr_type, "c*");
        assert_eq!(n[0].pattern, "(A)13{-1}(AG)8");
        let a = ssr(1, 12, "A");
        let b = ssr(13, 28, "AG");
        let m = build_loci(vec![a, b], seq, Some(0), Mode::Misa);
        assert_eq!(
            (m[0].ssr_type.as_str(), m[0].pattern.as_str()),
            ("c*", "(A)12(AG)8*")
        );
    }
}

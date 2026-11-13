//! In-silico PCR: find every product that a set of primer pairs can
//! amplify from a set of template sequences.
//!
//! * Both strands are searched and every primer can act as forward or reverse,
//!   so single-primer products (F/F, R/R) and swapped products (R/F) are found
//!   as well as the expected F/R product.
//! * A binding site may have up to `max_mismatches` mismatches, of which at
//!   most `max_three_prime_mismatches` within the `three_prime_len` bases at the
//!   primer's 3' end. Primers may contain IUPAC codes; template bases other than
//!   A/C/G/T (e.g. `N`) never match. Lowercase (soft-masked) bases are searched.
//! * Sites are found with an exact pigeonhole filter: a primer with at most `m`
//!   mismatches contains at least one of `m + 1` disjoint segments without
//!   mismatches, so seeding on those segments is lossless.
//!
//! The engine is streaming: records are scanned one at a time (each long record
//! in parallel chunks) and the results are merged in input order.

use std::fmt::Write;

use serde::{Deserialize, Serialize};

use crate::motif::CODE;
use crate::params::ParamError;

/// in-silico PCR search parameters.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct EpcrParams {
    /// Maximum mismatches per primer binding site.
    pub max_mismatches: usize,
    /// Length of the 3'-end region where mismatches are restricted.
    pub three_prime_len: usize,
    /// Maximum mismatches inside the 3'-end region.
    pub max_three_prime_mismatches: usize,
    /// Minimum product size (bp; never smaller than the longer primer).
    pub min_product: usize,
    /// Maximum product size (bp).
    pub max_product: usize,
    /// Amplicons listed per primer pair (all are counted).
    pub max_amplicons: usize,
    /// Binding sites kept per primer and sequence; beyond this the pair is
    /// flagged as saturated and its counts are lower bounds.
    pub max_sites: usize,
    /// Store the amplicon sequences.
    pub with_sequence: bool,
}

impl Default for EpcrParams {
    fn default() -> Self {
        Self {
            max_mismatches: 2,
            three_prime_len: 5,
            max_three_prime_mismatches: 0,
            min_product: 0,
            max_product: 2000,
            max_amplicons: 100,
            max_sites: 100_000,
            with_sequence: false,
        }
    }
}

impl EpcrParams {
    pub fn validate(&self) -> Result<(), ParamError> {
        let e = |m: &str| Err(ParamError(m.to_string()));
        if self.max_mismatches > 10 {
            return e("in-silico PCR mismatches must be <= 10");
        }
        if self.max_three_prime_mismatches > self.max_mismatches {
            return e("in-silico PCR 3' mismatches cannot exceed the total mismatches");
        }
        if self.max_product == 0 || self.max_product > 1_000_000 {
            return e("in-silico PCR maximum product size must be in 1..=1000000");
        }
        if self.min_product > self.max_product {
            return e("in-silico PCR product size must satisfy min <= max");
        }
        if self.max_sites == 0 {
            return e("in-silico PCR maximum sites must be >= 1");
        }
        Ok(())
    }
}

/// The product a pair is designed to amplify (1-based, inclusive).
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Expected {
    pub seq_id: String,
    pub start: usize,
    pub end: usize,
}

/// A primer pair to test.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PairSpec {
    pub id: String,
    /// Forward (left) primer, 5'->3'.
    pub forward: String,
    /// Reverse (right) primer, 5'->3'.
    pub reverse: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expected: Option<Expected>,
}

/// One primer binding site of an amplicon (1-based, inclusive, + strand).
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Binding {
    /// `F` (forward primer of the pair) or `R` (reverse primer).
    pub primer: char,
    pub start: usize,
    pub end: usize,
    pub mismatches: u8,
}

/// A predicted PCR product.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Amplicon {
    pub seq_id: String,
    pub start: usize,
    pub end: usize,
    pub size: usize,
    /// Primer extending on the + strand (5' end of the product).
    pub fwd: Binding,
    /// Primer extending on the − strand (3' end of the product).
    pub rev: Binding,
    /// The designed product of the pair.
    pub expected: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sequence: Option<String>,
}

/// Specificity verdict for a pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Status {
    /// Exactly one product (the expected one, when known).
    Specific,
    /// More than one product.
    Nonspecific,
    /// No product.
    NoProduct,
    /// A single product that is not the expected one.
    Unexpected,
}

impl Status {
    pub fn name(self) -> &'static str {
        match self {
            Status::Specific => "specific",
            Status::Nonspecific => "nonspecific",
            Status::NoProduct => "no_product",
            Status::Unexpected => "unexpected",
        }
    }
}

/// Accumulated result of one pair.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct PairResult {
    /// Total number of products.
    pub amplicons: u64,
    /// Products with no mismatch in either primer.
    pub perfect: u64,
    /// Whether the expected product was found.
    pub expected_found: bool,
    /// Some primer had more than `max_sites` sites in a sequence.
    pub saturated: bool,
    /// The first `max_amplicons` products in input order.
    #[serde(skip)]
    pub list: Vec<Amplicon>,
}

impl PairResult {
    pub fn status(&self, has_expected: bool) -> Status {
        match self.amplicons {
            0 => Status::NoProduct,
            1 if has_expected && !self.expected_found => Status::Unexpected,
            1 => Status::Specific,
            _ => Status::Nonspecific,
        }
    }
}

/// Hits of one record, sparse over pairs.
#[derive(Debug, Clone, Default)]
pub struct RecordHits {
    pub pairs: Vec<(usize, PairResult)>,
}

/// Template base as a 4-bit mask (A=1, C=2, G=4, T=8; anything else 0).
static TMASK: [u8; 256] = {
    let mut t = [0u8; 256];
    t[b'A' as usize] = 1;
    t[b'a' as usize] = 1;
    t[b'C' as usize] = 2;
    t[b'c' as usize] = 2;
    t[b'G' as usize] = 4;
    t[b'g' as usize] = 4;
    t[b'T' as usize] = 8;
    t[b't' as usize] = 8;
    t
};

/// IUPAC primer base as a 4-bit mask (0 = invalid).
fn iupac_mask(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'M' => 3,
        b'R' => 5,
        b'W' => 9,
        b'S' => 6,
        b'Y' => 10,
        b'K' => 12,
        b'V' => 7,
        b'H' => 11,
        b'D' => 13,
        b'B' => 14,
        b'N' => 15,
        _ => 0,
    }
}

fn comp_mask(m: u8) -> u8 {
    ((m & 1) << 3) | ((m & 2) << 1) | ((m & 4) >> 1) | ((m & 8) >> 3)
}

/// Largest seed length (the seed table has 4^K + 1 entries).
const K_MAX: usize = 10;
/// Most exact k-mers a single degenerate seed may expand to.
const MAX_EXPANSION: usize = 4096;
/// Records longer than this are scanned in parallel chunks.
const CHUNK: usize = 1 << 20;

/// A primer in one orientation on the + strand.
#[derive(Debug, Clone)]
struct Pattern {
    masks: Vec<u8>,
    pair: usize,
    /// 0 = forward primer of the pair, 1 = reverse primer.
    role: u8,
    /// `true` when the pattern is the reverse complement of the primer (the
    /// primer extends leftwards and its 3' end is the pattern's first base).
    minus: bool,
    /// Segment length for the pigeonhole filter.
    seg: usize,
}

#[derive(Debug, Clone, Copy)]
struct Entry {
    pat: u32,
    seg: u16,
    off: u16,
}

#[derive(Debug, Clone, Copy)]
struct Site {
    pat: u32,
    start: usize,
    mm: u8,
}

/// Compiled primer set, ready to scan sequences.
#[derive(Debug, Clone)]
pub struct Epcr {
    params: EpcrParams,
    pairs: Vec<PairSpec>,
    patterns: Vec<Pattern>,
    k: usize,
    offsets: Vec<u32>,
    entries: Vec<Entry>,
    max_len: usize,
}

fn clean_primer(s: &str) -> String {
    s.chars()
        .filter(|c| !c.is_whitespace())
        .collect::<String>()
        .to_ascii_uppercase()
}

impl Epcr {
    /// Compile a primer set.
    pub fn new(pairs: Vec<PairSpec>, params: EpcrParams) -> Result<Self, ParamError> {
        params.validate()?;
        let m1 = params.max_mismatches + 1;
        let mut patterns = Vec::with_capacity(pairs.len() * 4);
        let mut pairs = pairs;
        for (pi, p) in pairs.iter_mut().enumerate() {
            p.forward = clean_primer(&p.forward);
            p.reverse = clean_primer(&p.reverse);
            for (role, seq) in [(0u8, &p.forward), (1u8, &p.reverse)] {
                let masks: Vec<u8> = seq.bytes().map(iupac_mask).collect();
                if masks.contains(&0) || masks.is_empty() {
                    return Err(ParamError(format!(
                        "pair {:?}: invalid primer sequence {seq:?}",
                        p.id
                    )));
                }
                if masks.len() < 8 || masks.len() > 255 {
                    return Err(ParamError(format!(
                        "pair {:?}: primer length must be 8..=255 (got {})",
                        p.id,
                        masks.len()
                    )));
                }
                if masks.len() < m1 * 2 {
                    return Err(ParamError(format!(
                        "pair {:?}: primer {seq:?} is too short for {} mismatches",
                        p.id, params.max_mismatches
                    )));
                }
                let seg = masks.len() / m1;
                let rc: Vec<u8> = masks.iter().rev().map(|&m| comp_mask(m)).collect();
                patterns.push(Pattern {
                    masks,
                    pair: pi,
                    role,
                    minus: false,
                    seg,
                });
                patterns.push(Pattern {
                    masks: rc,
                    pair: pi,
                    role,
                    minus: true,
                    seg,
                });
            }
        }
        let k = patterns
            .iter()
            .map(|p| p.seg)
            .min()
            .unwrap_or(K_MAX)
            .min(K_MAX);
        let max_len = patterns.iter().map(|p| p.masks.len()).max().unwrap_or(0);

        // Seed table (CSR): every exact k-mer compatible with the first k bases
        // of each segment of each pattern.
        let nkeys = 1usize << (2 * k);
        let mut raw: Vec<(u32, Entry)> = Vec::new();
        for (pi, pat) in patterns.iter().enumerate() {
            for j in 0..m1 {
                let off = j * pat.seg;
                let prefix = &pat.masks[off..off + k];
                let n: usize = prefix.iter().map(|m| m.count_ones() as usize).product();
                if n > MAX_EXPANSION {
                    return Err(ParamError(format!(
                        "pair {:?}: primer too degenerate for the in-silico PCR seed filter",
                        pairs[pat.pair].id
                    )));
                }
                let mut keys = vec![0u32];
                for &m in prefix {
                    let mut next = Vec::with_capacity(keys.len() * m.count_ones() as usize);
                    for &key in &keys {
                        for c in 0..4u32 {
                            if m & (1 << c) != 0 {
                                next.push((key << 2) | c);
                            }
                        }
                    }
                    keys = next;
                }
                for key in keys {
                    raw.push((
                        key,
                        Entry {
                            pat: pi as u32,
                            seg: j as u16,
                            off: off as u16,
                        },
                    ));
                }
            }
        }
        let mut offsets = vec![0u32; nkeys + 1];
        for (key, _) in &raw {
            offsets[*key as usize + 1] += 1;
        }
        for i in 0..nkeys {
            offsets[i + 1] += offsets[i];
        }
        let mut fill = offsets.clone();
        let mut entries = vec![
            Entry {
                pat: 0,
                seg: 0,
                off: 0
            };
            raw.len()
        ];
        for (key, e) in raw {
            let slot = &mut fill[key as usize];
            entries[*slot as usize] = e;
            *slot += 1;
        }
        Ok(Self {
            params,
            pairs,
            patterns,
            k,
            offsets,
            entries,
            max_len,
        })
    }

    pub fn pairs(&self) -> &[PairSpec] {
        &self.pairs
    }

    pub fn params(&self) -> &EpcrParams {
        &self.params
    }

    /// Seed length used by the filter.
    pub fn seed_len(&self) -> usize {
        self.k
    }

    /// Mismatches of `pat` at `s`, or `None` if over the limits or if an
    /// earlier seed of the pattern also matches exactly (that seed reports it).
    #[inline]
    fn verify(&self, seq: &[u8], s: usize, pat: &Pattern, seg_j: usize) -> Option<u8> {
        let k = self.k;
        for i in 0..seg_j {
            let o = s + i * pat.seg;
            if (0..k).all(|t| TMASK[seq[o + t] as usize] & pat.masks[i * pat.seg + t] != 0) {
                return None;
            }
        }
        let len = pat.masks.len();
        let t3 = self.params.three_prime_len.min(len);
        let (lo3, hi3) = if pat.minus { (0, t3) } else { (len - t3, len) };
        let (m, m3) = (
            self.params.max_mismatches,
            self.params.max_three_prime_mismatches,
        );
        let mut mm = 0;
        let mut mm3 = 0;
        for (i, &pm) in pat.masks.iter().enumerate() {
            if TMASK[seq[s + i] as usize] & pm == 0 {
                mm += 1;
                if mm > m {
                    return None;
                }
                if i >= lo3 && i < hi3 {
                    mm3 += 1;
                    if mm3 > m3 {
                        return None;
                    }
                }
            }
        }
        Some(mm as u8)
    }

    /// Sites starting in `[cs, ce)`.
    fn scan_range(&self, seq: &[u8], cs: usize, ce: usize, out: &mut Vec<Site>) {
        let n = seq.len();
        let k = self.k;
        let kmask = (1usize << (2 * k)) - 1;
        let qend = (ce + self.max_len).min(n);
        let mut key = 0usize;
        let mut valid = 0usize;
        for pos in cs..qend {
            let c = CODE[seq[pos] as usize];
            if c >= 4 {
                valid = 0;
                continue;
            }
            key = ((key << 2) | c as usize) & kmask;
            valid += 1;
            if valid < k {
                continue;
            }
            let q = pos + 1 - k;
            let (a, b) = (self.offsets[key] as usize, self.offsets[key + 1] as usize);
            for e in &self.entries[a..b] {
                let off = e.off as usize;
                if q < off {
                    continue;
                }
                let s = q - off;
                let pat = &self.patterns[e.pat as usize];
                if s < cs || s >= ce || s + pat.masks.len() > n {
                    continue;
                }
                if let Some(mm) = self.verify(seq, s, pat, e.seg as usize) {
                    out.push(Site {
                        pat: e.pat,
                        start: s,
                        mm,
                    });
                }
            }
        }
    }

    fn sites(&self, seq: &[u8]) -> Vec<Site> {
        let n = seq.len();
        if self.patterns.is_empty() || n < self.k {
            return Vec::new();
        }
        #[cfg(feature = "parallel")]
        {
            if n > CHUNK {
                use rayon::prelude::*;
                let chunks: Vec<Vec<Site>> = (0..n.div_ceil(CHUNK))
                    .into_par_iter()
                    .map(|c| {
                        let mut v = Vec::new();
                        self.scan_range(seq, c * CHUNK, ((c + 1) * CHUNK).min(n), &mut v);
                        v
                    })
                    .collect();
                return chunks.into_iter().flatten().collect();
            }
        }
        let mut v = Vec::new();
        for c in 0..n.div_ceil(CHUNK) {
            self.scan_range(seq, c * CHUNK, ((c + 1) * CHUNK).min(n), &mut v);
        }
        v
    }

    /// Scan one template sequence.
    pub fn scan(&self, seq_id: &str, seq: &[u8]) -> RecordHits {
        let mut sites = self.sites(seq);
        // group by pair, then strand, then position
        sites.sort_unstable_by_key(|s| {
            let p = &self.patterns[s.pat as usize];
            (p.pair, p.minus, s.start, p.role)
        });
        let mut out = RecordHits::default();
        let mut i = 0;
        while i < sites.len() {
            let pair = self.patterns[sites[i].pat as usize].pair;
            let mut j = i;
            while j < sites.len() && self.patterns[sites[j].pat as usize].pair == pair {
                j += 1;
            }
            let r = self.pair_amplicons(seq_id, seq, pair, &sites[i..j]);
            if r.amplicons > 0 || r.saturated {
                out.pairs.push((pair, r));
            }
            i = j;
        }
        out
    }

    fn pair_amplicons(&self, seq_id: &str, seq: &[u8], pair: usize, sites: &[Site]) -> PairResult {
        let mut res = PairResult::default();
        let len = |s: &Site| self.patterns[s.pat as usize].masks.len();
        let role = |s: &Site| {
            if self.patterns[s.pat as usize].role == 0 {
                'F'
            } else {
                'R'
            }
        };
        let (mut plus, mut minus): (Vec<Site>, Vec<Site>) = sites
            .iter()
            .partition(|s| !self.patterns[s.pat as usize].minus);
        for v in [&mut plus, &mut minus] {
            for r in 0..2u8 {
                let n = v
                    .iter()
                    .filter(|s| self.patterns[s.pat as usize].role == r)
                    .count();
                if n > self.params.max_sites {
                    res.saturated = true;
                    let mut kept = 0;
                    v.retain(|s| {
                        if self.patterns[s.pat as usize].role != r {
                            return true;
                        }
                        kept += 1;
                        kept <= self.params.max_sites
                    });
                }
            }
        }
        // minus sites sorted by their end (the 3' end of the product)
        let end = |s: &Site| s.start + len(s) - 1;
        minus.sort_by_key(|s| (end(s), s.start, self.patterns[s.pat as usize].role));
        let ends: Vec<usize> = minus.iter().map(end).collect();
        let mut perfect_prefix = vec![0u64; minus.len() + 1];
        for (i, s) in minus.iter().enumerate() {
            perfect_prefix[i + 1] = perfect_prefix[i] + u64::from(s.mm == 0);
        }
        let spec = &self.pairs[pair];
        let longest = spec.forward.len().max(spec.reverse.len());
        let pmin = self.params.min_product.max(longest);
        let pmax = self.params.max_product;
        if pmin > pmax {
            return res;
        }
        let exp = spec
            .expected
            .as_ref()
            .filter(|e| e.seq_id == seq_id)
            .map(|e| (e.start - 1, e.end - 1));
        for a in &plus {
            let lo = ends.partition_point(|&e| e < a.start + pmin - 1);
            let hi = ends.partition_point(|&e| e < a.start + pmax);
            if lo >= hi {
                continue;
            }
            res.amplicons += (hi - lo) as u64;
            if a.mm == 0 {
                res.perfect += perfect_prefix[hi] - perfect_prefix[lo];
            }
            if let Some((es, ee)) = exp {
                if a.start == es {
                    let x = lo + ends[lo..hi].partition_point(|&e| e < ee);
                    if x < hi && ends[x] == ee {
                        res.expected_found = true;
                    }
                }
            }
            for b in &minus[lo..hi] {
                if res.list.len() >= self.params.max_amplicons {
                    break;
                }
                let (s, e) = (a.start, end(b));
                res.list.push(Amplicon {
                    seq_id: seq_id.to_string(),
                    start: s + 1,
                    end: e + 1,
                    size: e + 1 - s,
                    fwd: Binding {
                        primer: role(a),
                        start: a.start + 1,
                        end: a.start + len(a),
                        mismatches: a.mm,
                    },
                    rev: Binding {
                        primer: role(b),
                        start: b.start + 1,
                        end: e + 1,
                        mismatches: b.mm,
                    },
                    expected: exp == Some((s, e)),
                    sequence: self
                        .params
                        .with_sequence
                        .then(|| String::from_utf8_lossy(&seq[s..=e]).into_owned()),
                });
            }
        }
        res
    }
}

/// Results of a whole run, one entry per pair (same order as the pair specs).
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct EpcrReport {
    pub sequences: u64,
    pub bp: u64,
    pub results: Vec<PairResult>,
}

impl EpcrReport {
    pub fn new(n_pairs: usize) -> Self {
        Self {
            sequences: 0,
            bp: 0,
            results: vec![PairResult::default(); n_pairs],
        }
    }

    /// Merge the hits of the next record (records must be added in input order).
    pub fn add(&mut self, len: usize, hits: RecordHits, max_amplicons: usize) {
        self.sequences += 1;
        self.bp += len as u64;
        for (i, h) in hits.pairs {
            let r = &mut self.results[i];
            r.amplicons += h.amplicons;
            r.perfect += h.perfect;
            r.expected_found |= h.expected_found;
            r.saturated |= h.saturated;
            let room = max_amplicons.saturating_sub(r.list.len());
            r.list.extend(h.list.into_iter().take(room));
        }
    }

    /// Number of pairs per status.
    pub fn status_counts(
        &self,
        pairs: &[PairSpec],
    ) -> std::collections::BTreeMap<&'static str, u64> {
        let mut m = std::collections::BTreeMap::new();
        for (r, p) in self.results.iter().zip(pairs) {
            *m.entry(r.status(p.expected.is_some()).name()).or_default() += 1;
        }
        m
    }
}

/// Header of the amplicon table.
pub fn amplicons_header(with_sequence: bool) -> String {
    let mut h = String::from(
        "pair_id\tseq_id\tstart\tend\tsize\tfwd_primer\tfwd_start\tfwd_end\tfwd_mismatches\trev_primer\trev_start\trev_end\trev_mismatches\texpected",
    );
    if with_sequence {
        h.push_str("\tsequence");
    }
    h.push('\n');
    h
}

/// Append every listed amplicon.
pub fn write_amplicons(
    out: &mut String,
    report: &EpcrReport,
    pairs: &[PairSpec],
    with_sequence: bool,
) {
    for (r, p) in report.results.iter().zip(pairs) {
        for a in &r.list {
            let _ = write!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                p.id,
                a.seq_id,
                a.start,
                a.end,
                a.size,
                a.fwd.primer,
                a.fwd.start,
                a.fwd.end,
                a.fwd.mismatches,
                a.rev.primer,
                a.rev.start,
                a.rev.end,
                a.rev.mismatches,
                if a.expected { "yes" } else { "no" }
            );
            if with_sequence {
                out.push('\t');
                out.push_str(a.sequence.as_deref().unwrap_or(""));
            }
            out.push('\n');
        }
    }
}

/// Per-pair summary table.
pub fn write_pairs(out: &mut String, report: &EpcrReport, pairs: &[PairSpec]) {
    out.push_str("pair_id\tforward\treverse\tstatus\tamplicons\tperfect_amplicons\texpected_found\tsaturated\n");
    for (r, p) in report.results.iter().zip(pairs) {
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            p.id,
            p.forward,
            p.reverse,
            r.status(p.expected.is_some()).name(),
            r.amplicons,
            r.perfect,
            match (&p.expected, r.expected_found) {
                (None, _) => "NA",
                (Some(_), true) => "yes",
                (Some(_), false) => "no",
            },
            if r.saturated { "yes" } else { "no" }
        );
    }
}

fn is_primer_like(s: &str) -> bool {
    s.len() >= 8 && s.bytes().all(|b| iupac_mask(b) != 0)
}

/// Parse a primer list. Accepted layouts:
///
/// * a table with a header line naming the primer columns
///   (`left_primer`/`right_primer`, `forward`/`reverse`, `fwd`/`rev`, …),
///   e.g. nextSSR `.tsv` or `.primers.tsv` files; the identifier is taken from
///   `pair_id`, `id`, `name` or `seq_id`+`ssr_nr`(+`pair_rank`), else the row number;
/// * whitespace-separated lines `ID FORWARD REVERSE [...]` (in-silico PCR/min-silico PCR STS
///   files) or `FORWARD REVERSE`.
///
/// Empty lines, `#` comments and rows without both primers are skipped.
pub fn parse_primer_list(text: &str) -> Result<Vec<PairSpec>, ParamError> {
    let lines: Vec<&str> = text
        .lines()
        .map(|l| l.trim_end_matches('\r'))
        .filter(|l| !l.trim().is_empty() && !l.trim_start().starts_with('#'))
        .collect();
    let Some(first) = lines.first() else {
        return Ok(Vec::new());
    };
    let norm = |s: &str| s.trim().to_ascii_lowercase().replace([' ', '-'], "_");
    let head: Vec<String> = first.split('\t').map(norm).collect();
    let find = |names: &[&str]| head.iter().position(|h| names.contains(&h.as_str()));
    let fcol = find(&[
        "left_primer",
        "forward_primer",
        "forward",
        "fwd",
        "fwd_primer",
        "primer_left",
        "primer_left_sequence",
        "f",
        "left",
    ]);
    let rcol = find(&[
        "right_primer",
        "reverse_primer",
        "reverse",
        "rev",
        "rev_primer",
        "primer_right",
        "primer_right_sequence",
        "r",
        "right",
    ]);
    let mut out = Vec::new();
    if let (Some(fc), Some(rc)) = (fcol, rcol) {
        let idc = find(&["pair_id", "id", "name", "marker", "locus", "sts"]);
        let seqc = find(&["seq_id"]);
        let nrc = find(&["ssr_nr"]);
        let rankc = find(&["pair_rank"]);
        for (row, l) in lines[1..].iter().enumerate() {
            let f: Vec<&str> = l.split('\t').map(str::trim).collect();
            let get = |c: Option<usize>| c.and_then(|c| f.get(c).copied()).unwrap_or("");
            let (fw, rv) = (get(Some(fc)), get(Some(rc)));
            if fw.is_empty() || rv.is_empty() {
                continue;
            }
            let id = if !get(idc).is_empty() {
                get(idc).to_string()
            } else if !get(seqc).is_empty() && !get(nrc).is_empty() {
                let mut id = format!("{}_{}", get(seqc), get(nrc));
                if !get(rankc).is_empty() {
                    id.push('_');
                    id.push_str(get(rankc));
                }
                id
            } else {
                format!("pair{}", row + 1)
            };
            out.push(PairSpec {
                id,
                forward: fw.to_string(),
                reverse: rv.to_string(),
                expected: None,
            });
        }
        return Ok(out);
    }
    for (row, l) in lines.iter().enumerate() {
        let f: Vec<&str> = l.split_whitespace().collect();
        let (id, fw, rv) = if f.len() >= 3 && is_primer_like(f[1]) && is_primer_like(f[2]) {
            (f[0].to_string(), f[1], f[2])
        } else if f.len() >= 2 && is_primer_like(f[0]) && is_primer_like(f[1]) {
            (format!("pair{}", row + 1), f[0], f[1])
        } else {
            return Err(ParamError(format!(
                "primer list line {}: expected `ID FORWARD REVERSE` or a header with forward/reverse columns: {l:?}",
                row + 1
            )));
        };
        out.push(PairSpec {
            id,
            forward: fw.to_string(),
            reverse: rv.to_string(),
            expected: None,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::motif::revcomp;

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

    fn pair(id: &str, f: &[u8], r: &[u8]) -> PairSpec {
        PairSpec {
            id: id.into(),
            forward: String::from_utf8(f.to_vec()).unwrap(),
            reverse: String::from_utf8(r.to_vec()).unwrap(),
            expected: None,
        }
    }

    /// Brute force reference: every site with the same mismatch rules.
    fn brute_sites(seq: &[u8], primer: &[u8], p: &EpcrParams) -> Vec<(bool, usize, u8)> {
        let mut out = Vec::new();
        let fw: Vec<u8> = primer.iter().map(|&b| iupac_mask(b)).collect();
        let rc: Vec<u8> = fw.iter().rev().map(|&m| comp_mask(m)).collect();
        let l = fw.len();
        let t3 = p.three_prime_len.min(l);
        for (minus, pat) in [(false, &fw), (true, &rc)] {
            for s in 0..=seq.len().saturating_sub(l) {
                let mut mm = 0;
                let mut mm3 = 0;
                for i in 0..l {
                    if TMASK[seq[s + i] as usize] & pat[i] == 0 {
                        mm += 1;
                        let in3 = if minus { i < t3 } else { i >= l - t3 };
                        if in3 {
                            mm3 += 1;
                        }
                    }
                }
                if mm <= p.max_mismatches && mm3 <= p.max_three_prime_mismatches {
                    out.push((minus, s, mm as u8));
                }
            }
        }
        out.sort();
        out
    }

    #[test]
    fn finds_expected_product_on_both_strands() {
        let mut g = rnd_seq(3000, 3);
        let f = g[1000..1020].to_vec();
        let r = revcomp(&g[1180..1200]);
        let spec = PairSpec {
            expected: Some(Expected {
                seq_id: "chr".into(),
                start: 1001,
                end: 1200,
            }),
            ..pair("p1", &f, &r)
        };
        let e = Epcr::new(vec![spec.clone()], EpcrParams::default()).unwrap();
        let hits = e.scan("chr", &g);
        let mut rep = EpcrReport::new(1);
        rep.add(g.len(), hits, 100);
        let r0 = &rep.results[0];
        assert_eq!(r0.amplicons, 1);
        assert!(r0.expected_found);
        assert_eq!(r0.status(true), Status::Specific);
        let a = &r0.list[0];
        assert_eq!((a.start, a.end, a.size), (1001, 1200, 200));
        assert_eq!((a.fwd.primer, a.rev.primer), ('F', 'R'));

        // the same locus on the other strand of a second sequence
        g = revcomp(&g);
        let hits = e.scan("other", &g);
        rep.add(g.len(), hits, 100);
        let r0 = &rep.results[0];
        assert_eq!(r0.amplicons, 2);
        assert_eq!(r0.status(true), Status::Nonspecific);
        let a = &r0.list[1];
        assert_eq!((a.seq_id.as_str(), a.size), ("other", 200));
        assert_eq!((a.fwd.primer, a.rev.primer), ('R', 'F'));
        assert!(!a.expected);
    }

    #[test]
    fn mismatch_rules_match_brute_force() {
        let g = {
            let mut g = rnd_seq(20_000, 11);
            // plant mutated copies of a primer
            let primer = g[500..522].to_vec();
            for (k, pos) in [3000usize, 7000, 9000, 15000].into_iter().enumerate() {
                let mut c = primer.clone();
                for i in 0..k {
                    let j = 3 + 5 * i;
                    c[j] = if c[j] == b'A' { b'C' } else { b'A' };
                }
                if k == 3 {
                    c = revcomp(&c);
                }
                g[pos..pos + c.len()].copy_from_slice(&c);
            }
            // a mismatch inside the 3' region
            let mut c = primer.clone();
            c[20] = if c[20] == b'A' { b'C' } else { b'A' };
            g[12000..12022].copy_from_slice(&c);
            g
        };
        let primer = g[500..522].to_vec();
        for (m, m3) in [(0, 0), (1, 0), (2, 0), (2, 1), (3, 1)] {
            let p = EpcrParams {
                max_mismatches: m,
                max_three_prime_mismatches: m3,
                ..EpcrParams::default()
            };
            let e = Epcr::new(vec![pair("x", &primer, b"ACGTACGTACGTACGTAC")], p.clone()).unwrap();
            let mut got: Vec<(bool, usize, u8)> = e
                .sites(&g)
                .into_iter()
                .filter(|s| e.patterns[s.pat as usize].role == 0)
                .map(|s| (e.patterns[s.pat as usize].minus, s.start, s.mm))
                .collect();
            got.sort();
            assert_eq!(got, brute_sites(&g, &primer, &p), "m={m} m3={m3}");
        }
    }

    #[test]
    fn chunked_scan_matches_single_pass() {
        // sites straddling chunk borders must be reported exactly once
        let mut g = rnd_seq(CHUNK * 2 + 5000, 5);
        let primer = rnd_seq(20, 77);
        for pos in [CHUNK - 10, CHUNK - 20, 2 * CHUNK - 1, 100] {
            g[pos..pos + 20].copy_from_slice(&primer);
        }
        let e = Epcr::new(
            vec![pair("x", &primer, &revcomp(&primer))],
            EpcrParams::default(),
        )
        .unwrap();
        let mut all = e.sites(&g);
        let mut single = Vec::new();
        e.scan_range(&g, 0, g.len(), &mut single);
        let key = |s: &Site| (s.pat, s.start, s.mm);
        all.sort_by_key(key);
        single.sort_by_key(key);
        assert_eq!(all.len(), single.len());
        assert!(all.iter().zip(&single).all(|(a, b)| key(a) == key(b)));
        assert!(all.iter().filter(|s| s.mm == 0).count() >= 6);
    }

    #[test]
    fn product_limits_and_counting() {
        let mut g = rnd_seq(5000, 9);
        let f = rnd_seq(20, 1);
        let r = rnd_seq(20, 2);
        let rc = revcomp(&r);
        // F ... R at 300 bp, and F ... R at 2500 bp (too long by default)
        g[100..120].copy_from_slice(&f);
        g[380..400].copy_from_slice(&rc);
        g[1000..1020].copy_from_slice(&f);
        g[3480..3500].copy_from_slice(&rc);
        let e = Epcr::new(vec![pair("x", &f, &r)], EpcrParams::default()).unwrap();
        let h = e.scan("s", &g);
        let mut rep = EpcrReport::new(1);
        rep.add(g.len(), h, 100);
        // 101..400 (300 bp) and 1001..3500 is too long; 101..3500 too
        assert_eq!(rep.results[0].amplicons, 1);
        assert_eq!(rep.results[0].list[0].size, 300);
        let e = Epcr::new(
            vec![pair("x", &f, &r)],
            EpcrParams {
                max_product: 5000,
                ..EpcrParams::default()
            },
        )
        .unwrap();
        let mut rep = EpcrReport::new(1);
        rep.add(g.len(), e.scan("s", &g), 2);
        // 101..400, 101..3500 and 1001..3500
        assert_eq!(rep.results[0].amplicons, 3);
        assert_eq!(rep.results[0].perfect, 3);
        assert_eq!(rep.results[0].list.len(), 2);
    }

    #[test]
    fn iupac_and_soft_masking() {
        let mut g = rnd_seq(2000, 13);
        let f = g[200..220].to_vec();
        let r = revcomp(&g[400..420]);
        g[200..220].make_ascii_lowercase();
        let mut fd = f.clone();
        fd[2] = b'N';
        fd[5] = match f[5] {
            b'A' | b'G' => b'R',
            _ => b'Y',
        };
        let e = Epcr::new(
            vec![pair("x", &fd, &r)],
            EpcrParams {
                max_mismatches: 0,
                ..EpcrParams::default()
            },
        )
        .unwrap();
        let mut rep = EpcrReport::new(1);
        rep.add(g.len(), e.scan("s", &g), 10);
        assert_eq!(rep.results[0].amplicons, 1);
        // template N never matches
        g[205] = b'N';
        let mut rep = EpcrReport::new(1);
        rep.add(g.len(), e.scan("s", &g), 10);
        assert_eq!(rep.results[0].amplicons, 0);
        assert!(Epcr::new(vec![pair("x", b"ACGTXACGTACGT", &r)], EpcrParams::default()).is_err());
    }

    #[test]
    fn primer_lists() {
        let t = "# comment\nm1 ACGTACGTACGTACGTAA TTGCATGCATGCATGCAT 150\nACGTACGTACGTACGTAA\tTTGCATGCATGCATGCAT\n";
        let v = parse_primer_list(t).unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].id, "m1");
        assert_eq!(v[1].id, "pair2");
        let t = "seq_id\tssr_nr\tpair_rank\tleft_primer\tright_primer\nchr1\t3\t1\tACGTACGTACGTACGTAA\tTTGCATGCATGCATGCAT\nchr1\t4\t1\t\t\n";
        let v = parse_primer_list(t).unwrap();
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].id, "chr1_3_1");
        assert!(parse_primer_list("hello world\n").is_err());
    }
}

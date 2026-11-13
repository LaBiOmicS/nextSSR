//! Run statistics and the MISA `.statistics` report.

use std::collections::BTreeMap;

use serde::Serialize;

use crate::locus::Locus;
use crate::motif::misa_class;
use crate::params::{Definition, Mode};

/// Aggregated statistics (mergeable across threads / batches).
#[derive(Debug, Clone, Default, Serialize, PartialEq)]
pub struct Stats {
    /// Number of sequences examined.
    pub sequences: u64,
    /// Total size of examined sequences (bp).
    pub total_bp: u64,
    /// Number of individual SSRs (compound members counted separately).
    pub ssrs: u64,
    /// Number of reported loci (single SSRs + compound SSRs).
    pub loci: u64,
    /// Number of compound SSRs.
    pub compounds: u64,
    /// Number of SSRs that are members of a compound SSR.
    pub ssrs_in_compounds: u64,
    /// MISA's counter "Number of SSRs present in compound formation"
    /// (counts members - 1 per compound, a quirk of misa.pl).
    pub misa_ssr_in_compound: u64,
    /// Number of imperfect SSRs.
    pub imperfect: u64,
    /// Sequences containing at least one SSR.
    pub seqs_with_ssr: u64,
    /// Sequences containing exactly one SSR.
    pub seqs_with_one_ssr: u64,
    /// unit size -> number of SSRs
    pub per_unit: BTreeMap<usize, u64>,
    /// motif -> repeats -> count
    pub per_motif: BTreeMap<String, BTreeMap<usize, u64>>,
    /// Largest number of repeats observed (MISA starts at 1).
    pub max_repeats: usize,
    /// Number of loci with a designed primer pair.
    pub primers_designed: u64,
    /// Primer design failure reason -> count
    pub primer_failures: BTreeMap<String, u64>,
    /// Base counts of the examined sequences: A, C, G, T, other.
    pub composition: [u64; 5],
    /// unit size -> total length of SSRs (bp)
    pub per_unit_bp: BTreeMap<usize, u64>,
    /// unit size -> SSR length (bp) -> count
    pub per_unit_len: BTreeMap<usize, BTreeMap<usize, u64>>,
    /// motif -> total length of SSRs (bp)
    pub per_motif_bp: BTreeMap<String, u64>,
    /// Total length of all reported loci (bp).
    pub loci_bp: u64,
    /// Compound SSRs with adjacent or overlapping members (`c*`).
    pub compounds_contiguous: u64,
    /// Sum of the purity of imperfect SSRs.
    pub imperfect_purity_sum: f64,
    /// Best-pair product size (bp) -> count.
    pub primer_products: BTreeMap<usize, u64>,
    /// Sum of the Tm of designed primers (both primers of the best pair).
    pub primer_tm_sum: f64,
    /// Per-sequence figures, in input order (the first [`MAX_PER_SEQUENCE`]).
    #[serde(skip)]
    pub per_sequence: Vec<SeqStat>,
}

/// Sequences kept in [`Stats::per_sequence`].
pub const MAX_PER_SEQUENCE: usize = 100_000;

/// Figures of one sequence.
#[derive(Debug, Clone, Default, Serialize, PartialEq)]
pub struct SeqStat {
    pub id: String,
    pub length: u64,
    /// G+C over A+C+G+T (%).
    pub gc: f64,
    pub ssrs: u64,
    pub loci: u64,
    pub ssr_bp: u64,
}

impl Stats {
    pub fn new() -> Self {
        Self {
            max_repeats: 1,
            ..Default::default()
        }
    }

    /// Account for one analysed sequence.
    pub fn add_record(&mut self, id: &str, seq: &[u8], loci: &[Locus]) {
        self.sequences += 1;
        self.total_bp += seq.len() as u64;
        let mut comp = [0u64; 5];
        for &b in seq {
            comp[crate::motif::CODE[b as usize] as usize] += 1;
        }
        for (a, c) in self.composition.iter_mut().zip(comp) {
            *a += c;
        }
        let mut n_ssr = 0u64;
        let mut ssr_bp = 0u64;
        for l in loci {
            self.loci += 1;
            self.loci_bp += l.size().max(0) as u64;
            if l.is_compound() {
                self.compounds += 1;
                if l.ssr_type == "c*" {
                    self.compounds_contiguous += 1;
                }
                self.ssrs_in_compounds += l.members.len() as u64;
                self.misa_ssr_in_compound += l.members.len() as u64 - 1;
            }
            if let Some(p) = &l.primers {
                match p {
                    crate::primer::PrimerOutcome::Designed { pairs } => {
                        self.primers_designed += 1;
                        let b = &pairs[0];
                        *self.primer_products.entry(b.product_size).or_default() += 1;
                        self.primer_tm_sum += b.left.tm + b.right.tm;
                    }
                    crate::primer::PrimerOutcome::Failed { reason } => {
                        *self.primer_failures.entry(reason.clone()).or_default() += 1
                    }
                }
            }
            for s in &l.members {
                n_ssr += 1;
                let len = s.len();
                ssr_bp += len as u64;
                if let Some(i) = &s.imperfect {
                    self.imperfect += 1;
                    self.imperfect_purity_sum += i.purity;
                }
                *self.per_unit.entry(s.unit).or_default() += 1;
                *self.per_unit_bp.entry(s.unit).or_default() += len as u64;
                *self
                    .per_unit_len
                    .entry(s.unit)
                    .or_default()
                    .entry(len)
                    .or_default() += 1;
                *self.per_motif_bp.entry(s.motif.clone()).or_default() += len as u64;
                *self
                    .per_motif
                    .entry(s.motif.clone())
                    .or_default()
                    .entry(s.repeats)
                    .or_default() += 1;
                self.max_repeats = self.max_repeats.max(s.repeats);
            }
        }
        self.ssrs += n_ssr;
        if self.per_sequence.len() < MAX_PER_SEQUENCE {
            let acgt = comp[..4].iter().sum::<u64>();
            self.per_sequence.push(SeqStat {
                id: id.to_string(),
                length: seq.len() as u64,
                gc: if acgt > 0 {
                    100.0 * (comp[1] + comp[2]) as f64 / acgt as f64
                } else {
                    0.0
                },
                ssrs: n_ssr,
                loci: loci.len() as u64,
                ssr_bp,
            });
        }
        if n_ssr > 0 {
            self.seqs_with_ssr += 1;
            if n_ssr == 1 {
                self.seqs_with_one_ssr += 1;
            }
        }
    }

    /// Merge statistics from another batch.
    pub fn merge(&mut self, o: &Stats) {
        self.sequences += o.sequences;
        self.total_bp += o.total_bp;
        self.ssrs += o.ssrs;
        self.loci += o.loci;
        self.compounds += o.compounds;
        self.ssrs_in_compounds += o.ssrs_in_compounds;
        self.misa_ssr_in_compound += o.misa_ssr_in_compound;
        self.imperfect += o.imperfect;
        self.seqs_with_ssr += o.seqs_with_ssr;
        self.seqs_with_one_ssr += o.seqs_with_one_ssr;
        for (k, v) in &o.per_unit {
            *self.per_unit.entry(*k).or_default() += v;
        }
        for (m, reps) in &o.per_motif {
            let e = self.per_motif.entry(m.clone()).or_default();
            for (r, c) in reps {
                *e.entry(*r).or_default() += c;
            }
        }
        self.max_repeats = self.max_repeats.max(o.max_repeats);
        self.primers_designed += o.primers_designed;
        for (k, v) in &o.primer_failures {
            *self.primer_failures.entry(k.clone()).or_default() += v;
        }
        for (a, c) in self.composition.iter_mut().zip(o.composition) {
            *a += c;
        }
        for (k, v) in &o.per_unit_bp {
            *self.per_unit_bp.entry(*k).or_default() += v;
        }
        for (u, lens) in &o.per_unit_len {
            let e = self.per_unit_len.entry(*u).or_default();
            for (l, c) in lens {
                *e.entry(*l).or_default() += c;
            }
        }
        for (k, v) in &o.per_motif_bp {
            *self.per_motif_bp.entry(k.clone()).or_default() += v;
        }
        self.loci_bp += o.loci_bp;
        self.compounds_contiguous += o.compounds_contiguous;
        self.imperfect_purity_sum += o.imperfect_purity_sum;
        for (k, v) in &o.primer_products {
            *self.primer_products.entry(*k).or_default() += v;
        }
        self.primer_tm_sum += o.primer_tm_sum;
        let room = MAX_PER_SEQUENCE.saturating_sub(self.per_sequence.len());
        self.per_sequence
            .extend(o.per_sequence.iter().take(room).cloned());
    }

    /// G+C content over A+C+G+T (%).
    pub fn gc_percent(&self) -> f64 {
        let c = &self.composition;
        let acgt = c[0] + c[1] + c[2] + c[3];
        if acgt == 0 {
            0.0
        } else {
            100.0 * (c[1] + c[2]) as f64 / acgt as f64
        }
    }

    /// SSR counts grouped by MISA class (`AG/CT`), repeats -> count.
    pub fn per_class(&self) -> BTreeMap<String, BTreeMap<usize, u64>> {
        let mut out: BTreeMap<String, BTreeMap<usize, u64>> = BTreeMap::new();
        for (m, reps) in &self.per_motif {
            let e = out.entry(misa_class(m)).or_default();
            for (r, c) in reps {
                *e.entry(*r).or_default() += c;
            }
        }
        out
    }
}

fn by_len_then_lex(a: &String, b: &String) -> std::cmp::Ordering {
    a.len().cmp(&b.len()).then_with(|| a.cmp(b))
}

/// Render the `.statistics` report in the exact layout of `misa.pl` v2.0.
///
/// In [`Mode::Misa`] the "SSRs present in compound formation" counter
/// reproduces misa.pl (members - 1 per compound); in [`Mode::Native`] it is
/// the true number of SSRs belonging to compounds.
pub fn misa_statistics(
    st: &Stats,
    def: &Definition,
    max_interruption: Option<i64>,
    source: &str,
    mode: Mode,
) -> String {
    let mut o = String::new();
    o.push_str("Specifications\n==============\n\n");
    o.push_str(&format!("Sequence source file: \"{source}\"\n\n"));
    o.push_str("Definement of microsatellites (unit size / minimum number of repeats):\n");
    for (u, r) in &def.min_repeats {
        o.push_str(&format!("({u}/{r}) "));
    }
    o.push('\n');
    if let Some(d) = max_interruption {
        if d > 0 {
            o.push_str(&format!(
                "\nMaximal number of bases interrupting 2 SSRs in a compound microsatellite:  {d}\n"
            ));
        }
    }
    o.push_str("\n\n\n");

    let total = if st.per_unit.is_empty() {
        String::new()
    } else {
        st.per_unit.values().sum::<u64>().to_string()
    };
    let in_compound = match mode {
        Mode::Misa => st.misa_ssr_in_compound,
        Mode::Native => st.ssrs_in_compounds,
    };
    o.push_str("RESULTS OF MICROSATELLITE SEARCH\n================================\n\n");
    o.push_str(&format!(
        "Total number of sequences examined:              {}\n",
        st.sequences
    ));
    o.push_str(&format!(
        "Total size of examined sequences (bp):           {}\n",
        st.total_bp
    ));
    o.push_str(&format!(
        "Total number of identified SSRs:                 {total}\n"
    ));
    o.push_str(&format!(
        "Number of SSR containing sequences:              {}\n",
        st.seqs_with_ssr
    ));
    o.push_str(&format!(
        "Number of sequences containing more than 1 SSR:  {}\n",
        st.seqs_with_ssr - st.seqs_with_one_ssr
    ));
    o.push_str(&format!(
        "Number of SSRs present in compound formation:    {in_compound}\n\n\n"
    ));

    o.push_str("Distribution to different repeat type classes\n---------------------------------------------\n\n");
    o.push_str("Unit size\tNumber of SSRs\n");
    for (u, c) in &st.per_unit {
        o.push_str(&format!("{u}\t{c}\n"));
    }
    o.push('\n');

    // misa.pl only sets $min_repeats while scanning sequences.
    let min_r = if st.sequences > 0 {
        def.min_threshold()
    } else {
        1000
    };
    let max_r = st.max_repeats;
    let header = |o: &mut String| {
        o.push_str("Repeats");
        for i in min_r..=max_r {
            o.push_str(&format!("\t{i}"));
        }
        o.push_str("\ttotal\n");
    };
    let row = |o: &mut String, name: &str, unit: usize, counts: &BTreeMap<usize, u64>| {
        o.push_str(name);
        let thr = def.get(unit).unwrap_or(0);
        let mut total = 0u64;
        for j in min_r..=max_r {
            if j < thr {
                o.push_str("\t-");
                continue;
            }
            match counts.get(&j) {
                Some(c) if *c > 0 => {
                    o.push_str(&format!("\t{c}"));
                    total += c;
                }
                _ => o.push('\t'),
            }
        }
        o.push_str(&format!("\t{total}\n"));
    };

    o.push_str("Frequency of identified SSR motifs\n----------------------------------\n\n");
    header(&mut o);
    let mut motifs: Vec<&String> = st.per_motif.keys().collect();
    motifs.sort_by(|a, b| by_len_then_lex(a, b));
    for m in motifs {
        row(&mut o, m, m.len(), &st.per_motif[m]);
    }
    o.push('\n');

    o.push_str("Frequency of classified repeat types (considering sequence complementary)\n");
    o.push_str("-------------------------------------------------------------------------\n\n");
    header(&mut o);
    let classes = st.per_class();
    let mut names: Vec<&String> = classes.keys().collect();
    names.sort_by(|a, b| by_len_then_lex(a, b));
    for c in names {
        let unit = (c.len() - 1) / 2;
        row(&mut o, c, unit, &classes[c]);
    }
    o
}

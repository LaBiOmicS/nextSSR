//! Primer table (`.primers.tsv`): one row per designed pair, optionally with
//! the in-silico PCR specificity of each pair.

use std::fmt::Write;

use super::fmt_f;
use crate::epcr::{EpcrReport, Expected, PairSpec, Status};
use crate::pipeline::RecordResult;
use crate::primer::PrimerPair;

/// A designed pair together with its locus.
#[derive(Debug, Clone, PartialEq)]
pub struct PrimerRow {
    pub seq_id: String,
    pub nr: usize,
    pub ssr_type: String,
    pub pattern: String,
    pub ssr_start: usize,
    pub ssr_end: usize,
    /// 1 = best pair of the locus.
    pub rank: usize,
    pub pair: PrimerPair,
}

impl PrimerRow {
    /// `<seq_id>_<ssr_nr>_<rank>`.
    pub fn id(&self) -> String {
        format!("{}_{}_{}", self.seq_id, self.nr, self.rank)
    }

    /// in-silico PCR specification whose expected product is the designed one.
    pub fn spec(&self) -> PairSpec {
        PairSpec {
            id: self.id(),
            forward: self.pair.left.seq.clone(),
            reverse: self.pair.right.seq.clone(),
            expected: Some(Expected {
                seq_id: self.seq_id.clone(),
                start: self.pair.left.start,
                end: self.pair.right.end,
            }),
        }
    }
}

/// Rows of every designed pair of a record.
pub fn rows(r: &RecordResult) -> Vec<PrimerRow> {
    let mut out = Vec::new();
    for l in &r.loci {
        let Some(o) = &l.primers else { continue };
        for (i, p) in o.pairs().iter().enumerate() {
            out.push(PrimerRow {
                seq_id: r.id.clone(),
                nr: l.nr,
                ssr_type: l.ssr_type.clone(),
                pattern: l.pattern.clone(),
                ssr_start: l.start,
                ssr_end: l.end,
                rank: i + 1,
                pair: p.clone(),
            });
        }
    }
    out
}

pub fn header(with_epcr: bool) -> String {
    let mut h = String::from(
        "pair_id\tseq_id\tssr_nr\tssr_type\tpattern\tssr_start\tssr_end\tpair_rank\tleft_primer\tleft_start\tleft_end\tleft_tm\tleft_gc\tright_primer\tright_start\tright_end\tright_tm\tright_gc\tproduct_size\tpair_penalty",
    );
    if with_epcr {
        h.push_str("\tepcr_status\tepcr_amplicons\tepcr_perfect_amplicons\trecommended");
    }
    h.push('\n');
    h
}

/// Write `rows`; `epcr` must hold one result per row, in the same order.
pub fn write(out: &mut String, rows: &[PrimerRow], epcr: Option<&EpcrReport>) {
    // recommended pair per locus: the best-ranked specific pair
    let recommended: Vec<bool> = match epcr {
        Some(rep) => {
            let mut v = vec![false; rows.len()];
            let mut i = 0;
            while i < rows.len() {
                let mut j = i;
                while j < rows.len() && rows[j].seq_id == rows[i].seq_id && rows[j].nr == rows[i].nr
                {
                    j += 1;
                }
                if let Some(k) = (i..j).find(|&k| rep.results[k].status(true) == Status::Specific) {
                    v[k] = true;
                }
                i = j;
            }
            v
        }
        None => Vec::new(),
    };
    for (k, r) in rows.iter().enumerate() {
        let p = &r.pair;
        let _ = write!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.id(),
            r.seq_id,
            r.nr,
            r.ssr_type,
            r.pattern,
            r.ssr_start,
            r.ssr_end,
            r.rank,
            p.left.seq,
            p.left.start,
            p.left.end,
            fmt_f(p.left.tm),
            fmt_f(p.left.gc),
            p.right.seq,
            p.right.start,
            p.right.end,
            fmt_f(p.right.tm),
            fmt_f(p.right.gc),
            p.product_size,
            fmt_f(p.penalty)
        );
        if let Some(rep) = epcr {
            let x = &rep.results[k];
            let _ = write!(
                out,
                "\t{}\t{}\t{}\t{}",
                x.status(true).name(),
                x.amplicons,
                x.perfect,
                if recommended[k] { "yes" } else { "no" }
            );
        }
        out.push('\n');
    }
}

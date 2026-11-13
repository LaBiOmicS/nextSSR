//! Rich tab-separated output.

use std::fmt::Write;

use super::fmt_f;
use crate::locus::Locus;
use crate::motif::canonical;
use crate::pipeline::RecordResult;
use crate::primer::PrimerOutcome;

/// Which optional column groups are written.
#[derive(Debug, Clone, Copy, Default)]
pub struct TsvColumns {
    pub primers: bool,
    pub flanks: bool,
}

/// Header line.
pub fn header(c: TsvColumns) -> String {
    let mut h = String::from(
        "seq_id\tssr_nr\tssr_type\tpattern\tmotif\tcanonical_motif\tunit\trepeats\tsize\tstart\tend\tpurity\tinterruptions\tn_members",
    );
    if c.primers {
        h.push_str("\tprimer_status\tleft_primer\tleft_start\tleft_end\tleft_tm\tleft_gc\tright_primer\tright_start\tright_end\tright_tm\tright_gc\tproduct_size\tpair_penalty");
    }
    if c.flanks {
        h.push_str("\tflank_5p\tflank_3p");
    }
    h.push('\n');
    h
}

fn join<T: ToString>(it: impl Iterator<Item = T>) -> String {
    it.map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}

fn purity(l: &Locus) -> String {
    if l.is_compound() {
        return String::new();
    }
    match &l.members[0].imperfect {
        Some(i) => fmt_f(i.purity),
        None => "1".into(),
    }
}

fn interruptions(l: &Locus) -> usize {
    if l.is_compound() {
        l.members.len() - 1
    } else {
        l.members[0]
            .imperfect
            .as_ref()
            .map_or(0, |i| i.interruptions)
    }
}

/// Append the rows of one record.
pub fn write_record(out: &mut String, r: &RecordResult, c: TsvColumns) {
    for l in &r.loci {
        let m = &l.members;
        let _ = write!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.id,
            l.nr,
            l.ssr_type,
            l.pattern,
            join(m.iter().map(|s| s.motif.clone())),
            join(m.iter().map(|s| canonical(&s.motif))),
            join(m.iter().map(|s| s.unit)),
            join(m.iter().map(|s| s.repeats)),
            l.size(),
            l.start,
            l.end,
            purity(l),
            interruptions(l),
            m.len()
        );
        if c.primers {
            match &l.primers {
                Some(PrimerOutcome::Designed(p)) => {
                    let _ = write!(
                        out,
                        "\tok\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
                }
                Some(PrimerOutcome::Failed { reason }) => {
                    let _ = write!(out, "\t{reason}{}", "\t".repeat(12));
                }
                None => {
                    let _ = write!(out, "\tnot_requested{}", "\t".repeat(12));
                }
            }
        }
        if c.flanks {
            let _ = write!(
                out,
                "\t{}\t{}",
                l.flank5.as_deref().unwrap_or(""),
                l.flank3.as_deref().unwrap_or("")
            );
        }
        out.push('\n');
    }
}

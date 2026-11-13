//! GFF3 output (specification 1.26) with Sequence Ontology feature types.
//!
//! * single SSRs: `monomeric_repeat` (unit 1),
//!   `dinucleotide_repeat_microsatellite_feature`,
//!   `trinucleotide_repeat_microsatellite_feature`,
//!   `tetranucleotide_repeat_microsatellite_feature` and `microsatellite` (unit >= 5),
//!   the same SO terms used by MISA-web;
//! * compound SSRs: a `repeat_region` parent with one child feature per SSR;
//! * primers: a `PCR_product` with two `primer_binding_site` children.
//!
//! Reserved characters in column 1 and attribute values are percent-encoded.

use std::fmt::Write;

use super::fmt_f;
use crate::locus::{Locus, Ssr};
use crate::motif::canonical;
use crate::pipeline::RecordResult;
use crate::primer::PrimerOutcome;

/// File header.
pub fn header() -> String {
    format!("##gff-version 3\n#!processor nextSSR {}\n", crate::VERSION)
}

/// SO feature type for a single SSR.
pub fn so_type(unit: usize) -> &'static str {
    match unit {
        1 => "monomeric_repeat",
        2 => "dinucleotide_repeat_microsatellite_feature",
        3 => "trinucleotide_repeat_microsatellite_feature",
        4 => "tetranucleotide_repeat_microsatellite_feature",
        _ => "microsatellite",
    }
}

/// Escape a seqid (column 1): `[a-zA-Z0-9.:^*$@!+_?-|]`, everything else encoded.
pub fn escape_seqid(s: &str) -> String {
    let mut o = String::with_capacity(s.len());
    for &b in s.as_bytes() {
        if b.is_ascii_alphanumeric() || b".:^*$@!+_?-|".contains(&b) {
            o.push(b as char);
        } else {
            let _ = write!(o, "%{b:02X}");
        }
    }
    o
}

/// Escape an attribute value (`;`, `=`, `&`, `,`, `%`, tabs and control characters).
pub fn escape_value(s: &str) -> String {
    let mut o = String::with_capacity(s.len());
    for c in s.chars() {
        if matches!(c, ';' | '=' | '&' | ',' | '%') || c.is_control() {
            let mut buf = [0u8; 4];
            for b in c.encode_utf8(&mut buf).bytes() {
                let _ = write!(o, "%{b:02X}");
            }
        } else {
            o.push(c);
        }
    }
    o
}

fn ssr_attrs(s: &Ssr) -> String {
    let mut a = format!(
        "motif={};canonical_motif={};unit={};repeats={}",
        escape_value(&s.motif),
        escape_value(&canonical(&s.motif)),
        s.unit,
        s.repeats
    );
    if let Some(i) = &s.imperfect {
        let _ = write!(
            a,
            ";imperfect=true;purity={};interruptions={};pattern={}",
            fmt_f(i.purity),
            i.interruptions,
            escape_value(&i.pattern)
        );
    }
    a
}

fn write_primers(out: &mut String, seqid: &str, parent: &str, l: &Locus) {
    let Some(PrimerOutcome::Designed(p)) = &l.primers else {
        return;
    };
    let pcr = format!("{parent}.pcr");
    let _ = writeln!(
        out,
        "{seqid}\tnextSSR\tPCR_product\t{}\t{}\t.\t.\t.\tID={pcr};Parent={parent};product_size={};pair_penalty={}",
        p.left.start,
        p.right.end,
        p.product_size,
        fmt_f(p.penalty)
    );
    for (pr, strand, tag) in [(&p.left, '+', "left"), (&p.right, '-', "right")] {
        let _ = writeln!(
            out,
            "{seqid}\tnextSSR\tprimer_binding_site\t{}\t{}\t.\t{strand}\t.\tID={pcr}.{tag};Parent={pcr};sequence={};tm={};gc={}",
            pr.start,
            pr.end,
            pr.seq,
            fmt_f(pr.tm),
            fmt_f(pr.gc)
        );
    }
}

/// Append the features of one record (nothing for records without SSRs).
pub fn write_record(out: &mut String, r: &RecordResult) {
    if r.loci.is_empty() {
        return;
    }
    let seqid = escape_seqid(&r.id);
    let _ = writeln!(out, "##sequence-region {seqid} 1 {}", r.length);
    let base = escape_value(&r.id);
    for l in &r.loci {
        let id = format!("{base}.ssr{}", l.nr);
        if l.is_compound() {
            let _ = writeln!(
                out,
                "{seqid}\tnextSSR\trepeat_region\t{}\t{}\t.\t.\t.\tID={id};Name={};ssr_type={};compound=true;n_members={}",
                l.start,
                l.end,
                escape_value(&l.pattern),
                escape_value(&l.ssr_type),
                l.members.len()
            );
            for (k, s) in l.members.iter().enumerate() {
                let _ = writeln!(
                    out,
                    "{seqid}\tnextSSR\t{}\t{}\t{}\t.\t.\t.\tID={id}.{};Parent={id};{}",
                    so_type(s.unit),
                    s.start,
                    s.end,
                    k + 1,
                    ssr_attrs(s)
                );
            }
        } else {
            let s = &l.members[0];
            let _ = writeln!(
                out,
                "{seqid}\tnextSSR\t{}\t{}\t{}\t.\t.\t.\tID={id};Name={};ssr_type={};{}",
                so_type(s.unit),
                l.start,
                l.end,
                escape_value(&l.pattern),
                escape_value(&l.ssr_type),
                ssr_attrs(s)
            );
        }
        write_primers(out, &seqid, &id, l);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escaping() {
        assert_eq!(escape_seqid("chr 1;x"), "chr%201%3Bx");
        assert_eq!(escape_value("a=b;c,d&e%"), "a%3Db%3Bc%2Cd%26e%25");
        assert_eq!(escape_value("(A)13{-1}(AG)8"), "(A)13{-1}(AG)8");
        assert_eq!(escape_value("x y\tz"), "x y%09z");
    }
}

//! Primer3 (`primer3_core` ≥ 2.0) Boulder-IO input, the successor of MISA's
//! `p3_in.pl`. Each locus becomes one record with a template window of
//! `product_max` bp on each side, the SSR as `SEQUENCE_TARGET`, the other SSR
//! loci of the window as `SEQUENCE_EXCLUDED_REGION` (unless primers may overlap
//! them) and the same size, Tm, GC, poly-X, clamp and reaction conditions as the
//! built-in designer. Positions are 0-based (`PRIMER_FIRST_BASE_INDEX=0`).

use std::fmt::Write;

use crate::pipeline::RecordResult;
use crate::primer::PrimerParams;

pub fn write_record(out: &mut String, r: &RecordResult, seq: &[u8], p: &PrimerParams) {
    for l in &r.loci {
        let a = l.start - 1;
        let lo = a.saturating_sub(p.product_max);
        let hi = (l.end + p.product_max).min(seq.len());
        let template = String::from_utf8_lossy(&seq[lo..hi]).to_ascii_uppercase();
        let _ = write!(
            out,
            "SEQUENCE_ID={}_{}\nSEQUENCE_TEMPLATE={}\nSEQUENCE_TARGET={},{}\n",
            r.id,
            l.nr,
            template,
            a - lo,
            l.end - a,
        );
        if p.avoid_other_ssrs {
            let excluded: Vec<String> = r
                .loci
                .iter()
                .filter(|o| o.nr != l.nr && o.end > lo && o.start - 1 < hi)
                .map(|o| {
                    let s = (o.start - 1).max(lo);
                    let e = o.end.min(hi);
                    format!("{},{}", s - lo, e - s)
                })
                .collect();
            if !excluded.is_empty() {
                let _ = writeln!(out, "SEQUENCE_EXCLUDED_REGION={}", excluded.join(" "));
            }
        }
        let t = &p.thermo;
        let _ = write!(
            out,
            "PRIMER_FIRST_BASE_INDEX=0\nPRIMER_NUM_RETURN={}\nPRIMER_PRODUCT_SIZE_RANGE={}-{}\n\
             PRIMER_MIN_SIZE={}\nPRIMER_OPT_SIZE={}\nPRIMER_MAX_SIZE={}\n\
             PRIMER_MIN_TM={}\nPRIMER_OPT_TM={}\nPRIMER_MAX_TM={}\nPRIMER_PAIR_MAX_DIFF_TM={}\n\
             PRIMER_MIN_GC={}\nPRIMER_MAX_GC={}\nPRIMER_MAX_POLY_X={}\nPRIMER_GC_CLAMP={}\n\
             PRIMER_SALT_MONOVALENT={}\nPRIMER_SALT_DIVALENT={}\nPRIMER_DNTP_CONC={}\nPRIMER_DNA_CONC={}\n=\n",
            p.num_return,
            p.product_min,
            p.product_max,
            p.min_size,
            p.opt_size,
            p.max_size,
            p.min_tm,
            p.opt_tm,
            p.max_tm,
            p.max_tm_diff,
            p.min_gc,
            p.max_gc,
            p.max_poly_x,
            p.gc_clamp,
            t.mv_conc,
            t.dv_conc,
            t.dntp_conc,
            t.dna_conc
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::parse_all;
    use crate::params::Config;
    use crate::pipeline::analyze_record;

    #[test]
    fn excluded_regions_are_relative_to_the_template() {
        let mut fa = b">s\n".to_vec();
        fa.extend(b"ACGTTGCA".repeat(40));
        fa.extend_from_slice(b"AGAGAGAGAGAGAGAG");
        fa.extend(b"CCGATTGA".repeat(10));
        fa.extend_from_slice(b"AAAAAAAAAAAA");
        fa.extend(b"GTCAGTTC".repeat(40));
        let rec = &parse_all(&fa)[0];
        let cfg = Config {
            max_interruption: None,
            ..Config::default()
        };
        let r = analyze_record(rec, &cfg);
        assert_eq!(r.loci.len(), 2);
        let mut out = String::new();
        write_record(&mut out, &r, &rec.seq, &PrimerParams::default());
        // locus 1 = (AG)8 at 321..336, locus 2 = (A)13 at 416..428
        // template of locus 1 starts at 20: target 300,16; excluded 395,13
        assert!(
            out.contains("SEQUENCE_TARGET=300,16\nSEQUENCE_EXCLUDED_REGION=395,13\n"),
            "{out}"
        );
        assert!(out.contains("SEQUENCE_ID=s_2\n"));
        assert_eq!(out.matches("\n=\n").count(), 2);
    }
}

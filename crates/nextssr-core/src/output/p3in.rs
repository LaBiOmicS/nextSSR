//! Primer3 (`primer3_core` ≥ 2.0) Boulder-IO input, the successor of MISA's
//! `p3_in.pl`. Each locus becomes one record with a template window of
//! `product_max` bp on each side and the SSR as `SEQUENCE_TARGET`.

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
            "SEQUENCE_ID={}_{}\nSEQUENCE_TEMPLATE={}\nSEQUENCE_TARGET={},{}\n\
             PRIMER_PRODUCT_SIZE_RANGE={}-{}\nPRIMER_MIN_SIZE={}\nPRIMER_OPT_SIZE={}\nPRIMER_MAX_SIZE={}\n\
             PRIMER_MIN_TM={}\nPRIMER_OPT_TM={}\nPRIMER_MAX_TM={}\nPRIMER_MIN_GC={}\nPRIMER_MAX_GC={}\n=\n",
            r.id,
            l.nr,
            template,
            a - lo,
            l.end - a,
            p.product_min,
            p.product_max,
            p.min_size,
            p.opt_size,
            p.max_size,
            p.min_tm,
            p.opt_tm,
            p.max_tm,
            p.min_gc,
            p.max_gc
        );
    }
}

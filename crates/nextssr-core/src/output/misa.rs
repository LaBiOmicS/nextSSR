//! `misa.pl`-compatible `.misa` table.

use std::fmt::Write;

use crate::pipeline::RecordResult;

/// Header line written by misa.pl.
pub const HEADER: &str = "ID\tSSR nr.\tSSR type\tSSR\tsize\tstart\tend\n";

/// Append the lines of one record.
pub fn write_record(out: &mut String, r: &RecordResult) {
    for l in &r.loci {
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.misa_id,
            l.nr,
            l.ssr_type,
            l.pattern,
            l.size(),
            l.start,
            l.end
        );
    }
}

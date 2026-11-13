//! BED6 output (0-based, half-open).

use std::fmt::Write;

use crate::pipeline::RecordResult;

/// Append the intervals of one record. The name column is
/// `<seq_id>.ssr<nr>|<type>|<pattern>` (truncated to 255 characters), the
/// score column holds the number of repeats of the first member (capped at 1000).
pub fn write_record(out: &mut String, r: &RecordResult) {
    for l in &r.loci {
        let mut name = format!("{}.ssr{}|{}|{}", r.id, l.nr, l.ssr_type, l.pattern);
        if name.len() > 255 {
            let mut cut = 255;
            while !name.is_char_boundary(cut) {
                cut -= 1;
            }
            name.truncate(cut);
        }
        let name: String = name
            .chars()
            .map(|c| if c.is_whitespace() { '_' } else { c })
            .collect();
        let score = l.members[0].repeats.min(1000);
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t.",
            r.id,
            l.start - 1,
            l.end,
            name,
            score
        );
    }
}

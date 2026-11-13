//! Incremental (push-based) FASTA parser.
//!
//! Bytes can be fed in arbitrary chunks, which allows the same parser to be used
//! by the native CLI (reading from files/stdin/gzip) and by the WebAssembly
//! build (reading browser `File` streams).
//!
//! Sequence cleaning follows `misa.pl`: digits, whitespace and `>` characters
//! are removed; every other character (including `N`, IUPAC codes and lowercase
//! letters) is kept verbatim.

/// One FASTA record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// 0-based index of the record in the input.
    pub index: usize,
    /// Full header line without the leading `>` and without the line terminator.
    pub header: String,
    /// Cleaned sequence bytes.
    pub seq: Vec<u8>,
}

impl Record {
    /// Identifier: first whitespace-delimited token of the header.
    pub fn id(&self) -> &str {
        self.header.split_ascii_whitespace().next().unwrap_or("")
    }

    /// Identifier as written by `misa.pl` in `.misa` files: header trimmed and
    /// every whitespace character replaced by `_`.
    pub fn misa_id(&self) -> String {
        misa_id(&self.header)
    }
}

/// `$id =~ s/^\s*//g; $id =~ s/\s*$//g; $id =~ s/\s/_/g;`
pub fn misa_id(header: &str) -> String {
    let is_ws = |c: char| matches!(c, ' ' | '\t' | '\n' | '\r' | '\x0b' | '\x0c');
    header
        .trim_matches(is_ws)
        .chars()
        .map(|c| if is_ws(c) { '_' } else { c })
        .collect()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    /// Before the first header: content is ignored (as in misa.pl).
    Preamble,
    Header,
    Sequence,
}

/// Push-based FASTA parser.
#[derive(Debug)]
pub struct FastaParser {
    state: State,
    at_line_start: bool,
    header: Vec<u8>,
    seq: Vec<u8>,
    count: usize,
}

impl Default for FastaParser {
    fn default() -> Self {
        Self::new()
    }
}

impl FastaParser {
    pub fn new() -> Self {
        Self {
            state: State::Preamble,
            at_line_start: true,
            header: Vec::new(),
            seq: Vec::new(),
            count: 0,
        }
    }

    /// Number of records emitted so far.
    pub fn records_emitted(&self) -> usize {
        self.count
    }

    /// Bytes buffered for the record currently being read.
    pub fn buffered_bytes(&self) -> usize {
        self.seq.len() + self.header.len()
    }

    fn emit(&mut self, out: &mut Vec<Record>) {
        let mut header = std::mem::take(&mut self.header);
        while header.last() == Some(&b'\r') {
            header.pop();
        }
        let header = String::from_utf8_lossy(&header).into_owned();
        let seq = std::mem::take(&mut self.seq);
        out.push(Record {
            index: self.count,
            header,
            seq,
        });
        self.count += 1;
    }

    /// Feed a chunk of bytes; complete records are appended to `out`.
    pub fn push(&mut self, chunk: &[u8], out: &mut Vec<Record>) {
        for &b in chunk {
            match self.state {
                State::Header => {
                    if b == b'\n' {
                        self.state = State::Sequence;
                        self.at_line_start = true;
                    } else {
                        self.header.push(b);
                    }
                }
                State::Preamble | State::Sequence => {
                    if b == b'>' && self.at_line_start {
                        if self.state == State::Sequence {
                            self.emit(out);
                        }
                        self.state = State::Header;
                        continue;
                    }
                    self.at_line_start = b == b'\n';
                    if self.state == State::Sequence && keep_seq_byte(b) {
                        self.seq.push(b);
                    }
                }
            }
        }
    }

    /// Signal end of input; the last record (if any) is appended to `out`.
    ///
    /// A final header line without a terminating newline is still reported as a
    /// record with an empty sequence.
    pub fn finish(&mut self, out: &mut Vec<Record>) {
        match self.state {
            State::Sequence | State::Header => self.emit(out),
            State::Preamble => {}
        }
        self.state = State::Preamble;
        self.at_line_start = true;
    }
}

#[inline]
fn keep_seq_byte(b: u8) -> bool {
    !(b.is_ascii_digit() || matches!(b, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c | b'>'))
}

/// Convenience: parse a complete in-memory FASTA text.
pub fn parse_all(data: &[u8]) -> Vec<Record> {
    let mut p = FastaParser::new();
    let mut out = Vec::new();
    p.push(data, &mut out);
    p.finish(&mut out);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_records_across_chunks() {
        let data = b"junk\n>seq1 desc here\r\nACGT 12\nacgtN\n>seq2\n\n>empty\n>s3\nAA";
        let whole = parse_all(data);
        for split in 1..data.len() {
            let mut p = FastaParser::new();
            let mut out = Vec::new();
            p.push(&data[..split], &mut out);
            p.push(&data[split..], &mut out);
            p.finish(&mut out);
            assert_eq!(out, whole, "split at {split}");
        }
        assert_eq!(whole.len(), 4);
        assert_eq!(whole[0].header, "seq1 desc here");
        assert_eq!(whole[0].id(), "seq1");
        assert_eq!(whole[0].misa_id(), "seq1_desc_here");
        assert_eq!(whole[0].seq, b"ACGTacgtN");
        assert_eq!(whole[1].seq, b"");
        assert_eq!(whole[2].id(), "empty");
        assert_eq!(whole[3].seq, b"AA");
    }

    #[test]
    fn empty_header_is_allowed() {
        let r = parse_all(b">\nACGT\n");
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].id(), "");
        assert_eq!(r[0].seq, b"ACGT");
    }
}

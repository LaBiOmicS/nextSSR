//! Output formats.
//!
//! Every writer appends text to a `String`, so the same code serves the CLI
//! (which streams the strings to files) and the WebAssembly build.

pub mod bed;
pub mod gff3;
pub mod misa;
pub mod p3in;
pub mod primers;
pub mod report;
pub mod tsv;

/// Output formats supported by nextSSR.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Format {
    /// Tab-separated table with all attributes (`.tsv`).
    Tsv,
    /// GFF3 with Sequence Ontology feature types (`.gff3`).
    Gff3,
    /// BED6 (`.bed`).
    Bed,
    /// misa.pl-compatible table (`.misa`).
    Misa,
    /// misa.pl-compatible statistics (`.statistics`).
    Statistics,
    /// Machine-readable summary (`.summary.json`).
    Json,
    /// Primer3 Boulder-IO input (`.p3in`).
    P3in,
    /// One row per designed primer pair (`.primers.tsv`).
    Primers,
    /// Self-contained HTML report (`.report.html`).
    Report,
}

impl Format {
    pub const ALL: [Format; 9] = [
        Format::Tsv,
        Format::Gff3,
        Format::Bed,
        Format::Misa,
        Format::Statistics,
        Format::Json,
        Format::P3in,
        Format::Primers,
        Format::Report,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Format::Tsv => "tsv",
            Format::Gff3 => "gff3",
            Format::Bed => "bed",
            Format::Misa => "misa",
            Format::Statistics => "statistics",
            Format::Json => "json",
            Format::P3in => "p3in",
            Format::Primers => "primers",
            Format::Report => "report",
        }
    }

    pub fn extension(self) -> &'static str {
        match self {
            Format::Json => "summary.json",
            Format::Primers => "primers.tsv",
            Format::Report => "report.html",
            other => other.name(),
        }
    }

    pub fn parse(s: &str) -> Option<Format> {
        Format::ALL.into_iter().find(|f| {
            f.name().eq_ignore_ascii_case(s.trim())
                || (s.trim() == "gff" && *f == Format::Gff3)
                || (s.trim() == "stats" && *f == Format::Statistics)
                || (s.trim() == "html" && *f == Format::Report)
        })
    }
}

/// Format a float without trailing noise (up to 3 decimals).
pub(crate) fn fmt_f(x: f64) -> String {
    let s = format!("{x:.3}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s.is_empty() || s == "-0" {
        "0".into()
    } else {
        s.into()
    }
}

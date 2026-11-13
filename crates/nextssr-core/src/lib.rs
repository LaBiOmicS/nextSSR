//! # nextssr-core
//!
//! Core engine of **nextSSR**, a fast microsatellite (SSR) finder designed as a
//! drop-in, improved replacement for MISA (Thiel et al. 2003; Beier et al. 2017).
//!
//! Two detection semantics are provided:
//!
//! * [`Mode::Misa`] reproduces `misa.pl` v2.0 exactly, including the masking
//!   behaviour of its regular expression scan and the quirks of its compound
//!   notation, so `.misa` / `.statistics` outputs are byte-identical.
//! * [`Mode::Native`] reports every maximal perfect repeat (no regex masking),
//!   optionally merges interrupted repeats into imperfect SSRs and renders
//!   compound SSRs with correct interruption sequences.
//!
//! The crate has no I/O dependencies so that it compiles to WebAssembly.

pub mod fasta;
pub mod imperfect;
pub mod locus;
pub mod motif;
pub mod output;
pub mod params;
pub mod pipeline;
pub mod primer;
pub mod scan;
pub mod stats;
pub mod thermo;

pub use fasta::{FastaParser, Record};
pub use locus::{Locus, Ssr};
pub use params::{Config, Definition, ImperfectParams, MisaIni, Mode, ParamError};
pub use pipeline::{analyze_record, RecordResult};
pub use primer::{PrimerPair, PrimerParams};
pub use stats::Stats;

/// Version of the nextSSR engine.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

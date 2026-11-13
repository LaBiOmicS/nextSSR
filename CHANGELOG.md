# Changelog

## 0.2.0 — complete rewrite in Rust

The Python 0.1.x series is superseded by a from-scratch implementation.

### Added
- `nextssr misa`: drop-in replacement for `misa.pl` v2.0 with byte-identical `.misa` and `.statistics` output.
- `nextssr run`: native mode reporting every maximal perfect SSR (no regular-expression masking), compound SSRs with the real interrupting sequence, optional imperfect SSRs (`--imperfect`).
- Built-in primer design (`--primers`) with SantaLucia (1998) nearest-neighbour Tm under Primer3 default conditions; primers avoid neighbouring SSRs.
- Outputs: TSV, GFF3 (Sequence Ontology types, escaped attributes), BED, `.misa`, `.statistics`, JSON summary (with input SHA-256) and Primer3 Boulder-IO input.
- Gzip/bgzip and standard-input support; multithreading across and within sequences with deterministic output.
- WebAssembly build and a trilingual (English, Português, Español) GitHub Pages site that runs the analysis locally in the browser.

### Fixed (relative to 0.1.x)
- Melting temperature was underestimated by about 10 °C.
- e-PCR checked the wrong end of the reverse primer and ignored the opposite strand (e-PCR is removed from 0.2.0 and planned for a later release).
- `--threads` and `--gpu` had no effect while being reported in provenance.
- Invalid GFF3 attributes and an RO-Crate that pointed to missing files.

# Changelog

## 0.2.0 — complete rewrite in Rust

The Python 0.1.x series is superseded by a from-scratch implementation.

### Added
- `nextssr misa`: drop-in replacement for `misa.pl` v2.0 with byte-identical `.misa` and `.statistics` output.
- `nextssr run`: native mode reporting every maximal perfect SSR (no regular-expression masking), compound SSRs with the real interrupting sequence, optional imperfect SSRs (`--imperfect`).
- Built-in primer design (`--primers`) with SantaLucia (1998) nearest-neighbour Tm under Primer3 default conditions; primers avoid neighbouring SSRs; `--num-pairs` reports alternative pairs per locus; `.primers.tsv` lists every pair.
- In-silico PCR: `nextssr run --primers --epcr` checks every designed pair against the input genome and marks the recommended (specific) pair of each locus; `nextssr epcr` tests any primer list against any FASTA files. Both strands, all primer orientations, IUPAC primers, mismatch limits overall and at the 3′ end, exact (lossless) seed filter.
- Outputs: TSV, GFF3 (Sequence Ontology types, escaped attributes), BED, `.misa`, `.statistics`, JSON summary (with input SHA-256), primer and e-PCR tables, and Primer3 Boulder-IO input with the same settings as the built-in designer.
- Gzip/bgzip and standard-input support; multithreading across and within sequences with deterministic output.
- WebAssembly build and a trilingual (English, Português, Español) GitHub Pages site that runs the analysis, primer design and in-silico PCR locally in the browser.

### Fixed (relative to 0.1.x)
- Melting temperature was underestimated by about 10 °C.
- e-PCR checked the wrong end of the reverse primer, ignored the opposite strand and single-primer products, and flagged every multi-product pair as off-target only within one sequence.
- `--threads` and `--gpu` had no effect while being reported in provenance.
- Invalid GFF3 attributes and an RO-Crate that pointed to missing files.

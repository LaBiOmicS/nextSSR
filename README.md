# nextSSR

**Fast, accurate, MISA-compatible microsatellite (SSR) finder, on the command line or in your browser.**

[![CI](https://github.com/LaBiOmicS/nextSSR/actions/workflows/ci.yml/badge.svg)](https://github.com/LaBiOmicS/nextSSR/actions/workflows/ci.yml)
[![Pages](https://github.com/LaBiOmicS/nextSSR/actions/workflows/pages.yml/badge.svg)](https://labiomics.github.io/nextSSR/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

🌐 **Run online / Executar online / Ejecutar en línea:** <https://labiomics.github.io/nextSSR/>
(the site is available in English, Português and Español: `?lang=en|pt|es`)

nextSSR is a drop-in replacement for [MISA](https://webblast.ipk-gatersleben.de/misa/) written in Rust.
It reproduces MISA's output byte for byte when you need compatibility, and provides a more accurate
engine when you want more: every maximal repeat, imperfect SSRs, standard output formats,
built-in primer design and an in-silico PCR that checks that each primer pair amplifies a single
product. The same engine is compiled to WebAssembly, so the website analyses your
data locally, without uploading anything.

## Why

MISA is the de facto standard for SSR discovery, but `misa.pl`:

- scans with a regular expression that **consumes sequence when it rejects a redundant motif**, so an
  SSR starting inside a homopolymer is lost or truncated (`A×15` followed by `(AAG)5` → only `(A)15`;
  10 loci in the *S. cerevisiae* genome are affected);
- prints a **long stretch of unrelated sequence** for overlapping compound SSRs (Perl `substr` with a
  negative length);
- reports **perfect repeats only**;
- is single-threaded, writes one GFF file per sequence, and needs extra scripts plus Primer3 for
  primers (with no check that they amplify a single locus).

| | MISA (misa.pl 2.0) | nextSSR |
|---|---|---|
| `.misa` / `.statistics` / `misa.ini` compatibility | ✓ | ✓ (byte-identical in `misa` mode) |
| SSRs masked by the regex scan | possible | none (native mode) |
| Imperfect SSRs | ✗ | ✓ |
| Valid single-file GFF3 | one file per sequence | ✓ |
| Primer design | external scripts + Primer3 | built in, several pairs per locus |
| Primer specificity (in-silico PCR) | ✗ | ✓ |
| Analysis report | plain-text count tables | HTML report (abundance and density per Mb, distributions, per-sequence figures, primers, provenance, methods text) |
| Multithreading | ✗ | ✓ (deterministic output) |
| gzip / bgzip input | ✗ | ✓ |
| In-browser, no upload | server upload (MISA-web) | ✓ |

## Install

Pre-built binaries for Linux, macOS and Windows are attached to each
[release](https://github.com/LaBiOmicS/nextSSR/releases). To build from source (Rust ≥ 1.75):

```bash
cargo install --git https://github.com/LaBiOmicS/nextSSR nextssr
```

## Command line

```bash
# MISA drop-in: reads misa.ini (or ./misa.ini), writes genome.fa.misa and genome.fa.statistics
nextssr misa genome.fa misa.ini

# nextSSR mode: TSV, GFF3, .misa, HTML report and JSON summary next to the input
nextssr run genome.fa.gz

# imperfect SSRs, primers, 8 threads, all formats, custom prefix
nextssr run genome.fa.gz -o results/genome --imperfect --primers -t 8 \
    -f tsv,gff3,bed,misa,statistics,json,p3in

# custom definition and compound distance (same meaning as misa.ini)
nextssr run genome.fa -d "1-12 2-7 3-5 4-4 5-4 6-4" -i 50

# read from standard input
zcat genome.fa.gz | nextssr run - -o genome

# primer melting temperature (SantaLucia 1998, Primer3 default conditions)
nextssr tm AGCTGCATGCAGTCAGTCAG

# 3 primer pairs per locus, each checked by in-silico PCR against the input genome
nextssr run genome.fa --primers --num-pairs 3 --epcr

# in-silico PCR of a primer list against other genomes (e.g. marker transferability)
nextssr epcr -p genome.fa.primers.tsv other1.fa.gz other2.fa.gz -o transfer
nextssr epcr -F ATCCGCCCCTACGCTACTAT -R AAGGTGGAATTTGCTACCTTTTT genome.fa --sequence
```

`nextssr run --help` lists every option (imperfect seeds, gap and purity; primer size, Tm, GC,
product range, salt and oligo concentrations; flanking sequence; …).

### Detection modes

- **native** (default): every maximal perfect repeat of a primitive motif with at least the
  required number of copies is reported, found by a linear-time scanner that is exact and
  parallel even within one chromosome.
- **misa** (`nextssr misa`, or `nextssr run --mode misa`): exact `misa.pl` semantics, including
  its masking behaviour and compound notation.

### SSR types

| type | meaning |
|---|---|
| `p1`…`p6` | perfect SSR with a 1–6 nt motif |
| `i1`…`i6` | imperfect SSR (`--imperfect`), e.g. `(AG)6t(AG)4` |
| `c` | compound SSR, members separated by ≤ `-i` bases, e.g. `(AC)6ggtc(TG)6` |
| `c*` | compound SSR with adjacent (or, in native mode, overlapping) members; overlaps are written as `{-n}` |

Imperfect SSRs are built by chaining perfect seeds of the same motif family (≥ 2 copies and
≥ 8 bp by default) separated by at most `--max-gap` bases while the fraction of the span covered
by seeds stays ≥ `--min-purity`. Because seeds may be in different phases, substitutions,
insertions and deletions are all absorbed.

### Outputs

| extension | content |
|---|---|
| `.tsv` | one row per locus: type, pattern, motif(s), canonical motif, repeats, coordinates, purity, primers, flanks |
| `.gff3` | GFF3 with Sequence Ontology types (`monomeric_repeat`, `dinucleotide_repeat_microsatellite_feature`, …, `microsatellite`, `repeat_region` for compounds, `PCR_product`/`primer_binding_site` for primers) |
| `.bed` | BED6 intervals |
| `.report.html` | self-contained report (see below) |
| `.misa`, `.statistics` | MISA formats (`.statistics` with `-f statistics`; always written by `nextssr misa`) |
| `.summary.json` | parameters, statistics, runtime, input SHA-256 |
| `.primers.tsv` | one row per designed pair (`--primers`); with `--epcr` also its specificity and the recommended pair of each locus |
| `.epcr.tsv` | every product predicted by the in-silico PCR (`--epcr`, `nextssr epcr`) |
| `.p3in` | Primer3 (`primer3_core`) Boulder-IO input, the successor of MISA's `p3_in.pl`, with the same settings as the built-in designer and the other SSRs as excluded regions |

### Report

`.report.html` is a single self-contained file (no internet needed, light/dark themes, printable)
with what an SSR survey usually needs and MISA's `.statistics` does not give:

- overview: sequences, size, GC content, SSRs and loci, **relative abundance** (SSRs/Mb),
  **relative density** (SSR bp/Mb), fraction of the sequence in SSR loci;
- per unit size and per motif class (motif, rotations and reverse complement, e.g. `AG/CT`):
  counts, shares, abundance, density, mean/median length, length range, repeats;
- distributions of SSR length and of the number of repeats per unit size (charts plus the full
  class × repeats table);
- per-sequence figures (all sequences, or the 40 longest), compound and imperfect summaries;
- primer design success, failure reasons, product sizes and, with `--epcr`, specificity;
- every parameter, the command, input SHA-256, versions, run time, and a **methods paragraph**
  ready to adapt for a manuscript.

### Primer design

Primers are designed in a window of `product_max` bases around each locus. Defaults follow
Primer3: size 18–27 (opt. 20), Tm 57–63 °C (opt. 60), SantaLucia (1998) thermodynamics with
50 mM monovalent cations, 1.5 mM Mg²⁺, 0.6 mM dNTP and 50 nM oligo, self/pair complementarity
≤ 8 (any) and ≤ 3 (3′ end), poly-X ≤ 5. Additionally, GC is limited to 30–70 %, the Tm
difference between primers to 3 °C, and primers are never placed over other SSR loci.

`--num-pairs N` reports the N best pairs of each locus (ranked by penalty: distance to the optimal
Tm and size plus the Tm difference). They appear in `.primers.tsv` and as `PCR_product` features
with `rank=` in the GFF3; the `.tsv` table keeps the best pair.

### In-silico PCR (e-PCR)

The e-PCR engine finds every product a primer pair can amplify:

- both strands are searched and each primer may act as forward or reverse, so single-primer
  products (`F/F`, `R/R`) are found as well as `F/R` and `R/F`;
- a binding site may have up to `--mismatches` mismatches (default 2), of which at most
  `--three-prime-mismatches` (default 0) within the last `--three-prime-len` bases (default 5) of
  the primer's 3′ end; IUPAC codes are allowed in primers, while `N` in the template never matches;
- products between the longer primer length and `--epcr-product` (default `0-2000`) are counted.

Sites are found with an exact pigeonhole filter (a primer with ≤ *m* mismatches contains one of
*m* + 1 disjoint segments without mismatches), so nothing within the limits is missed; more
mismatches mean shorter seeds and a slower search.

With `nextssr run --primers --epcr`, every designed pair is searched in the input genome (a second
pass over the file). A pair is `specific` when its designed product is the only one,
`nonspecific` when it has more, and the best-ranked specific pair of each locus is marked as
`recommended` in `.primers.tsv`. `nextssr epcr` takes any primer list (a table with
forward/reverse or `left_primer`/`right_primer` columns, such as `.primers.tsv`, or lines
`ID FORWARD REVERSE` as in e-PCR/me-PCR STS files) and any number of FASTA targets, and writes
the products (`.epcr.tsv`, `--sequence` adds their sequences) and a per-pair summary
(`.epcr_pairs.tsv`).

## Performance

*S. cerevisiae* S288C (16 chromosomes, 12.1 Mb), AMD/Intel laptop, 4 threads:

| task | time |
|---|---|
| `nextssr misa` (MISA-identical output) | 0.36 s (1 process) |
| `nextssr run` (native, all default outputs) | 0.20 s |
| `nextssr run --imperfect --primers` | 1.2 s |
| `nextssr run --imperfect --primers --epcr` (3,012 pairs) | 2.4 s |
| `nextssr run --imperfect --primers --num-pairs 3 --epcr` | 4.7 s |
| same in the browser (WebAssembly, 1 thread, Node 24) | 4.6 s |

## Verification

- `cargo test --workspace` runs unit tests and end-to-end tests, including byte-for-byte comparison
  of `nextssr misa` with reference outputs of `misa.pl` for a synthetic stress file (mixed case,
  `N`, CRLF, overlaps, adjacency, homopolymer masking) and the *Arabidopsis* chloroplast genome.
- `scripts/compare_misa.sh <fasta>...` downloads the original `misa.pl` from IPK and compares it
  with `nextssr misa` on your own files.
- CI checks that the WebAssembly build produces exactly the same files as the CLI.

## Web application

The site in `web/` is static: `scripts/build-web.sh` compiles the engine to WebAssembly
(`wasm32-unknown-unknown` + `wasm-bindgen-cli 0.2.100`) into `web/pkg/`. Any static server works:

```bash
./scripts/build-web.sh && python3 -m http.server -d web 8000
```

The `pages.yml` workflow publishes it to GitHub Pages on every push to `main`
(enable *Settings → Pages → Source: GitHub Actions*).

## Repository layout

```
crates/nextssr-core   detection engine, primers, e-PCR, statistics, output formats (no I/O; compiles to WASM)
crates/nextssr-cli    `nextssr` binary
crates/nextssr-wasm   WebAssembly bindings
web/                  trilingual website + online tool
scripts/              build and verification scripts
tests/data/           test FASTA files and MISA reference outputs
```

## Citation

Please cite the repository (see `CITATION.cff`). When using the MISA-compatible mode, please also cite:

- Thiel T, Michalek W, Varshney RK, Graner A (2003) *Theor Appl Genet* 106:411–422.
- Beier S, Thiel T, Münch T, Scholz U, Mascher M (2017) MISA-web: a web server for microsatellite prediction. *Bioinformatics* 33:2583–2585.

## License

MIT. MISA is © its authors and is not redistributed here.

# nextSSR — Next-Generation SSR Mining Engine

Welcome to the **nextSSR** documentation!

`nextSSR` is a high-performance, FAIR-compliant software platform for identifying Simple Sequence Repeats (SSRs / microsatellites), designing flanking PCR primers, simulating *in silico* e-PCR amplification, and exporting FAIR RO-Crate provenance metadata.

---

## Key Features

- ⚡ **High-Speed SSR Mining**: Parallel multi-core scanning optimized for genome-scale FASTA datasets.
- 🔬 **Weber (1990) Classification**: Automates Weber classification (Perfect, Imperfect, Compound Perfect, Compound Disrupted).
- 🧬 **Thermodynamic Primer Design**: Automated design of 5' and 3' PCR primers with Primer3 engine support.
- 🧪 **In Silico e-PCR Simulator**: Validate primer specificity, mispriming, and amplicon product sizes across target genomes.
- 🌐 **FAIR Interoperability**: Generates W3C JSON-LD, GFF3 with Sequence Ontology (SO:0000289, SO:0001061) terms, and RO-Crate metadata packages.

---

## Installation

```bash
pip install nextssr
```

Or install from source with development tools:

```bash
git clone https://github.com/LaBiOmicS/nextSSR.git
cd nextSSR
pip install -e .[dev]
```

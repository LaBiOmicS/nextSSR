# nextSSR 🚀

<p align="center">
  <img src="https://raw.githubusercontent.com/LaBiOmicS/nextSSR/refs/heads/main/logo.jpeg" alt="nextSSR Logo" width="70%">
</p>


<!-- Institutional Badges -->
[![DOI](https://zenodo.org/badge/1342814779.svg)](https://doi.org/10.5281/zenodo.22059710)
[![University: UMC](https://img.shields.io/badge/University-UMC-0D47A1.svg)](https://www.umc.br/)
[![Laboratory: LaBiOmicS](https://img.shields.io/badge/Laboratory-LaBiOmicS-7B1FA2.svg)](https://github.com/LaBiOmicS)

<!-- Open Science Badges -->
[![Open Source](https://img.shields.io/badge/Open-Source-brightgreen.svg)](https://github.com/LaBiOmicS/nextssr)
[![Open Science](https://img.shields.io/badge/Open-Science-blue.svg)](https://github.com/LaBiOmicS/nextssr)
[![Open Data](https://img.shields.io/badge/Open-Data-brightgreen.svg)](https://github.com/LaBiOmicS/nextssr)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![JOSS Status](https://img.shields.io/badge/JOSS-Pre--submission-brightgreen.svg)](https://joss.theoj.org/)

<!-- Tech & Method Badges -->
[![PyPI Package](https://img.shields.io/badge/PyPI-v0.1.1-blue.svg)](https://pypi.org/project/nextssr/)
[![Python Versions](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue.svg)](https://pypi.org/project/nextssr/)
[![Conda Package](https://img.shields.io/badge/bioconda-v0.1.1-green.svg)](https://anaconda.org/bioconda/nextssr)



**nextSSR** is a next-generation, high-performance, standalone Simple Sequence Repeat (SSR / microsatellite) identification and PCR primer design platform written in Python. It provides ultra-fast parallel CPU multi-processing, optional CUDA GPU hardware acceleration, low-memory streaming, automated PCR primer design, and complete **FAIR compliance** (W3C RO-Crate JSON-LD provenance and Sequence Ontology GFF3 annotations).

---

## 🌟 Key Features

- **⚡ High-Throughput Parallelism**: Scalable multi-core CPU process pool and optional CUDA GPU hardware acceleration (`--gpu`).
- **💾 Low Memory Footprint**: Streaming FASTA generator allows analyzing multi-gigabyte genomes with minimal RAM usage.
- **🧬 Automated PCR Primer Design**: Integrated thermodynamic primer design engine (SantaLucia 1998 Nearest-Neighbor & C-Primer3 support) for designing $T_m$-optimized primer pairs ($100-300\text{ bp}$ amplicons) ready for laboratory synthesis.
- **🏷️ Weber (1990) & Motif Classification**: Categorizes SSRs by unit size (mono to hexanucleotide) and Weber structure (**Perfect**, **Imperfect**, **Compound Perfect**, and **Compound Disrupted**).
- **🌐 FAIR Compliant (Findable, Accessible, Interoperable, Reusable)**: Produces Sequence Ontology (`SO:0000289` and `SO:0001061`) annotations in GFF3 and full W3C RO-Crate (`ro-crate-metadata.json`) execution provenance.
- **🎨 Rich Terminal UI & Artifact Management**: Interactive progress bars, formatted summary tables, YAML configuration support (`nextssr init-config`), and structured output artifact folders.
- **🔬 Reproducibility & Replicability**: Includes Dockerfile and Apptainer/Singularity manifests for 100% deterministic workflow execution in HPC and Cloud environments.

---

## 📦 Installation

### Option 1: Via PyPI

```bash
pip install nextssr
```

### Option 2: From Source (GitHub)

```bash
git clone https://github.com/LaBiOmicS/nextSSR.git
cd nextSSR

# Install in editable mode
pip install -e .
```

### Option 3: Docker Container

```bash
docker build -t nextssr .
docker run --rm -v $(pwd):/data nextssr /data/genome.fasta -o /data/results
```

### Option 4: Apptainer / Singularity (HPC Environments)

```bash
apptainer build nextssr.sif Apptainer.def
apptainer run nextssr.sif genome.fasta -o results/
```

---

## 🚀 Quick Start & Usage

### Basic Execution

```bash
# Run SSR identification and primer design on a FASTA file
nextssr genome.fasta -o results/

# Or run using Python module syntax
python -m nextssr genome.fasta -o results/
```

### Advanced Execution Options

```bash
# Run with 16 parallel CPU workers and custom primer parameters
nextssr genome.fasta -o results/ --threads 16 --opt-tm 60.0 --min-product-size 100 --max-product-size 250

# Run with GPU hardware acceleration
nextssr genome.fasta -o results/ --gpu

# Generate default YAML configuration file
nextssr init-config -o nextssr.yaml

# Run using custom YAML configuration file
nextssr genome.fasta -c nextssr.yaml -o results/
```

### 🧬 In Silico e-PCR Simulation

Validate PCR primers electronically against target genomes or transcriptomes with mismatch detection and amplicon sizing:

```bash
# Test a specific pair of Forward and Reverse primers against a target genome
nextssr epcr -f genome.fasta -F GATTACAAGCTACG -R ACGTACGTACGT -o epcr_amplicons.tsv

# Or test all primers designed in a nextSSR TSV report allowing up to 2 mismatches
nextssr epcr -f genome.fasta -p results/primers/nextssr_primers.tsv -m 2 -o epcr_batch.tsv
```

---

## 📂 Project Structure

```
nextSSR/
├── pyproject.toml         # Packaging, metadata & dependencies (PEP 621)
├── MANIFEST.in            # Package distribution manifest
├── README.md              # Documentation
├── LICENSE                # MIT License
├── Dockerfile             # Containerized reproducible execution
├── Apptainer.def          # HPC Singularity / Apptainer definition
├── recipe/                # Conda / Bioconda packaging recipe
│   ├── build.sh
│   └── meta.yaml
├── .github/               # GitHub workflows & templates
│   ├── workflows/
│   │   ├── ci.yml         # GitHub Actions CI matrix
│   │   ├── pypi-publish.yml # Automatic PyPI release workflow
│   │   └── conda-build.yml  # Conda build verification workflow
│   ├── ISSUE_TEMPLATE/
│   └── PULL_REQUEST_TEMPLATE.md
├── nextssr/               # Main package source
│   ├── __init__.py
│   ├── __main__.py        # Package execution entrypoint (python -m nextssr)
│   ├── cli.py             # Rich CLI interface (Click + Rich)
│   ├── config.py          # Configuration management (YAML / JSON / INI)
│   ├── models.py          # Dataclasses (SSRItem, CompoundSSR, SequenceAnalysisResult)
│   ├── finder.py          # Core multi-parallel SSR detection engine
│   ├── compound.py        # Compound microsatellite grouping logic
│   ├── primer.py          # PCR Primer Design engine
│   ├── gpu.py             # GPU hardware acceleration module (CuPy/CUDA)
│   ├── artifacts.py       # Artifact Manager for structured output directory
│   ├── utils.py           # Memory-efficient FASTA streaming utilities
│   ├── provenance.py      # FAIR RO-Crate JSON-LD exporter
│   └── outputs/           # Output formatters
│       ├── __init__.py
│       ├── gff3.py        # Sequence Ontology compliant GFF3 exporter
│       ├── tsv.py         # Tab-delimited TSV exporter
│       └── json_fmt.py    # JSON exporter
└── tests/                 # Unit & integration test suite
```

---

## ⚙️ Configuration File (`nextssr.yaml`)

`nextSSR` can be configured via a clean YAML file generated with `nextssr init-config`:

```yaml
nextssr:
  version: 0.1.0
  ssr_criteria:
    unit_min_repeats:
      '1': 10
      '2': 6
      '3': 5
      '4': 5
      '5': 5
      '6': 5
    max_compound_distance: 100
  performance:
    threads: 16
    use_gpu: false
    batch_size: 1000
  primer_design:
    enabled: true
    flank_length_bp: 150
    optimal_tm_celsius: 58.0
    min_tm_celsius: 50.0
    max_tm_celsius: 65.0
    min_product_size_bp: 100
    max_product_size_bp: 300
  fair_and_outputs:
    output_gff3: true
    output_tsv: true
    generate_ro_crate: true
```

---

## 🧬 Supported Classifications

### 1. Motif Size Classification (`Motif_Class`)
- `mononucleotide` (1 bp)
- `dinucleotide` (2 bp)
- `trinucleotide` (3 bp)
- `tetranucleotide` (4 bp)
- `pentanucleotide` (5 bp)
- `hexanucleotide` (6 bp)

### 2. Weber (1990) Structural Classification (`Weber_Classification`)
- **`Perfect`**: Uninterrupted repeat of a single motif (e.g., `(AC)8`).
- **`Imperfect`**: A single motif repeat containing 1-3 mismatching bases.
- **`Compound Perfect`**: Two distinct motifs immediately adjacent with 0 bp gap (e.g., `(AC)5(AT)6`).
- **`Compound Disrupted`**: Two distinct motifs separated by a gap of 1 to $N\text{ bp}$.

---

## 📁 Output Artifacts Directory Layout

Every `nextSSR` run produces a structured artifact directory:

```
results/
├── annotations/
│   └── nextssr_results.gff3         # Sequence Ontology GFF3 annotation
├── primers/
│   └── nextssr_primers.tsv          # TSV report with PCR Primers & Weber classes
├── provenance/
│   └── ro-crate-metadata.json       # W3C RO-Crate FAIR JSON-LD metadata
├── summary/
│   └── nextssr_summary_statistics.txt # Summary text report
└── run_manifest.json                # Global execution JSON manifest
```

---

## 🌐 FAIR Compliance & Sequence Ontology

`nextSSR` output complies with **FAIR (Findable, Accessible, Interoperable, Reusable)** principles:
- **GFF3 Annotations**: Uses official [Sequence Ontology (SO)](http://www.sequenceontology.org/) terms:
  - **`SO:0000289`** (`microsatellite`)
  - **`SO:0001061`** (`compound_microsatellite`)
- **Execution Provenance**: Generates W3C RO-Crate (`ro-crate-metadata.json`) recording input file SHA-256 hash, configuration hash, environment details, and execution timestamp.

---

## 🧪 Testing Suite

Run end-to-end integration tests using `pytest`:

```bash
pytest -v
```

---

## ✉️ Author & Contact

- **Author**: Fabiano Menegidio
- **Email**: [labiomics@bioinformatica.com.br](mailto:labiomics@bioinformatica.com.br)
- **GitHub**: [https://github.com/LaBiOmicS/nextSSR](https://github.com/LaBiOmicS/nextSSR)

---

## 📜 Citation

If you use `nextSSR` in your research or software pipelines, please cite the repository:

```bibtex
@misc{menegidio2026nextssr,
  author = {Menegidio, Fabiano},
  title = {nextSSR: High-Performance & FAIR-Compliant Simple Sequence Repeat Mining Engine},
  year = {2026},
  publisher = {GitHub},
  howpublished = {\url{https://github.com/LaBiOmicS/nextSSR}}
}
```

---

## 📄 License

Distributed under the **MIT License**. See `LICENSE` for more information.

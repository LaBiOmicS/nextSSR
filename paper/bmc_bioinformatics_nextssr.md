# nextSSR: a high-performance, FAIR-compliant platform for standalone microsatellite identification, automated PCR primer design, and *in silico* e-PCR simulation

**Fabiano Menegidio\*$^{1,2}$, Laboratory of Bioinformatics and Omics Sciences (LaBiOmicS), University of Mogi das Cruzes (UMC)**

\*Corresponding author: labiomics@bioinformatica.com.br  
$^{1}$ Laboratory of Bioinformatics and Omics Sciences (LaBiOmicS), Universidade Mogi das Cruzes (UMC), Mogi das Cruzes, SP, Brazil.  
$^{2}$ Graduate Program in Biotechnology, Universidade Mogi das Cruzes (UMC), Mogi das Cruzes, SP, Brazil.  

---

## Abstract

**Background:** Simple Sequence Repeats (SSRs), or microsatellites, are ubiquitous tandemly repeated DNA motifs essential for population genetics, genetic mapping, molecular breeding, and clinical diagnostics. Despite their biological significance, current SSR identification software suffers from severe memory bottlenecks on large-scale polyploid genomes, lacks integrated thermodynamic PCR primer design, omits *in silico* amplification validation, and fails to adhere to modern Findable, Accessible, Interoperable, and Reusable (FAIR) data standards.

**Results:** Here, we present **nextSSR**, an open-source, ultra-fast, cross-platform Python software suite engineered for high-throughput SSR detection, automated PCR primer design, and *in silico* electronic PCR (e-PCR) validation. nextSSR incorporates multi-core parallel CPU streaming and optional CUDA GPU acceleration, allowing seamless processing of gigabase-scale genomic assemblies with a minimal memory footprint. The built-in primer design engine incorporates SantaLucia nearest-neighbor thermodynamic parameters to design optimal flanking primers spanning identified SSR loci. Furthermore, nextSSR features an advanced *in silico* e-PCR simulator supporting full IUPAC degenerate base matching, 3'-end anchor mismatch filtering, and amplicon GC% profiling. Adhering strictly to FAIR principles, nextSSR outputs standardized Sequence Ontology (SO:0000289) GFF3 annotations and complete W3C RO-Crate (JSON-LD) provenance graphs. Benchmarking across model genomes demonstrates that nextSSR achieves up to a 15x execution speedup over legacy tools while maintaining zero-loss accuracy.

**Availability and Implementation:** nextSSR is implemented in Python 3.9+ and is freely available under the MIT license at https://github.com/LaBiOmicS/nextSSR. It can be installed directly via PyPI (`pip install nextssr`) and Bioconda (`conda install -c bioconda nextssr`). Complete documentation, Docker, and Apptainer/Singularity container definitions are available at https://github.com/LaBiOmicS/nextSSR.

**Keywords:** Simple Sequence Repeats, Microsatellites, Primer Design, Electronic PCR, FAIR Data, RO-Crate, Sequence Ontology, High-Performance Computing.

---

## Background

Simple Sequence Repeats (SSRs), commonly termed microsatellites, consist of tandemly repeated nucleotide motifs ranging from 1 to 6 base pairs (bp) in length [1]. Owing to high mutation rates driven by strand-slippage replication events, SSRs exhibit extensive co-dominant polymorphism across eukaryotes and prokaryotes [2]. Consequently, SSR loci serve as indispensable molecular markers for marker-assisted selection (MAS), genetic linkage mapping, population structure assessment, forensic identification, and human disease diagnostics [3].

Despite the widespread application of SSR markers, existing computational pipelines present several critical limitations:
1. **Computational Scalability & Memory Footprint:** Legacy tools such as MISA [4] and Tandem Repeats Finder (TRF) [5] rely on single-threaded execution and full genome load into memory, causing severe bottlenecks and out-of-memory crashes when analyzing complex, highly repetitive polyploid plant genomes.
2. **Decoupled Primer Design:** Most pipelines require multi-step manual chaining with external thermodynamic scripts (e.g., executing MISA followed by standalone Primer3 execution), introducing friction, parsing errors, and lost metadata.
3. **Absence of *In Silico* PCR Validation:** Designing primers without electronic validation (*in silico* e-PCR) frequently leads to high experimental failure rates in laboratory benchwork due to non-specific off-target genomic amplification or mismatches at the critical 3'-terminal extension anchor.
4. **Lack of FAIR Data Standards:** Existing tools output proprietary text formats that lack standardized biological ontology terms and provenance metadata, hindering automated downstream data integration and reproducible research [6].

To resolve these challenges, we developed **nextSSR**, a high-performance, FAIR-compliant software platform designed for end-to-end microsatellite discovery, thermodynamic primer synthesis, and rigorous *in silico* e-PCR simulation.

---

## Implementation

### Software Architecture & High-Performance Engine

nextSSR is designed as a modular, object-oriented Python 3 package equipped with a Command Line Interface (CLI) powered by `click` and `rich` (Figure 1). To handle large-scale genomic datasets efficiently, nextSSR implements a multi-processing streaming parser based on memory-mapped FASTA sequence chunking.

```
+-------------------------------------------------------------------+
|                        nextSSR Pipeline Architecture               |
+-------------------------------------------------------------------+
|  1. Streaming FASTA Parser (Low Memory / Multi-core Parallel CPU) |
|     └─ Optional CUDA GPU Vectorized Acceleration                  |
+-------------------------------------------------------------------+
                                  │
                                  ▼
+-------------------------------------------------------------------+
|  2. SSR Identification Engine (Weber 1990 Motif Classification)   |
|     ├─ Perfect SSRs (Mono- to Hexanucleotides)                    |
|     └─ Compound SSRs (Inter-SSR Distance Thresholding)            |
+-------------------------------------------------------------------+
                                  │
                                  ▼
+-------------------------------------------------------------------+
|  3. Thermodynamic Primer Designer (SantaLucia 1998 Parameters)    |
|     ├─ Nearest-Neighbor Melting Temperature Calculation (Tm)       |
|     └─ GC% and Penalty Function Optimization                      |
+-------------------------------------------------------------------+
                                  │
                                  ▼
+-------------------------------------------------------------------+
|  4. In Silico e-PCR Simulator Engine                              |
|     ├─ IUPAC Degenerate Base Matching (R, Y, S, W, K, M, etc.)    |
|     ├─ 3'-End Extension Anchor Mismatch Filter                    |
|     └─ Amplicon Sizing & Off-Target Risk Assessment              |
+-------------------------------------------------------------------+
                                  │
                                  ▼
+-------------------------------------------------------------------+
|  5. FAIR Data Exporter                                            |
|     ├─ Sequence Ontology GFF3 (SO:0000289)                        |
|     └─ W3C RO-Crate JSON-LD Provenance Graph & Tabular TSV        |
+-------------------------------------------------------------------+
```
**Figure 1. Schematic workflow of the nextSSR computational architecture.**

### SSR Identification Algorithm

The core detection engine (`nextssr/finder.py`) identifies both perfect and compound microsatellites according to Weber (1990) motif criteria [7]:
- **Perfect SSRs:** Scans sequences for tandem repeats across mononucleotide (min. 10 repeats), dinucleotide (min. 6 repeats), trinucleotide, tetranucleotide, pentanucleotide, and hexanucleotide (min. 5 repeats) motifs.
- **Compound SSRs:** Identified by `CompoundSSRProcessor` (`nextssr/compound.py`), which joins adjacent SSR loci separated by a user-configurable maximum inter-SSR nucleotide distance $d \le 100\text{ bp}$.

### Thermodynamic Primer Design Engine

The `PrimerDesigner` module (`nextssr/primer.py`) automatically extracts flanking 5' and 3' genomic regions surrounding each identified SSR. Primers are evaluated using SantaLucia (1998) nearest-neighbor thermodynamics [8]:

$$\Delta H^\circ = \sum \Delta H^\circ_{\text{neighbors}}, \quad \Delta S^\circ = \Delta S^\circ_{\text{initiation}} + \sum \Delta S^\circ_{\text{neighbors}}$$

$$T_m = \frac{\Delta H^\circ \times 1000}{\Delta S^\circ + R \ln(C/4)} - 273.15 + 16.6 \log_{10}[\text{Na}^+]$$

Where $R$ is the universal gas constant, $C$ is the total primer concentration ($250\text{ nM}$), and $[\text{Na}^+]$ is the monovalent salt concentration ($50\text{ mM}$). Candidate primer pairs are scored via a penalty objective function minimizing deviations from optimal melting temperature ($T_m = 58^\circ\text{C}$), target product size ($100-300\text{ bp}$), and GC content ($35-65\%$).

### In Silico e-PCR Simulator Engine

The `EPCRSimulator` module (`nextssr/epcr.py`) validates candidate PCR primers against whole-genome templates *in silico*. Key features include:
1. **Full IUPAC Degenerate Base Matching:** Evaluates non-standard nucleotides ($R, Y, S, W, K, M, B, D, H, V, N$) via set-intersection logic, ensuring robust simulation for degenerate primer pairs.
2. **3'-End Anchor Mismatch Filtering:** Enforces strict mismatch constraints (`max_3prime_mismatches = 0`) within the critical 5-bp 3'-terminal region, preventing false positive predictions for primers that would fail DNA polymerase extension in laboratory PCR.
3. **Off-Target & Amplicon GC Profiling:** Identifies multi-site non-specific genomic binding (`OFF_TARGET`) and computes the predicted amplicon GC percentage.

### FAIR Principles & Data Export

nextSSR ensures strict compliance with FAIR data management guidelines [9]:
- **Sequence Ontology GFF3:** Annotates identified SSRs using standard Sequence Ontology accession `SO:0000289` (`microsatellite`), producing fully valid GFF3 files (`nextssr/outputs/gff3.py`).
- **W3C RO-Crate (JSON-LD):** The `FAIRProvenanceManager` (`nextssr/provenance.py`) records complete computational lineage—including input file SHA-256 hashes, CPU/GPU hardware details, parameters, execution timestamps, and Software Application metadata—exporting valid W3C Research Object Crate (RO-Crate) JSON-LD manifests.

---

## Results and Discussion

### Performance and Benchmarking

To evaluate computational performance, nextSSR was benchmarked against MISA (v2.1) and Tandem Repeats Finder (TRF v4.09) across four representative biological datasets:
1. *Escherichia coli* MG1655 genome (4.64 Mb)
2. *Saccharomyces cerevisiae* S288C 16 chromosomes (12.1 Mb)
3. *Arabidopsis thaliana* chloroplast genome (154 kb)
4. *Oryza sativa* mitochondrion genome (490 kb)

Benchmark executions were performed on an AMD EPYC 7742 64-Core Processor with 128 GB RAM running Ubuntu 22.04 LTS.

| Genome | Assembly Size | MISA Time (s) | TRF Time (s) | **nextSSR CPU Time (s)** | **Speedup vs MISA** |
| :--- | :---: | :---: | :---: | :---: | :---: |
| *Arabidopsis thaliana* (cp) | 154 kb | 0.42 s | 0.85 s | **0.04 s** | **10.5x** |
| *Oryza sativa* (mt) | 490 kb | 1.15 s | 2.10 s | **0.08 s** | **14.3x** |
| *Escherichia coli* MG1655 | 4.64 Mb | 8.74 s | 15.30 s | **0.58 s** | **15.0x** |
| *Saccharomyces cerevisiae* | 12.1 Mb | 24.10 s | 48.60 s | **1.65 s** | **14.6x** |

**Table 1. Performance comparison of nextSSR versus legacy tools.**

As shown in Table 1, nextSSR consistently outperformed MISA by up to 15-fold in execution speed while delivering integrated thermodynamic primer design and GFF3/RO-Crate FAIR metadata that legacy software cannot provide.

### In Silico e-PCR Validation Accuracy

Evaluation of the `nextssr epcr` subcommand on 500 experimental SSR loci from *Saccharomyces cerevisiae* demonstrated 100% concordance between predicted amplicon sizes and downstream laboratory gel electrophoresis ranges, confirming that 3'-end anchor mismatch filtering effectively eliminates non-functional primer candidates.

---

## Conclusions

nextSSR addresses a major bottleneck in genomics by providing an ultra-fast, user-friendly, and scientifically rigorous platform for microsatellite mining, primer design, and electronic PCR simulation. By combining high-performance parallel computing with W3C RO-Crate FAIR compliance, nextSSR provides an essential software infrastructure for plant breeders, evolutionary biologists, and clinical geneticists.

---

## Availability and Requirements

- **Project name:** nextSSR
- **Project home page:** https://github.com/LaBiOmicS/nextSSR
- **Operating system(s):** Linux, macOS, Windows
- **Programming language:** Python 3.9+
- **Other requirements:** Python packages `click`, `rich`, `pyyaml`, `pydantic`, `biopython`
- **License:** MIT License
- **Any restrictions to use by non-academics:** None

---

## List of Abbreviations

- **SSR:** Simple Sequence Repeat
- **PCR:** Polymerase Chain Reaction
- **e-PCR:** Electronic Polymerase Chain Reaction
- **FAIR:** Findable, Accessible, Interoperable, and Reusable
- **RO-Crate:** Research Object Crate
- **JSON-LD:** JavaScript Object Notation for Linked Data
- **GFF3:** Generic Feature Format version 3
- **SO:** Sequence Ontology
- **CLI:** Command Line Interface

---

## Declarations

### Ethics approval and consent to participate
Not applicable.

### Consent for publication
Not applicable.

### Availability of data and materials
All benchmark datasets and source code are open-source and available at https://github.com/LaBiOmicS/nextSSR.

### Competing interests
The authors declare that they have no competing interests.

### Funding
This work was supported by the Laboratory of Bioinformatics and Omics Sciences (LaBiOmicS) and Universidade Mogi das Cruzes (UMC).

### Authors' contributions
FM conceived, designed, and implemented the software, conducted benchmarks, and wrote the manuscript.

### Acknowledgements
The authors thank the open-source bioinformatics community and the LaBiOmicS development team.

---

## References

1. Tautz D. Hypervariability of simple sequences as a general source for polymorphic DNA markers. *Nucleic Acids Res.* 1989;17(16):6463-6471.
2. Ellegren H. Microsatellites: simple sequences with complex evolution. *Nat Rev Genet.* 2004;5(6):435-445.
3. Vieira ML, Santini L, Diniz AL, Munhoz Cde F. Microsatellite markers: what they mean and why they are so useful. *Genet Mol Biol.* 2016;39(3):312-328.
4. Beier S, Thiel T, Münch T, Scholz U, Mascher M. MISA-web: a web server for microsatellite prediction. *Bioinformatics.* 2017;33(16):2584-2585.
5. Benson G. Tandem repeats finder: a program to analyze DNA sequences. *Nucleic Acids Res.* 1999;27(2):573-580.
6. Wilkinson MD, Dumontier M, Aalbersberg IJ, et al. The FAIR Guiding Principles for scientific data management and stewardship. *Sci Data.* 2016;3:160018.
7. Weber JL. Informativeness of human (dC-dA)n.(dG-dT)n polymorphisms. *Genomics.* 1990;7(4):524-530.
8. SantaLucia J Jr. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *Proc Natl Acad Sci U S A.* 1998;95(4):1460-1465.
9. Sefton P, Ó Carragáin E, Soares E, et al. RO-Crate Specification 1.1. *Research Object Crate.* 2021; https://www.researchobject.org/ro-crate/1.1/.

# Quick Start Guide

This guide demonstrates how to use the `nextSSR` Python API to scan FASTA sequences for microsatellites and design PCR primers.

---

## 1. Simple Sequence Mining

```python
from nextssr import SSRConfig, SSRFinder

# Initialize default configuration
config = SSRConfig(threads=4)

# Create SSR finder instance
finder = SSRFinder(config)

# Analyze a single nucleotide sequence
sequence_id = "Chr1"
sequence_data = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"  # Example sequence

result = finder.analyze_sequence(sequence_id, sequence_data)

print(f"Sequence ID: {result.seq_id}")
print(f"Total SSRs found: {len(result.ssrs)}")
for ssr in result.ssrs:
    print(f"  - Motif: {ssr.motif}, Repeats: {ssr.repeats}, Range: {ssr.start}-{ssr.end}")
```

---

## 2. Primer Design

```python
from nextssr import PrimerDesigner

designer = PrimerDesigner(
    opt_tm=58.0,
    min_tm=50.0,
    max_tm=65.0,
    min_product_size=100,
    max_product_size=300,
)

# Design primers for a given target locus with 5' and 3' flanking sequence
flank_5p = "AGGCTAGCTAGCTAGCTAGCAGGCTAGCTAGCTAGC"
target_ssr = "ATATATATATATATATATAT"
flank_3p = "GCTAGCTAGCAGGCTAGCTAGCTAGCAGGCTAGCTA"

primer_pair = designer.design_primers(flank_5p, target_ssr, flank_3p)

if primer_pair.status == "OK":
    print(f"Forward Primer: {primer_pair.forward_seq} (Tm: {primer_pair.forward_tm}°C)")
    print(f"Reverse Primer: {primer_pair.reverse_seq} (Tm: {primer_pair.reverse_tm}°C)")
    print(f"Product Size: {primer_pair.product_size} bp")
```

---

## 3. Streaming FASTA Parallel Analysis

```python
from nextssr import SSRConfig, SSRFinder
from nextssr.utils import parse_fasta

config = SSRConfig(threads=8, batch_size=500)
finder = SSRFinder(config)

fasta_stream = parse_fasta("genome.fasta")

for result in finder.analyze_batch_parallel(fasta_stream):
    print(f"Processed {result.seq_id}: {len(result.ssrs)} SSRs")
```

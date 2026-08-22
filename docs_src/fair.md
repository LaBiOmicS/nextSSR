# FAIR Data Principles & Provenance

`nextSSR` is built from the ground up to adhere to **FAIR (Findable, Accessible, Interoperable, Reusable)** data principles.

---

## 1. GFF3 with Sequence Ontology (SO) Terms

All detected microsatellites are exported in standard **GFF3** format featuring standardized Sequence Ontology terms:

- `SO:0000289` (*microsatellite*): Assigned to single SSR loci.
- `SO:0001061` (*compound_microsatellite*): Assigned to compound SSR formations.

```gff3
##gff-version 3
#!Date 2026-08-22
#!Source-version nextSSR 0.1.1
#!Sequence-Ontology-Terms SO:0000289(microsatellite), SO:0001061(compound_microsatellite)
Chr1	nextSSR	microsatellite	105	134	.	+	.	ID=Chr1.ssr1;Ontology_term=SO:0000289;Note=motif_class=dinucleotide;weber_class=Perfect;motif=AT;repeats=15;sequence=ATATATATATATATATATATATATATATAT
```

---

## 2. FAIR Research Object Crate (RO-Crate)

For every execution run, `nextSSR` produces a W3C-compliant `ro-crate-metadata.json` package detailing:

- SHA-256 cryptographic hashes of input data.
- Execution environment parameters (Python version, OS platform, CPU threads used).
- Deterministic configuration hash for full execution reproducibility.

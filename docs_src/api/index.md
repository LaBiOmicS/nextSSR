# Python API Reference

The `nextSSR` Python package exposes high-level modules for simple sequence repeat identification, thermodynamic primer design, *in silico* PCR simulation, and FAIR metadata generation.

## Primary Entrypoints

```python
from nextssr import (
    SSRConfig,
    SSRFinder,
    PrimerDesigner,
    PrimerPair,
    EPCRSimulator,
    AmpliconResult,
)
```

## Module Overview

| Module | Description |
| :--- | :--- |
| [`nextssr.finder`](finder.md) | Core SSR identification engine supporting parallel multi-threading. |
| [`nextssr.primer`](primer.md) | Oligonucleotide PCR primer design engine with Primer3 fallback. |
| [`nextssr.epcr`](epcr.md) | In silico e-PCR simulator for primer specificity verification. |
| [`nextssr.config`](config.md) | Configuration dataclasses, YAML/JSON/INI loaders. |
| [`nextssr.models`](models.md) | Dataclasses for SSR items, compound microsatellites, and provenance. |
| [`nextssr.compound`](compound.md) | Processor for compound microsatellite clustering. |
| [`nextssr.provenance`](provenance.md) | W3C JSON-LD and FAIR RO-Crate metadata manager. |
| [`nextssr.outputs`](exporters.md) | Exporters for GFF3 annotations and TSV tabular outputs. |

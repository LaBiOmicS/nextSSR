"""nextSSR - Next-generation High-Performance & FAIR-compliant Simple Sequence Repeat Platform.

Provides high-throughput identification of simple sequence repeats (SSRs / microsatellites),
in silico PCR primer design, Weber (1990) classification, and FAIR-compliant RO-Crate provenance tracking.

Modules:
    config: Configuration models and YAML/JSON/INI loaders.
    finder: SSR identification engine supporting parallel multi-core execution.
    primer: Oligonucleotide primer design engine with optional Primer3 backend.
    epcr: In silico e-PCR simulator for primer specificity testing.
    models: Dataclasses for SSRs, primer pairs, and execution provenance.
    compound: Processor for compound microsatellite clustering.
    artifacts: Output directory and run manifest manager.
    provenance: FAIR RO-Crate metadata generator.
"""

from nextssr.config import SSRConfig
from nextssr.finder import SSRFinder
from nextssr.primer import PrimerDesigner, PrimerPair
from nextssr.epcr import EPCRSimulator, AmpliconResult

__version__ = "0.1.1"

__all__ = [
    "SSRConfig",
    "SSRFinder",
    "PrimerDesigner",
    "PrimerPair",
    "EPCRSimulator",
    "AmpliconResult",
]

from dataclasses import dataclass, field
from typing import List, Optional
import platform
from nextssr.primer import PrimerPair


@dataclass
class SSRItem:
    """Represents an individual SSR detection with Weber (1990) and motif size classifications."""

    seq_id: str
    motif: str
    motif_length: int
    repeats: int
    start: int  # 1-based indexing
    end: int  # 1-based indexing
    sequence: str
    motif_class: str = ""  # e.g., 'dinucleotide', 'trinucleotide'
    weber_class: str = "Perfect"  # 'Perfect', 'Imperfect' (Weber 1990)
    flank_5p: str = ""
    flank_3p: str = ""
    primer_pair: Optional[PrimerPair] = None
    so_term: str = "SO:0000289"  # Sequence Ontology: microsatellite


@dataclass
class CompoundSSR:
    """Represents a compound SSR with Weber (1990) classification."""

    seq_id: str
    ssrs: List[SSRItem]
    start: int
    end: int
    compound_type: str
    full_pattern: str
    weber_class: str = (
        "Compound Perfect"  # 'Compound Perfect', 'Compound Disrupted' (Weber 1990)
    )
    so_term: str = "SO:0001061"  # Sequence Ontology: compound_microsatellite


@dataclass
class SequenceAnalysisResult:
    """Container for SSRs found in a single FASTA sequence."""

    seq_id: str
    seq_length: int
    ssrs: List[SSRItem] = field(default_factory=list)
    compounds: List[CompoundSSR] = field(default_factory=list)
    checksum: Optional[str] = None


@dataclass
class ExecutionProvenance:
    """FAIR-compliant execution provenance tracking for reproducibility."""

    tool_name: str = "nextSSR"
    tool_version: str = "0.1.0"
    python_version: str = platform.python_version()
    platform_info: str = platform.platform()
    execution_time: Optional[str] = None
    input_file_hash: Optional[str] = None
    parameters_hash: Optional[str] = None
    total_sequences: int = 0
    total_ssrs: int = 0
    total_compounds: int = 0
    total_primers_designed: int = 0
    device_used: str = "CPU"
    threads_used: int = 1

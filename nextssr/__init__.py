"""
nextSSR - Next-generation Simple Sequence Repeat Identification Tool.
"""

from nextssr.config import SSRConfig
from nextssr.finder import SSRFinder
from nextssr.primer import PrimerDesigner, PrimerPair
from nextssr.epcr import EPCRSimulator, AmpliconResult

__version__ = "0.1.0"

__all__ = [
    "SSRConfig",
    "SSRFinder",
    "PrimerDesigner",
    "PrimerPair",
    "EPCRSimulator",
    "AmpliconResult",
]

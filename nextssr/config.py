from dataclasses import dataclass, field
from typing import Dict, Optional
import os
import json
import yaml
import hashlib


@dataclass
class SSRConfig:
    """Configuration class for nextSSR identification and primer design.

    Attributes:
        unit_min_repeats (Dict[int, int]): Mapping of motif unit size to minimum repeat count threshold.
        max_compound_distance (int): Maximum base-pair distance between adjacent SSRs to be considered compound.
        threads (int): Number of parallel CPU worker threads.
        use_gpu (bool): Whether CUDA GPU acceleration is enabled.
        gpu_device_id (int): CUDA GPU device ordinal index.
        chunk_size_mb (int): Memory chunk size in megabytes for streaming inputs.
        batch_size (int): Number of FASTA records per processing batch.
        design_primers (bool): Flag indicating whether to perform PCR primer design.
        flank_len (int): Flanking region length (bp) extracted for primer design.
        opt_tm (float): Optimal melting temperature (°C) for designed primers.
        min_tm (float): Minimum allowable melting temperature (°C).
        max_tm (float): Maximum allowable melting temperature (°C).
        min_product_size (int): Minimum PCR amplicon product size (bp).
        max_product_size (int): Maximum PCR amplicon product size (bp).
        output_gff (bool): Enable GFF3 output generation.
        output_tsv (bool): Enable TSV tabular output generation.
        output_json_ld (bool): Enable JSON-LD metadata generation.
        seed (int): Random seed for reproducibility.
        generate_ro_crate (bool): Enable FAIR RO-Crate metadata generation.
    """

    # Motifs definition: unit_size -> min_repeats
    unit_min_repeats: Dict[int, int] = field(
        default_factory=lambda: {1: 10, 2: 6, 3: 5, 4: 5, 5: 5, 6: 5}
    )
    max_compound_distance: int = 100

    # Parallelism & Acceleration
    threads: int = os.cpu_count() or 4
    use_gpu: bool = False
    gpu_device_id: int = 0
    chunk_size_mb: int = 64
    batch_size: int = 1000

    # Primer Design
    design_primers: bool = True
    flank_len: int = 150
    opt_tm: float = 58.0
    min_tm: float = 50.0
    max_tm: float = 65.0
    min_product_size: int = 100
    max_product_size: int = 300

    # Output Formats & FAIR Standards
    output_gff: bool = True
    output_tsv: bool = True
    output_json_ld: bool = True

    # Reproducibility
    seed: int = 42
    generate_ro_crate: bool = True

    def get_hash(self) -> str:
        """Returns a deterministic SHA256 hash of configuration parameters for FAIR provenance.

        Returns:
            str: SHA256 hex digest representing the configuration signature.
        """
        config_str = json.dumps(
            {
                "unit_min_repeats": self.unit_min_repeats,
                "max_compound_distance": self.max_compound_distance,
                "seed": self.seed,
                "opt_tm": self.opt_tm,
            },
            sort_keys=True,
        )
        return hashlib.sha256(config_str.encode("utf-8")).hexdigest()

    @classmethod
    def generate_default_config(cls, filepath: str = "nextssr.yaml") -> str:
        """Generate a documented YAML configuration file for nextSSR.

        Args:
            filepath (str): Target output path for the generated YAML config file. Defaults to "nextssr.yaml".

        Returns:
            str: Absolute path of the created configuration file.
        """
        default_dict = {
            "nextssr": {
                "version": "0.1.1",
                "ssr_criteria": {
                    "unit_min_repeats": {
                        "1": 10,
                        "2": 6,
                        "3": 5,
                        "4": 5,
                        "5": 5,
                        "6": 5,
                    },
                    "max_compound_distance": 100,
                },
                "performance": {
                    "threads": os.cpu_count() or 4,
                    "use_gpu": False,
                    "batch_size": 1000,
                },
                "primer_design": {
                    "enabled": True,
                    "flank_length_bp": 150,
                    "optimal_tm_celsius": 58.0,
                    "min_tm_celsius": 50.0,
                    "max_tm_celsius": 65.0,
                    "min_product_size_bp": 100,
                    "max_product_size_bp": 300,
                },
                "fair_and_outputs": {
                    "output_gff3": True,
                    "output_tsv": True,
                    "generate_ro_crate": True,
                },
            }
        }
        with open(filepath, "w", encoding="utf-8") as f:
            f.write("# ==========================================\n")
            f.write("# nextSSR Configuration File\n")
            f.write("# Next-Generation SSR & Primer Design Platform\n")
            f.write("# ==========================================\n\n")
            yaml.dump(default_dict, f, default_flow_style=False, sort_keys=False)

        return os.path.abspath(filepath)

    @classmethod
    def from_file(
        cls,
        config_path: str,
        threads: Optional[int] = None,
        use_gpu: Optional[bool] = None,
    ) -> "SSRConfig":
        """Parse nextSSR YAML, JSON, or INI configuration file.

        Args:
            config_path (str): Path to the configuration file.
            threads (Optional[int]): Override for thread count.
            use_gpu (Optional[bool]): Override for GPU usage flag.

        Returns:
            SSRConfig: Populated SSRConfig instance.

        Raises:
            FileNotFoundError: If `config_path` does not exist.
        """
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Config file not found: {config_path}")

        # Check for YAML / JSON
        if config_path.endswith((".yaml", ".yml", ".json")):
            with open(config_path, "r", encoding="utf-8") as f:
                data = (
                    yaml.safe_load(f)
                    if config_path.endswith((".yaml", ".yml"))
                    else json.load(f)
                )

            n_cfg = data.get("nextssr", {})
            criteria = n_cfg.get("ssr_criteria", {})
            perf = n_cfg.get("performance", {})
            primer_cfg = n_cfg.get("primer_design", {})
            fair_cfg = n_cfg.get("fair_and_outputs", {})

            unit_reps = {
                int(k): int(v) for k, v in criteria.get("unit_min_repeats", {}).items()
            }

            return cls(
                unit_min_repeats=unit_reps or {1: 10, 2: 6, 3: 5, 4: 5, 5: 5, 6: 5},
                max_compound_distance=criteria.get("max_compound_distance", 100),
                threads=(
                    threads
                    if threads is not None
                    else perf.get("threads", os.cpu_count() or 4)
                ),
                use_gpu=use_gpu if use_gpu is not None else perf.get("use_gpu", False),
                batch_size=perf.get("batch_size", 1000),
                design_primers=primer_cfg.get("enabled", True),
                flank_len=primer_cfg.get("flank_length_bp", 150),
                opt_tm=primer_cfg.get("optimal_tm_celsius", 58.0),
                min_tm=primer_cfg.get("min_tm_celsius", 50.0),
                max_tm=primer_cfg.get("max_tm_celsius", 65.0),
                min_product_size=primer_cfg.get("min_product_size_bp", 100),
                max_product_size=primer_cfg.get("max_product_size_bp", 300),
                output_gff=fair_cfg.get("output_gff3", True),
                output_tsv=fair_cfg.get("output_tsv", True),
                generate_ro_crate=fair_cfg.get("generate_ro_crate", True),
            )

        # Fallback to INI parser
        unit_repeats = {}
        max_dist = 100
        gff_flag = True

        with open(config_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith("definition") or line.startswith("def"):
                    parts = line.split(":", 1)
                    if len(parts) > 1:
                        pairs = parts[1].strip().split()
                        for pair in pairs:
                            if "-" in pair:
                                unit, min_rep = map(int, pair.split("-"))
                                unit_repeats[unit] = min_rep
                elif line.startswith("interruptions") or line.startswith("int"):
                    parts = line.split(":", 1)
                    if len(parts) > 1:
                        max_dist = int(parts[1].strip())
                elif line.startswith("GFF"):
                    parts = line.split(":", 1)
                    if len(parts) > 1:
                        gff_flag = parts[1].strip().lower() == "true"

        return cls(
            unit_min_repeats=(
                unit_repeats if unit_repeats else {1: 10, 2: 6, 3: 5, 4: 5, 5: 5, 6: 5}
            ),
            max_compound_distance=max_dist,
            output_gff=gff_flag,
            threads=threads if threads is not None else 4,
            use_gpu=use_gpu if use_gpu is not None else False,
        )

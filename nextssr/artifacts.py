import os
import json
import datetime
from typing import Dict, Any
from nextssr.models import SequenceAnalysisResult, ExecutionProvenance
from nextssr.outputs.gff3 import GFF3Exporter
from nextssr.outputs.tsv import TSVExporter
from nextssr.provenance import FAIRProvenanceManager
from nextssr.config import SSRConfig

class ArtifactManager:
    """Manages structured output directories and artifacts for nextSSR execution runs."""

    def __init__(self, base_output_dir: str):
        self.base_dir = os.path.abspath(base_output_dir)
        self.annotations_dir = os.path.join(self.base_dir, "annotations")
        self.primers_dir = os.path.join(self.base_dir, "primers")
        self.provenance_dir = os.path.join(self.base_dir, "provenance")
        self.summary_dir = os.path.join(self.base_dir, "summary")

    def initialize_dirs(self):
        """Create directory structure for artifacts."""
        for d in [self.annotations_dir, self.primers_dir, self.provenance_dir, self.summary_dir]:
            os.makedirs(d, exist_ok=True)

    def save_artifacts(
        self,
        results: list,
        config: SSRConfig,
        provenance: ExecutionProvenance,
        input_filepath: str
    ) -> Dict[str, str]:
        """Save all run artifacts into structured directories."""
        self.initialize_dirs()

        # 1. Annotations
        gff_path = os.path.join(self.annotations_dir, "nextssr_results.gff3")
        GFF3Exporter.export(results, gff_path)

        # 2. Primers & TSV
        tsv_path = os.path.join(self.primers_dir, "nextssr_primers.tsv")
        TSVExporter.export(results, tsv_path)

        # 3. FAIR Provenance RO-Crate
        crate_path = FAIRProvenanceManager.generate_ro_crate(provenance, config, input_filepath, self.provenance_dir)

        # 4. Summary statistics
        summary_path = os.path.join(self.summary_dir, "nextssr_summary_statistics.txt")
        self._write_summary(summary_path, provenance)

        # 5. Global Run Manifest
        manifest_path = os.path.join(self.base_dir, "run_manifest.json")
        manifest = {
            "execution_id": provenance.parameters_hash or "run_01",
            "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "input_file": os.path.basename(input_filepath),
            "artifacts": {
                "gff3_annotation": os.path.relpath(gff_path, self.base_dir),
                "tsv_primers": os.path.relpath(tsv_path, self.base_dir),
                "ro_crate": os.path.relpath(crate_path, self.base_dir),
                "summary": os.path.relpath(summary_path, self.base_dir)
            },
            "metrics": {
                "total_sequences": provenance.total_sequences,
                "total_ssrs": provenance.total_ssrs,
                "total_compounds": provenance.total_compounds,
                "total_primers_designed": provenance.total_primers_designed,
                "execution_time": provenance.execution_time
            }
        }
        with open(manifest_path, 'w', encoding='utf-8') as f:
            json.dump(manifest, f, indent=2)

        return {
            "gff": gff_path,
            "tsv": tsv_path,
            "ro_crate": crate_path,
            "summary": summary_path,
            "manifest": manifest_path
        }

    def _write_summary(self, filepath: str, prov: ExecutionProvenance):
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write("===============================================\n")
            f.write("        nextSSR EXECUTION SUMMARY REPORT        \n")
            f.write("===============================================\n\n")
            f.write(f"Tool Version:              {prov.tool_name} {prov.tool_version}\n")
            f.write(f"Platform:                  {prov.platform_info}\n")
            f.write(f"Device Used:               {prov.device_used}\n")
            f.write(f"Threads Used:              {prov.threads_used}\n")
            f.write(f"Execution Time:            {prov.execution_time}\n")
            f.write(f"Input File SHA-256:        {prov.input_file_hash}\n\n")
            f.write("RESULTS OVERVIEW\n")
            f.write("-----------------------------------------------\n")
            f.write(f"Total Sequences Analyzed:  {prov.total_sequences}\n")
            f.write(f"Total SSRs Identified:     {prov.total_ssrs}\n")
            f.write(f"Total Compound SSRs:       {prov.total_compounds}\n")
            f.write(f"PCR Primers Designed:      {prov.total_primers_designed}\n")
            f.write("===============================================\n")

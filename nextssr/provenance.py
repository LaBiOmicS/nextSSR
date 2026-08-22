import json
import datetime
import hashlib
import os
from nextssr.models import ExecutionProvenance
from nextssr.config import SSRConfig


class FAIRProvenanceManager:
    """Generates FAIR-compliant RO-Crate and JSON-LD metadata for reproducibility and provenance."""

    @staticmethod
    def compute_file_sha256(filepath: str) -> str:
        """Compute SHA256 checksum of input file for data provenance."""
        sha256 = hashlib.sha256()
        with open(filepath, "rb") as f:
            while chunk := f.read(65536):
                sha256.update(chunk)
        return sha256.hexdigest()

    @classmethod
    def generate_ro_crate(
        cls,
        prov: ExecutionProvenance,
        config: SSRConfig,
        input_filepath: str,
        output_dir: str,
    ):
        """Export RO-Crate metadata (ro-crate-metadata.json) conforming to W3C and FAIR principles."""
        ro_crate_path = os.path.join(output_dir, "ro-crate-metadata.json")

        crate = {
            "@context": "https://w3id.org/ro/crate/1.1/context",
            "@graph": [
                {
                    "@type": "CreativeWork",
                    "@id": "ro-crate-metadata.json",
                    "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"},
                    "about": {"@id": "./"},
                },
                {
                    "@id": "./",
                    "@type": "Dataset",
                    "name": "nextSSR Microsatellite Analysis Results",
                    "description": "FAIR-compliant SSR and compound microsatellite identification dataset.",
                    "datePublished": datetime.datetime.now(
                        datetime.timezone.utc
                    ).isoformat(),
                    "license": "https://spdx.org/licenses/MIT.html",
                    "hasPart": [
                        {"@id": os.path.basename(input_filepath)},
                        {"@id": "nextssr_results.gff3"},
                        {"@id": "nextssr_results.tsv"},
                    ],
                },
                {
                    "@id": os.path.basename(input_filepath),
                    "@type": "File",
                    "name": os.path.basename(input_filepath),
                    "sha256": prov.input_file_hash,
                },
                {
                    "@id": "nextSSR_software",
                    "@type": "SoftwareApplication",
                    "name": prov.tool_name,
                    "softwareVersion": prov.tool_version,
                    "programmingLanguage": "Python " + prov.python_version,
                    "executionEnvironment": prov.platform_info,
                    "threadsUsed": prov.threads_used,
                    "deviceUsed": prov.device_used,
                },
                {
                    "@id": "execution_run",
                    "@type": "CreateAction",
                    "name": "SSR Identification Run",
                    "instrument": {"@id": "nextSSR_software"},
                    "object": {"@id": os.path.basename(input_filepath)},
                    "endTime": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    "configHash": config.get_hash(),
                },
            ],
        }

        with open(ro_crate_path, "w", encoding="utf-8") as f:
            json.dump(crate, f, indent=2)

        return ro_crate_path

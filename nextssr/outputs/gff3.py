import datetime
from typing import List
from nextssr.models import SequenceAnalysisResult


class GFF3Exporter:
    """FAIR Interoperable GFF3 format exporter using Sequence Ontology (SO) terms and Weber classifications."""

    @staticmethod
    def export(
        results: List[SequenceAnalysisResult],
        output_path: str,
        tool_version: str = "0.1.0",
    ):
        with open(output_path, "w", encoding="utf-8") as f:
            f.write("##gff-version 3\n")
            f.write(f"#!Date {datetime.date.today().isoformat()}\n")
            f.write(f"#!Source-version nextSSR {tool_version}\n")
            f.write(
                "#!Sequence-Ontology-Terms SO:0000289(microsatellite), SO:0001061(compound_microsatellite)\n"
            )

            for res in results:
                f.write(f"##sequence-region {res.seq_id} 1 {res.seq_length}\n")

                # Single SSRs
                idx = 1
                for ssr in res.ssrs:
                    feature_type = "microsatellite"
                    note = (
                        f"motif_class={ssr.motif_class};weber_class={ssr.weber_class};"
                        f"motif={ssr.motif};repeats={ssr.repeats};sequence={ssr.sequence}"
                    )
                    attr = f"ID={res.seq_id}.ssr{idx};Ontology_term={ssr.so_term};Note={note}"
                    f.write(
                        f"{res.seq_id}\tnextSSR\t{feature_type}\t{ssr.start}\t{ssr.end}\t.\t+\t.\t{attr}\n"
                    )
                    idx += 1

                # Compound SSRs
                c_idx = 1
                for comp in res.compounds:
                    feature_type = "compound_microsatellite"
                    note = f"weber_class={comp.weber_class};pattern={comp.full_pattern}"
                    comp_attr = f"ID={res.seq_id}.compound{c_idx};Ontology_term={comp.so_term};Note={note}"
                    f.write(
                        f"{res.seq_id}\tnextSSR\t{feature_type}\t{comp.start}\t{comp.end}\t.\t+\t.\t{comp_attr}\n"
                    )
                    c_idx += 1

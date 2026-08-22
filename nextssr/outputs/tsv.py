import os
from typing import List
from nextssr.models import SequenceAnalysisResult

class TSVExporter:
    """Tabular TSV exporter including Motif Class, Weber (1990) Class, and Primers."""
    
    @staticmethod
    def export(results: List[SequenceAnalysisResult], output_path: str):
        with open(output_path, 'w', encoding='utf-8') as f:
            headers = [
                "Seq_ID", "SSR_Nr", "Motif_Class", "Weber_Classification", "Motif", "Repeats", "Size_bp", "Start", "End",
                "Forward_Primer", "Forward_Tm", "Forward_GC",
                "Reverse_Primer", "Reverse_Tm", "Reverse_GC",
                "Product_Size_bp", "Primer_Status",
                "Flank_5p", "Flank_3p", "Sequence"
            ]
            f.write("\t".join(headers) + "\n")

            for res in results:
                idx = 1
                for ssr in res.ssrs:
                    p = ssr.primer_pair
                    f_seq = p.forward_seq if p else ""
                    f_tm = str(p.forward_tm) if p else ""
                    f_gc = str(p.forward_gc) if p else ""
                    r_seq = p.reverse_seq if p else ""
                    r_tm = str(p.reverse_tm) if p else ""
                    r_gc = str(p.reverse_gc) if p else ""
                    p_size = str(p.product_size) if p else ""
                    p_status = p.status if p else "N/A"

                    row = [
                        res.seq_id,
                        f"SSR_{idx}",
                        ssr.motif_class,
                        ssr.weber_class,
                        ssr.motif,
                        str(ssr.repeats),
                        str(len(ssr.sequence)),
                        str(ssr.start),
                        str(ssr.end),
                        f_seq, f_tm, f_gc,
                        r_seq, r_tm, r_gc,
                        p_size, p_status,
                        ssr.flank_5p,
                        ssr.flank_3p,
                        ssr.sequence
                    ]
                    f.write("\t".join(row) + "\n")
                    idx += 1

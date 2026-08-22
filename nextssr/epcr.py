from dataclasses import dataclass
from typing import List, Tuple, Dict
from nextssr.utils import parse_fasta


@dataclass
class AmpliconResult:
    """Dataclass holding in silico e-PCR amplification results."""

    seq_id: str
    forward_primer: str
    reverse_primer: str
    forward_start: int
    forward_end: int
    forward_mismatches: int
    reverse_start: int
    reverse_end: int
    reverse_mismatches: int
    product_size: int
    amplicon_sequence: str
    status: str = "SPECIFIC"  # SPECIFIC, OFF_TARGET, MISMATCH


class EPCRSimulator:
    """High-performance in silico e-PCR Simulator for validating PCR primers on target FASTA genomes."""

    def __init__(
        self,
        max_mismatches: int = 2,
        min_product_size: int = 50,
        max_product_size: int = 1000,
    ):
        self.max_mismatches = max_mismatches
        self.min_product_size = min_product_size
        self.max_product_size = max_product_size

    @staticmethod
    def reverse_complement(seq: str) -> str:
        """Return the reverse complement of a DNA sequence."""
        complement = str.maketrans("ATCGNatcgn", "TAGCNtagcn")
        return seq.translate(complement)[::-1]

    def _find_matches(self, primer: str, sequence: str) -> List[Tuple[int, int, int]]:
        """Find primer binding sites in a sequence allowing up to max_mismatches.

        Returns list of (start_idx, end_idx, mismatch_count).
        """
        primer_len = len(primer)
        seq_len = len(sequence)
        matches = []

        primer_upper = primer.upper()
        seq_upper = sequence.upper()

        for i in range(seq_len - primer_len + 1):
            subseq = seq_upper[i : i + primer_len]
            mismatches = sum(1 for a, b in zip(primer_upper, subseq) if a != b)
            if mismatches <= self.max_mismatches:
                matches.append((i + 1, i + primer_len, mismatches))

        return matches

    def simulate_sequence(
        self, seq_id: str, sequence: str, forward_primer: str, reverse_primer: str
    ) -> List[AmpliconResult]:
        """Simulate in silico PCR for a single target sequence."""
        amplicons = []

        # Forward primer binds on positive strand (+)
        fwd_hits = self._find_matches(forward_primer, sequence)

        # Reverse primer binds on negative strand (-), so match reverse complement on positive strand
        rev_comp_primer = self.reverse_complement(reverse_primer)
        rev_hits = self._find_matches(rev_comp_primer, sequence)

        for f_start, f_end, f_mismatches in fwd_hits:
            for r_start, r_end, r_mismatches in rev_hits:
                # Forward primer must be upstream of Reverse primer
                if r_end > f_start:
                    prod_size = r_end - f_start + 1
                    if self.min_product_size <= prod_size <= self.max_product_size:
                        amplicon_seq = sequence[f_start - 1 : r_end]
                        status = "SPECIFIC"
                        if f_mismatches > 0 or r_mismatches > 0:
                            status = "MISMATCH"

                        amplicons.append(
                            AmpliconResult(
                                seq_id=seq_id,
                                forward_primer=forward_primer,
                                reverse_primer=reverse_primer,
                                forward_start=f_start,
                                forward_end=f_end,
                                forward_mismatches=f_mismatches,
                                reverse_start=r_start,
                                reverse_end=r_end,
                                reverse_mismatches=r_mismatches,
                                product_size=prod_size,
                                amplicon_sequence=amplicon_seq,
                                status=status,
                            )
                        )

        # Mark off-target if multiple amplicons found for same primer pair
        if len(amplicons) > 1:
            for amp in amplicons:
                amp.status = "OFF_TARGET"

        return amplicons

    def run_fasta(
        self, fasta_path: str, forward_primer: str, reverse_primer: str
    ) -> List[AmpliconResult]:
        """Run in silico e-PCR against a FASTA file for a specific primer pair."""
        results = []
        for seq_id, sequence in parse_fasta(fasta_path):
            amps = self.simulate_sequence(
                seq_id, sequence, forward_primer, reverse_primer
            )
            results.extend(amps)
        return results

    def run_primers_tsv(
        self, fasta_path: str, primers_tsv_path: str
    ) -> Dict[str, List[AmpliconResult]]:
        """Run in silico e-PCR for all primer pairs in a nextSSR primers TSV file."""
        all_results = {}

        # Parse TSV
        primer_pairs = []
        with open(primers_tsv_path, "r", encoding="utf-8") as f:
            headers = f.readline().strip().split("\t")
            if "Forward_Primer" in headers and "Reverse_Primer" in headers:
                fwd_idx = headers.index("Forward_Primer")
                rev_idx = headers.index("Reverse_Primer")
                ssr_id_idx = headers.index("SSR_ID") if "SSR_ID" in headers else 0

                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) > max(fwd_idx, rev_idx):
                        fwd = parts[fwd_idx]
                        rev = parts[rev_idx]
                        ssr_id = parts[ssr_id_idx]
                        if fwd and rev and fwd != "N/A" and rev != "N/A":
                            primer_pairs.append((ssr_id, fwd, rev))

        # Run e-PCR for each pair across FASTA
        for ssr_id, fwd, rev in primer_pairs:
            pair_key = f"{ssr_id} ({fwd} / {rev})"
            amps = self.run_fasta(fasta_path, fwd, rev)
            all_results[pair_key] = amps

        return all_results

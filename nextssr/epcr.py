from dataclasses import dataclass
from typing import List, Tuple, Dict, Set

IUPAC_MAP: Dict[str, Set[str]] = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


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
    amplicon_gc: float = 0.0
    status: str = "SPECIFIC"  # SPECIFIC, OFF_TARGET, MISMATCH, 3PRIME_MISMATCH


class EPCRSimulator:
    """High-performance in silico e-PCR Simulator for validating PCR primers on target FASTA genomes."""

    def __init__(
        self,
        max_mismatches: int = 2,
        max_3prime_mismatches: int = 0,
        three_prime_len: int = 5,
        min_product_size: int = 50,
        max_product_size: int = 1000,
    ):
        self.max_mismatches = max_mismatches
        self.max_3prime_mismatches = max_3prime_mismatches
        self.three_prime_len = three_prime_len
        self.min_product_size = min_product_size
        self.max_product_size = max_product_size

    @staticmethod
    def reverse_complement(seq: str) -> str:
        """Return the reverse complement of a DNA sequence including IUPAC bases."""
        comp_map = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "U": "A",
            "R": "Y",
            "Y": "R",
            "S": "S",
            "W": "W",
            "K": "M",
            "M": "K",
            "B": "V",
            "V": "B",
            "D": "H",
            "H": "D",
            "N": "N",
        }
        seq_upper = seq.upper()
        res = [comp_map.get(b, b) for b in reversed(seq_upper)]
        return "".join(res)

    @staticmethod
    def calculate_gc(sequence: str) -> float:
        """Calculate GC percentage of sequence."""
        if not sequence:
            return 0.0
        seq_u = sequence.upper()
        gc_count = seq_u.count("G") + seq_u.count("C")
        return round((gc_count / len(seq_u)) * 100, 2)

    def _match_bases(self, primer_base: str, target_base: str) -> bool:
        """Check if primer base matches target base under IUPAC rules."""
        p_u = primer_base.upper()
        t_u = target_base.upper()
        if t_u == "N" and p_u != "N":
            return False
        allowed = IUPAC_MAP.get(p_u, {p_u})
        target_allowed = IUPAC_MAP.get(t_u, {t_u})
        return bool(allowed.intersection(target_allowed))

    def _find_matches(
        self, primer: str, sequence: str
    ) -> List[Tuple[int, int, int, int]]:
        """Find primer binding sites in a sequence allowing up to max_mismatches.

        Returns list of (start_idx, end_idx, total_mismatches, 3prime_mismatches).
        """
        primer_len = len(primer)
        seq_len = len(sequence)
        matches = []

        primer_upper = primer.upper()
        seq_upper = sequence.upper()

        three_prime_start_idx = max(0, primer_len - self.three_prime_len)

        for i in range(seq_len - primer_len + 1):
            subseq = seq_upper[i : i + primer_len]
            total_mismatches = 0
            three_prime_mismatches = 0

            for idx, (p_base, s_base) in enumerate(zip(primer_upper, subseq)):
                if not self._match_bases(p_base, s_base):
                    total_mismatches += 1
                    if idx >= three_prime_start_idx:
                        three_prime_mismatches += 1

            if total_mismatches <= self.max_mismatches:
                matches.append(
                    (i + 1, i + primer_len, total_mismatches, three_prime_mismatches)
                )

        return matches

    def simulate_sequence(
        self, seq_id: str, sequence: str, forward_primer: str, reverse_primer: str
    ) -> List[AmpliconResult]:
        """Simulate in silico PCR for a single target sequence.

        Args:
            seq_id (str): Identifier header of target sequence.
            sequence (str): Nucleotide sequence template.
            forward_primer (str): 5' forward primer sequence.
            reverse_primer (str): 3' reverse primer sequence.

        Returns:
            List[AmpliconResult]: List of predicted amplicon result dataclasses.
        """
        amplicons = []

        # Forward primer binds on positive strand (+)
        fwd_hits = self._find_matches(forward_primer, sequence)

        # Reverse primer binds on negative strand (-), so match reverse complement on positive strand
        rev_comp_primer = self.reverse_complement(reverse_primer)
        rev_hits = self._find_matches(rev_comp_primer, sequence)

        for f_start, f_end, f_mismatches, f_3p_mismatches in fwd_hits:
            for r_start, r_end, r_mismatches, r_3p_mismatches in rev_hits:
                # Forward primer must be upstream of Reverse primer
                if r_end > f_start:
                    prod_size = r_end - f_start + 1
                    if self.min_product_size <= prod_size <= self.max_product_size:
                        amplicon_seq = sequence[f_start - 1 : r_end]
                        amplicon_gc = self.calculate_gc(amplicon_seq)

                        status = "SPECIFIC"
                        if (
                            f_3p_mismatches > self.max_3prime_mismatches
                            or r_3p_mismatches > self.max_3prime_mismatches
                        ):
                            status = "3PRIME_MISMATCH"
                        elif f_mismatches > 0 or r_mismatches > 0:
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
                                amplicon_gc=amplicon_gc,
                                status=status,
                            )
                        )

        # Mark off-target if multiple amplicons found for same primer pair
        if len(amplicons) > 1:
            for amp in amplicons:
                if amp.status != "3PRIME_MISMATCH":
                    amp.status = "OFF_TARGET"

        return amplicons

    def run_fasta(
        self, fasta_path: str, forward_primer: str, reverse_primer: str
    ) -> List[AmpliconResult]:
        """Run in silico e-PCR against a FASTA file for a specific primer pair.

        Args:
            fasta_path (str): Path to FASTA genome file.
            forward_primer (str): Forward primer sequence.
            reverse_primer (str): Reverse primer sequence.

        Returns:
            List[AmpliconResult]: List of all amplicons found across sequences in FASTA.
        """
        from nextssr.utils import parse_fasta

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
        """Run in silico e-PCR for all primer pairs listed in a nextSSR primers TSV file.

        Args:
            fasta_path (str): Path to FASTA genome file.
            primers_tsv_path (str): Path to nextSSR primers TSV file.

        Returns:
            Dict[str, List[AmpliconResult]]: Mapping of primer pair key to amplicon results.
        """
        all_results = {}

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

        for ssr_id, fwd, rev in primer_pairs:
            pair_key = f"{ssr_id} ({fwd} / {rev})"
            amps = self.run_fasta(fasta_path, fwd, rev)
            all_results[pair_key] = amps

        return all_results

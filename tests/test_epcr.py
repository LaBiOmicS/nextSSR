import os
import tempfile
from nextssr.epcr import EPCRSimulator


def test_epcr_simulator_basic():
    # Sequence containing Forward primer GATTACA at 10..16 and Reverse primer ACGTACGT at 150..157
    forward_primer = "GATTACAAGCTACG"
    reverse_primer = "ACGTACGTACGT"
    rev_comp_reverse = EPCRSimulator.reverse_complement(reverse_primer)

    target_seq = (
        "NNNNNNNNNN" + forward_primer + ("A" * 100) + rev_comp_reverse + "NNNNNNNNNN"
    )

    simulator = EPCRSimulator(
        max_mismatches=1, min_product_size=50, max_product_size=300
    )
    amplicons = simulator.simulate_sequence(
        "seq1", target_seq, forward_primer, reverse_primer
    )

    assert len(amplicons) == 1
    amp = amplicons[0]
    assert amp.seq_id == "seq1"
    assert amp.forward_mismatches == 0
    assert amp.reverse_mismatches == 0
    assert amp.product_size == len(forward_primer) + 100 + len(reverse_primer)
    assert amp.status in ["SPECIFIC", "MISMATCH"]


def test_epcr_fasta_file():
    forward_primer = "ATCGATCGAT"
    reverse_primer = "CGATCGATCG"
    rev_comp_reverse = EPCRSimulator.reverse_complement(reverse_primer)

    fasta_content = (
        ">target_chr1 Test Chromosome 1\n"
        "GCGCGC" + forward_primer + ("T" * 80) + rev_comp_reverse + "GCGCGC\n"
    )

    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as f:
        f.write(fasta_content)
        fasta_path = f.name

    try:
        simulator = EPCRSimulator(max_mismatches=0)
        amplicons = simulator.run_fasta(fasta_path, forward_primer, reverse_primer)
        assert len(amplicons) == 1
        assert amplicons[0].seq_id == "target_chr1"
        assert amplicons[0].product_size == 10 + 80 + 10
    finally:
        if os.path.exists(fasta_path):
            os.remove(fasta_path)

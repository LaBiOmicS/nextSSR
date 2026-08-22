from nextssr.epcr import EPCRSimulator


def test_epcr_simulator_basic():
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


def test_epcr_iupac_degenerate_matching():
    # Degenerate primer R = A or G, Y = C or T
    forward_primer = "GATTRCAAGCTACG"  # R matches A or G
    reverse_primer = "ACGTAYGTACGT"  # Y matches C or T
    rev_comp_reverse = EPCRSimulator.reverse_complement(reverse_primer)

    target_seq = "GCGCGC" + "GATTACAAGCTACG" + ("C" * 60) + rev_comp_reverse + "GCGCGC"

    simulator = EPCRSimulator(max_mismatches=0)
    amplicons = simulator.simulate_sequence(
        "seq_deg", target_seq, forward_primer, reverse_primer
    )

    assert len(amplicons) == 1
    assert amplicons[0].status == "SPECIFIC"
    assert amplicons[0].amplicon_gc > 0.0


def test_epcr_3prime_anchor_mismatch():
    forward_primer = "GATTACAAGCTACG"  # 3' end is TACG
    reverse_primer = "ACGTACGTACGT"
    rev_comp_reverse = EPCRSimulator.reverse_complement(reverse_primer)

    # Introduce a 1-base mismatch at the 3' end of forward primer
    mismatched_forward_target = "GATTACAAGCTACT"

    target_seq = (
        "GCGCGC" + mismatched_forward_target + ("A" * 60) + rev_comp_reverse + "GCGCGC"
    )

    # max_3prime_mismatches=0 will flag 3PRIME_MISMATCH
    simulator = EPCRSimulator(max_mismatches=1, max_3prime_mismatches=0)
    amplicons = simulator.simulate_sequence(
        "seq_3p", target_seq, forward_primer, reverse_primer
    )

    assert len(amplicons) == 1
    assert amplicons[0].status == "3PRIME_MISMATCH"

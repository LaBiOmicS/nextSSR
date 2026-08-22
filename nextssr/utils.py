from typing import Generator, Tuple, List


def parse_fasta(filepath: str) -> Generator[Tuple[str, str], None, None]:
    """Simple FASTA generator yielding (header_id, sequence)."""
    current_id = None
    current_seq: List[str] = []

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    yield current_id, "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            yield current_id, "".join(current_seq)

import gzip
from typing import Generator, Tuple, List


def parse_fasta(filepath: str) -> Generator[Tuple[str, str], None, None]:
    """Memory-efficient FASTA streaming generator yielding (header_id, sequence).
    
    Supports both uncompressed (.fasta, .fa, .fna) and gzip-compressed (.gz, .bgz) files seamlessly.
    """
    current_id = None
    current_seq: List[str] = []

    is_gz = filepath.endswith(".gz") or filepath.endswith(".bgz") or filepath.endswith(".gzip")
    open_fn = gzip.open if is_gz else open

    with open_fn(filepath, "rt", encoding="utf-8", errors="replace") as f:
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

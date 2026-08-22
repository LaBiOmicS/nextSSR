import re
from typing import List, Tuple, Generator, Dict
from concurrent.futures import ProcessPoolExecutor
from nextssr.models import SSRItem, SequenceAnalysisResult
from nextssr.config import SSRConfig
from nextssr.gpu import GPUAccelerator

MOTIF_CLASS_MAP = {
    1: "mononucleotide",
    2: "dinucleotide",
    3: "trinucleotide",
    4: "tetranucleotide",
    5: "pentanucleotide",
    6: "hexanucleotide",
}

_PERFECT_PATTERNS: Dict[Tuple[int, int], Tuple[re.Pattern, List[re.Pattern]]] = {}


def _get_perfect_patterns(
    motif_len: int, min_reps: int
) -> Tuple[re.Pattern, List[re.Pattern]]:
    key = (motif_len, min_reps)
    if key not in _PERFECT_PATTERNS:
        p = re.compile(f"(([ACGT]{{{motif_len}}})\\2{{{min_reps - 1},}})")
        sub_p_list = []
        for sub_len in range(1, motif_len):
            if motif_len % sub_len == 0:
                sub_p_list.append(
                    re.compile(
                        f"^([ACGT]{{{sub_len}}})\\1{{{motif_len // sub_len - 1}}}$"
                    )
                )
        _PERFECT_PATTERNS[key] = (p, sub_p_list)
    return _PERFECT_PATTERNS[key]


def _analyze_single_sequence_worker(
    args: Tuple[str, str, dict, int],
) -> SequenceAnalysisResult:
    """High-speed single sequence worker optimized for single-pass scanning."""
    seq_id, sequence, unit_min_repeats, flank_len = args
    cleaned_seq = re.sub(r"[\d\s>]", "", sequence).upper()
    seq_len = len(cleaned_seq)
    ssrs: List[SSRItem] = []

    for motif_len, min_reps in unit_min_repeats.items():
        motif_class = MOTIF_CLASS_MAP.get(motif_len, f"{motif_len}-mer")
        pattern, sub_patterns = _get_perfect_patterns(motif_len, min_reps)

        for match in pattern.finditer(cleaned_seq):
            full_match = match.group(1)
            motif = match.group(2)

            if sub_patterns and any(sp.match(motif) for sp in sub_patterns):
                continue

            repeats = len(full_match) // motif_len
            end_pos = match.end()
            start_pos = end_pos - len(full_match) + 1

            flank_5p_start = max(0, start_pos - 1 - flank_len)
            flank_5p = cleaned_seq[flank_5p_start : start_pos - 1]

            flank_3p_end = min(seq_len, end_pos + flank_len)
            flank_3p = cleaned_seq[end_pos:flank_3p_end]

            ssrs.append(
                SSRItem(
                    seq_id=seq_id,
                    motif=motif,
                    motif_length=motif_len,
                    repeats=repeats,
                    start=start_pos,
                    end=end_pos,
                    sequence=full_match,
                    motif_class=motif_class,
                    weber_class="Perfect",
                    flank_5p=flank_5p,
                    flank_3p=flank_3p,
                )
            )

    ssrs.sort(key=lambda x: x.start)
    return SequenceAnalysisResult(seq_id=seq_id, seq_length=seq_len, ssrs=ssrs)


class SSRFinder:
    """Core SSR identification engine supporting Weber (1990) classification and Motif Size mapping."""

    def __init__(self, config: SSRConfig, flank_len: int = 150):
        self.config = config
        self.flank_len = flank_len
        self.gpu_acc = GPUAccelerator(config.gpu_device_id) if config.use_gpu else None

    def analyze_sequence(self, seq_id: str, sequence: str) -> SequenceAnalysisResult:
        """Single sequence analyzer."""
        return _analyze_single_sequence_worker(
            (seq_id, sequence, self.config.unit_min_repeats, self.flank_len)
        )

    def analyze_batch_parallel(
        self, sequence_generator: Generator[Tuple[str, str], None, None]
    ) -> Generator[SequenceAnalysisResult, None, None]:
        """Stream sequences and process in parallel using multi-processing pool."""
        # Buffer first items to check total sequence count
        first_batch = []
        for item in sequence_generator:
            first_batch.append(item)
            if len(first_batch) >= self.config.batch_size:
                break

        # If sequence count is small or single sequence, execute inline to avoid process spawn overhead
        if len(first_batch) == 1 and not hasattr(sequence_generator, "__next__"):
            seq_id, seq_str = first_batch[0]
            yield _analyze_single_sequence_worker(
                (seq_id, seq_str, self.config.unit_min_repeats, self.flank_len)
            )
            return

        if self.config.threads == 1:
            for seq_id, seq_str in first_batch:
                yield _analyze_single_sequence_worker(
                    (seq_id, seq_str, self.config.unit_min_repeats, self.flank_len)
                )
            for seq_id, seq_str in sequence_generator:
                yield _analyze_single_sequence_worker(
                    (seq_id, seq_str, self.config.unit_min_repeats, self.flank_len)
                )
        else:
            with ProcessPoolExecutor(max_workers=self.config.threads) as executor:
                batch = [
                    (seq_id, seq_str, self.config.unit_min_repeats, self.flank_len)
                    for seq_id, seq_str in first_batch
                ]
                results = executor.map(
                    _analyze_single_sequence_worker,
                    batch,
                    chunksize=max(1, len(batch) // self.config.threads),
                )
                for res in results:
                    yield res

                batch.clear()
                for seq_id, seq_str in sequence_generator:
                    batch.append(
                        (seq_id, seq_str, self.config.unit_min_repeats, self.flank_len)
                    )
                    if len(batch) >= self.config.batch_size:
                        results = executor.map(
                            _analyze_single_sequence_worker,
                            batch,
                            chunksize=max(1, len(batch) // self.config.threads),
                        )
                        for res in results:
                            yield res
                        batch.clear()

                if batch:
                    results = executor.map(
                        _analyze_single_sequence_worker, batch, chunksize=1
                    )
                    for res in results:
                        yield res

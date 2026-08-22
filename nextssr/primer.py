from dataclasses import dataclass
from typing import Optional, List, Dict, Tuple
import logging
import math

logger = logging.getLogger("nextssr.primer")

@dataclass
class PrimerPair:
    """Dataclass holding designed primer pairs and quality parameters."""
    forward_seq: str
    forward_tm: float
    forward_gc: float
    forward_start: int
    forward_length: int
    
    reverse_seq: str
    reverse_tm: float
    reverse_gc: float
    reverse_start: int
    reverse_length: int

    product_size: int
    pair_penalty: float = 0.0
    status: str = "OK"
    error_reason: Optional[str] = None

class PrimerDesigner:
    """High-reliability Primer Designer with optional Primer3 C-engine and fallback thermodynamics."""

    def __init__(
        self,
        opt_tm: float = 58.0,
        min_tm: float = 50.0,
        max_tm: float = 65.0,
        opt_size: int = 20,
        min_size: int = 18,
        max_size: int = 25,
        min_gc: float = 35.0,
        max_gc: float = 65.0,
        min_product_size: int = 70,
        max_product_size: int = 300
    ):
        self.opt_tm = opt_tm
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.opt_size = opt_size
        self.min_size = min_size
        self.max_size = max_size
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.min_product_size = min_product_size
        self.max_product_size = max_product_size

        # Check for primer3-py
        self.use_primer3_lib = False
        try:
            import primer3
            self.primer3 = primer3
            self.use_primer3_lib = True
            logger.info("Primer3 C-engine backend detected and activated.")
        except ImportError:
            logger.info("primer3-py not installed. Using native nextSSR thermodynamic primer engine.")

    @staticmethod
    def calculate_gc(seq: str) -> float:
        """Calculate GC percentage of sequence."""
        if not seq:
            return 0.0
        gc_count = sum(1 for base in seq.upper() if base in ('G', 'C'))
        return round((gc_count / len(seq)) * 100.0, 2)

    @staticmethod
    def calculate_tm(seq: str) -> float:
        """Calculate Tm using nearest-neighbor thermodynamic formula for oligonucleotides."""
        seq = seq.upper()
        if len(seq) < 14:
            w_tm = (sum(1 for b in seq if b in ('A', 'T')) * 2) + (sum(1 for b in seq if b in ('G', 'C')) * 4)
            return float(w_tm)

        nn_values = {
            'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2),
            'AT': (-7.2, -20.4), 'TA': (-7.2, -21.3),
            'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
            'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
            'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0),
            'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
            'CG': (-10.6, -27.2), 'GC': (-9.8, -24.4),
            'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9)
        }
        
        delta_h = 0.0
        delta_s = -10.8  # Initiation entropy
        
        for i in range(len(seq) - 1):
            pair = seq[i:i+2]
            if pair in nn_values:
                h, s = nn_values[pair]
                delta_h += h
                delta_s += s

        R = 1.987
        conc = 0.00000025  # 250 nM primer concentration
        tm_kelvin = (delta_h * 1000) / (delta_s + R * math.log(conc / 4))
        tm_celsius = tm_kelvin - 273.15 + 16.6 * math.log10(0.05)
        return round(tm_celsius, 2)

    @staticmethod
    def reverse_complement(seq: str) -> str:
        """Return reverse complement of sequence."""
        trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        return seq.translate(trans)[::-1]

    def design_primers(
        self,
        flank_5p: str,
        target_seq: str,
        flank_3p: str
    ) -> PrimerPair:
        """Design forward and reverse primers spanning across the target SSR."""
        full_seq = (flank_5p + target_seq + flank_3p).upper()
        target_start_idx = len(flank_5p)
        target_end_idx = target_start_idx + len(target_seq)

        if 'N' in full_seq:
            return PrimerPair(
                forward_seq="", forward_tm=0.0, forward_gc=0.0, forward_start=0, forward_length=0,
                reverse_seq="", reverse_tm=0.0, reverse_gc=0.0, reverse_start=0, reverse_length=0,
                product_size=0, status="REJECTED", error_reason="Contains ambiguous N bases"
            )

        if len(flank_5p) < self.min_size or len(flank_3p) < self.min_size:
            return PrimerPair(
                forward_seq="", forward_tm=0.0, forward_gc=0.0, forward_start=0, forward_length=0,
                reverse_seq="", reverse_tm=0.0, reverse_gc=0.0, reverse_start=0, reverse_length=0,
                product_size=0, status="REJECTED", error_reason="Flanking regions too short"
            )

        if self.use_primer3_lib:
            try:
                res = self.primer3.bindings.design_primers(
                    seq_args={
                        'SEQUENCE_TEMPLATE': full_seq,
                        'SEQUENCE_TARGET': [target_start_idx, len(target_seq)]
                    },
                    global_args={
                        'PRIMER_OPT_SIZE': self.opt_size,
                        'PRIMER_MIN_SIZE': self.min_size,
                        'PRIMER_MAX_SIZE': self.max_size,
                        'PRIMER_OPT_TM': self.opt_tm,
                        'PRIMER_MIN_TM': self.min_tm,
                        'PRIMER_MAX_TM': self.max_tm,
                        'PRIMER_MIN_GC': self.min_gc,
                        'PRIMER_MAX_GC': self.max_gc,
                        'PRIMER_PRODUCT_SIZE_RANGE': [[self.min_product_size, self.max_product_size]]
                    }
                )
                if res.get('PRIMER_PAIR_NUM_RETURNED', 0) > 0:
                    f_seq = res['PRIMER_LEFT_0_SEQUENCE']
                    r_seq = res['PRIMER_RIGHT_0_SEQUENCE']
                    f_start, f_len = res['PRIMER_LEFT_0']
                    r_start, r_len = res['PRIMER_RIGHT_0']
                    prod_size = res['PRIMER_PAIR_0_PRODUCT_SIZE']
                    
                    return PrimerPair(
                        forward_seq=f_seq,
                        forward_tm=round(res['PRIMER_LEFT_0_TM'], 2),
                        forward_gc=round(res['PRIMER_LEFT_0_GC_PERCENT'], 2),
                        forward_start=f_start + 1,
                        forward_length=f_len,
                        reverse_seq=r_seq,
                        reverse_tm=round(res['PRIMER_RIGHT_0_TM'], 2),
                        reverse_gc=round(res['PRIMER_RIGHT_0_GC_PERCENT'], 2),
                        reverse_start=r_start + 1,
                        reverse_length=r_len,
                        product_size=prod_size,
                        status="OK"
                    )
            except Exception as e:
                logger.debug(f"primer3-py fallback triggered: {e}")

        candidates_f = []
        candidates_r = []

        for length in range(self.min_size, self.max_size + 1):
            for i in range(0, target_start_idx - length + 1):
                f_seq = full_seq[i : i + length]
                tm = self.calculate_tm(f_seq)
                gc = self.calculate_gc(f_seq)
                if self.min_tm <= tm <= self.max_tm and self.min_gc <= gc <= self.max_gc:
                    penalty = abs(tm - self.opt_tm) + abs(length - self.opt_size) * 0.5
                    candidates_f.append((i, length, f_seq, tm, gc, penalty))

        for length in range(self.min_size, self.max_size + 1):
            for j in range(target_end_idx, len(full_seq) - length + 1):
                r_template = full_seq[j : j + length]
                r_seq = self.reverse_complement(r_template)
                tm = self.calculate_tm(r_seq)
                gc = self.calculate_gc(r_seq)
                if self.min_tm <= tm <= self.max_tm and self.min_gc <= gc <= self.max_gc:
                    penalty = abs(tm - self.opt_tm) + abs(length - self.opt_size) * 0.5
                    candidates_r.append((j, length, r_seq, tm, gc, penalty))

        if not candidates_f or not candidates_r:
            return PrimerPair(
                forward_seq="", forward_tm=0.0, forward_gc=0.0, forward_start=0, forward_length=0,
                reverse_seq="", reverse_tm=0.0, reverse_gc=0.0, reverse_start=0, reverse_length=0,
                product_size=0, status="REJECTED", error_reason="No primers met Tm/GC constraints"
            )

        best_pair = None
        min_total_penalty = float('inf')

        for f in candidates_f:
            for r in candidates_r:
                f_start, f_len, f_seq, f_tm, f_gc, f_pen = f
                r_start, r_len, r_seq, r_tm, r_gc, r_pen = r
                
                prod_size = (r_start + r_len) - f_start
                if self.min_product_size <= prod_size <= self.max_product_size:
                    tm_diff = abs(f_tm - r_tm)
                    if tm_diff <= 4.0:
                        total_penalty = f_pen + r_pen + tm_diff
                        if total_penalty < min_total_penalty:
                            min_total_penalty = total_penalty
                            best_pair = PrimerPair(
                                forward_seq=f_seq,
                                forward_tm=f_tm,
                                forward_gc=f_gc,
                                forward_start=f_start + 1,
                                forward_length=f_len,
                                reverse_seq=r_seq,
                                reverse_tm=r_tm,
                                reverse_gc=r_gc,
                                reverse_start=r_start + 1,
                                reverse_length=r_len,
                                product_size=prod_size,
                                pair_penalty=round(total_penalty, 2),
                                status="OK"
                            )

        if best_pair:
            return best_pair

        return PrimerPair(
            forward_seq="", forward_tm=0.0, forward_gc=0.0, forward_start=0, forward_length=0,
            reverse_seq="", reverse_tm=0.0, reverse_gc=0.0, reverse_start=0, reverse_length=0,
            product_size=0, status="REJECTED", error_reason="No pair met product size or Tm diff limits"
        )

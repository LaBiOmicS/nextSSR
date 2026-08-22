from typing import List
from nextssr.models import SSRItem, CompoundSSR, SequenceAnalysisResult
from nextssr.config import SSRConfig

class CompoundSSRProcessor:
    """Group individual SSRs into compound formations with Weber (1990) Compound Perfect vs Compound Disrupted classification."""
    def __init__(self, config: SSRConfig):
        self.config = config

    def process(self, result: SequenceAnalysisResult) -> SequenceAnalysisResult:
        """Analyze SSR list and group compound microsatellites."""
        if not result.ssrs:
            return result

        compounds: List[CompoundSSR] = []
        i = 0
        n = len(result.ssrs)

        while i < n:
            current_group = [result.ssrs[i]]
            j = i
            is_disrupted = False

            while j + 1 < n:
                dist = result.ssrs[j+1].start - result.ssrs[j].end - 1
                if 0 <= dist <= self.config.max_compound_distance:
                    if dist > 0:
                        is_disrupted = True
                    current_group.append(result.ssrs[j+1])
                    j += 1
                else:
                    break

            if len(current_group) > 1:
                pattern = " + ".join([f"({s.motif}){s.repeats}" for s in current_group])
                weber = "Compound Disrupted" if is_disrupted else "Compound Perfect"

                compounds.append(CompoundSSR(
                    seq_id=result.seq_id,
                    ssrs=current_group,
                    start=current_group[0].start,
                    end=current_group[-1].end,
                    compound_type="compound",
                    full_pattern=pattern,
                    weber_class=weber
                ))
            i = j + 1

        result.compounds = compounds
        return result

from math import floor
from pathlib import Path

import pysam
from loguru import logger

from markadoros.summary.results_models import ContigStats


class Contigs:
    """A contigs FASTA file, with assembly statistics computed on demand."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self._stats: ContigStats | None = None
        self._stats_computed = False

    def _get_stats(self) -> ContigStats | None:
        """Return assembly statistics, computing them on first use."""
        if not self._stats_computed:
            self._stats = self._calculate_stats()
            self._stats_computed = True
        return self._stats

    def _calculate_stats(self) -> ContigStats | None:
        count = 0
        lengths: list[int] = []
        with pysam.FastxFile(str(self.path)) as asm:
            for record in asm:
                count += 1
                if record.sequence is not None:
                    lengths.append(len(record.sequence))

        if not lengths:
            logger.warning("No contigs found in assembly file!")
            return None

        asm_size = sum(lengths)
        longest = max(lengths)

        lengths.sort(reverse=True)
        cumulative_sum = 0
        n50 = 0
        for length in lengths:
            cumulative_sum += length
            if cumulative_sum >= floor(asm_size / 2):
                n50 = length
                break

        return ContigStats(n=count, size=asm_size, n50=n50, longest=longest)

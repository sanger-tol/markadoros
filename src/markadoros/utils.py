from pathlib import Path

from Bio import Seq


def get_simple_name(input: Path) -> str:
    """Strip file extension, handling compound extensions like .fastq.gz"""
    if input.name.lower().endswith(".gz"):
        return Path(input.name[:-3]).stem

    return input.stem


def get_canonical_sequence(seq) -> str:
    """
    Get the canonical form of a sequence (lexicographically smaller of seq and its reverse complement).
    """
    seq_str = str(seq).upper()
    rev_comp = str(Seq(seq_str).reverse_complement())

    return min(seq_str, rev_comp)

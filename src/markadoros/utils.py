from pathlib import Path


def get_simple_name(input: Path) -> str:
    """Strip file extension, handling compound extensions like .fastq.gz"""
    if input.name.lower().endswith(".gz"):
        return Path(input.name[:-3]).stem

    return input.stem

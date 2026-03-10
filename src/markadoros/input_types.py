"""Input type handling for different sequencing platforms."""

from enum import Enum


class InputType(str, Enum):
    """Normalized input types for the search pipeline."""

    SHORT_READ = "short_read"
    RNASEQ = "rnaseq"
    PACBIO = "pacbio"
    ONT = "ont"
    CONTIGS = "contigs"


# Mapping from CLI aliases to normalized input types
INPUT_TYPE_ALIASES: dict[str, InputType] = {
    # Short read aliases
    "sr": InputType.SHORT_READ,
    "short": InputType.SHORT_READ,
    "illumina": InputType.SHORT_READ,
    # RNA-seq (kept separate as it uses different assembler settings)
    "rnaseq": InputType.RNASEQ,
    # PacBio aliases
    "pb": InputType.PACBIO,
    "pacbio": InputType.PACBIO,
    "pacbio_hifi": InputType.PACBIO,
    # Oxford Nanopore aliases
    "ont": InputType.ONT,
    "oxford_nanopore": InputType.ONT,
    "nanopore": InputType.ONT,
    # Pre-assembled contigs
    "contigs": InputType.CONTIGS,
}


def normalize_input_type(input_type: str) -> InputType:
    """Normalize a CLI input type string to an InputType enum.

    Args:
        input_type: Raw input type string from CLI.

    Returns:
        Normalized InputType enum value.

    Raises:
        ValueError: If the input type is not recognized.
    """
    input_type_lower = input_type.lower()

    if input_type_lower not in INPUT_TYPE_ALIASES:
        valid_types = ", ".join(sorted(INPUT_TYPE_ALIASES.keys()))
        raise ValueError(
            f"Unsupported input type: '{input_type}'. Valid types are: {valid_types}"
        )

    return INPUT_TYPE_ALIASES[input_type_lower]


def get_valid_input_types() -> list[str]:
    """Return a sorted list of all valid input type strings for CLI."""
    return sorted(INPUT_TYPE_ALIASES.keys())

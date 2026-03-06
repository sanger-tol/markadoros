import json
import os
import shutil
import subprocess
from math import floor
from pathlib import Path

import pysam
from jsonschema import ValidationError, validate
from loguru import logger


def get_simple_name(input: Path) -> str:
    """Strip file extension, handling compound extensions like .fastq.gz"""
    if input.name.lower().endswith(".gz"):
        return Path(input.name[:-3]).stem

    return input.stem


def validate_and_load_index(database: Path) -> dict:
    """
    Validate a database JSON file against the expected schema.

    Args:
        database: Path to the database.json file

    Returns:
        dict: The validated database configuration

    Raises:
        ValueError: If the file doesn't exist, is invalid JSON, or doesn't conform to schema
        ValidationError: If the JSON doesn't match the expected schema
    """
    if not database.exists():
        raise ValueError(f"{database} is not a valid database!")

    try:
        with open(database, "r") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in {database}: {e}")

    db_schema = {
        "type": "object",
        "additionalProperties": {
            "type": "object",
            "properties": {
                "marker": {"type": "string"},
                "min_length": {"type": "integer", "minimum": 1},
                "min_seq_id": {"type": "number", "minimum": 0, "maximum": 1},
                "min_aln_len": {"type": "integer", "minimum": 1},
                "db": {"type": "string"},
            },
            "required": ["marker", "min_length", "min_seq_id", "min_aln_len", "db"],
        },
    }

    try:
        validate(instance=data, schema=db_schema)
    except ValidationError as e:
        raise ValueError(
            f"Database configuration in {database} does not match expected schema: {e.message}"
        ) from e

    return data


def _check_mmseqs_binary(mmseqs_path: str) -> None:
    """Check that mmseqs2 binary exists and is executable.

    Args:
        mmseqs_path: Path to the mmseqs2 binary

    Raises:
        RuntimeError: If binary cannot be executed
    """
    try:
        result = subprocess.run(
            [mmseqs_path],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"mmseqs2 at {mmseqs_path} returned exit code {result.returncode}"
            )
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"mmseqs2 binary check timed out at {mmseqs_path}")
    except FileNotFoundError:
        raise RuntimeError(f"mmseqs2 binary not found at {mmseqs_path}")
    except Exception as e:
        raise RuntimeError(f"Failed to check mmseqs2 binary at {mmseqs_path}: {e}")


def set_mmseqs_path() -> None:
    """
    Attempt to find mmseqs2 binary on PATH and set MMSEQS2_PATH environment variable.

    Raises RuntimeError if mmseqs2 cannot be found or is not executable.
    """
    # First check if MMSEQS2_PATH is already set
    existing_path = os.getenv("MMSEQS2_PATH")

    if existing_path and os.path.exists(existing_path):
        try:
            _check_mmseqs_binary(existing_path)
            logger.info(f"Using pre-configured mmseqs2 at {existing_path}")
        except RuntimeError as e:
            logger.error(f"Error with pre-configured mmseqs2: {e}")
            raise
        return

    # Try to find mmseqs2 on PATH
    mmseqs_binary = shutil.which("mmseqs")

    if mmseqs_binary:
        try:
            _check_mmseqs_binary(mmseqs_binary)
            # Set the environment variable for pymmseqs to use
            os.environ["MMSEQS2_PATH"] = mmseqs_binary
            logger.info(f"Found mmseqs2 at {mmseqs_binary}")
        except RuntimeError as e:
            logger.error(f"Error checking mmseqs2 binary: {e}")
            raise
        return

    # Fallback: raise error if we couldn't find mmseqs2
    raise RuntimeError(
        "Could not find mmseqs2 binary on PATH or MMSEQS2_PATH. "
        "Please install mmseqs2 or set the MMSEQS2_PATH environment variable."
    )


def calculate_contig_statistics(contigs: Path) -> dict[str, int]:
    """Calculate assembly statistics from contigs.

    Args:
        contigs: Path to contigs FASTA file

    Returns:
        Dictionary with assembly statistics (n, size, n50)
    """
    count = 0
    lengths = []
    with pysam.FastxFile(str(contigs)) as asm:
        for record in asm:
            count += 1
            if record.sequence is not None:
                lengths.append(len(record.sequence))

    asm_size = sum(lengths)

    lengths.sort(reverse=True)
    cumulative_sum = 0
    n50 = 0
    for length in lengths:
        cumulative_sum += length
        if cumulative_sum >= floor(asm_size / 2):
            n50 = length
            break

    return {"n": count, "size": asm_size, "n50": n50, "longest": max(lengths)}

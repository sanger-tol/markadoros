import json
from pathlib import Path

from jsonschema import ValidationError, validate


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

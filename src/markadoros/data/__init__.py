"""Data loading utilities for markadoros."""

import json
from pathlib import Path

# Path to the bold.json file
_DATA_DIR = Path(__file__).parent
BOLD_CONFIG_PATH = _DATA_DIR / "bold.json"
UNITE_CONFIG_PATH = _DATA_DIR / "unite.json"


def load_config(config_path: Path) -> dict:
    """Load and return the BOLD database configuration.

    Returns:
        dict: Dictionary containing BOLD marker configurations.
              Keys are marker names (e.g., "BOLD_COI", "BOLD_CYTB", etc.)
              Values are dictionaries with marker, min_length, min_seq_id, and min_aln_len.

    Raises:
        FileNotFoundError: If bold.json cannot be found.
        json.JSONDecodeError: If bold.json is malformed.
    """
    if not config_path.exists():
        raise FileNotFoundError(f"BOLD configuration file not found: {config_path}")

    with open(config_path, "r") as f:
        return json.load(f)


BOLD_CONFIG = load_config(BOLD_CONFIG_PATH)
UNITE_CONFIG = load_config(UNITE_CONFIG_PATH)

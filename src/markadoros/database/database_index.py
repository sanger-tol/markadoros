import json
from pathlib import Path


class DatabaseIndex:
    """Manages JSON metadata index for marker databases."""

    def __init__(self, index_path: Path):
        self._index_path = index_path

    def _load(self) -> dict:
        """Load the database index from disk."""
        try:
            with open(self._index_path, "r") as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            return {}

    def _save(self, data: dict) -> None:
        """Save the database index to disk."""
        with open(self._index_path, "w") as f:
            json.dump(data, f, indent=2)

    def _convert_paths(self, entry: dict) -> dict:
        """Recursively convert all Path objects to strings in a dictionary."""
        converted = {}
        for key, value in entry.items():
            if isinstance(value, Path):
                converted[key] = str(value.resolve())
            elif isinstance(value, dict):
                converted[key] = self._convert_paths(value)
            elif isinstance(value, list):
                converted[key] = [
                    str(v.resolve()) if isinstance(v, Path) else v for v in value
                ]
            else:
                converted[key] = value

        return converted

    def add_entry(self, marker_name: str, db_entry: dict) -> None:
        """Add or update an entry in the database index."""
        db_dict = self._load()
        db_dict[marker_name] = self._convert_paths(db_entry)
        self._save(db_dict)

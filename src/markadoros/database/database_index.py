from pathlib import Path

from pydantic import BaseModel, Field, RootModel, field_validator


class DatabaseEntry(BaseModel):
    """Metadata for one marker database recorded in the index."""

    marker: str
    db: Path
    taxon_db: Path
    n_seqs: int | None = Field(default=None, ge=0)
    built_from: Path | None = None
    deduplicated: bool | None = None
    clustered: bool | None = None
    bins_taxa: Path | None = None

    @field_validator("db", "taxon_db", "built_from", "bins_taxa")
    @classmethod
    def _resolve_paths(cls, value: Path | None) -> Path | None:
        return value.resolve() if value is not None else None


class DatabaseIndexFile(RootModel[dict[str, DatabaseEntry]]):
    """The full db.json index, keyed by database name (e.g. BOLD.COI)."""


class DatabaseIndex:
    """Manages the validated JSON metadata index for marker databases."""

    def __init__(self, index_path: Path):
        self._index_path = index_path

    def load(self) -> dict[str, DatabaseEntry]:
        """Load and validate the index. Raises if the file is missing or invalid."""
        return DatabaseIndexFile.model_validate_json(self._index_path.read_text()).root

    def add_entry(self, name: str, entry: DatabaseEntry) -> None:
        """Add or update an entry, keeping any existing entries."""
        entries = self.load() if self._index_path.exists() else {}
        entries[name] = entry
        self._index_path.write_text(
            DatabaseIndexFile(entries).model_dump_json(indent=2)
        )

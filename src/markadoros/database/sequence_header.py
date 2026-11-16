from pydantic import BaseModel


class SequenceHeader(BaseModel):
    """A processed database FASTA header: ``seq_id|marker|taxon|lineage``."""

    seq_id: str
    marker: str
    taxon: str
    lineage: str

    @classmethod
    def from_header(cls, header: str) -> "SequenceHeader":
        fields = header.split("|", 3)
        if len(fields) != 4:
            raise ValueError(
                f"Expected 4 pipe-separated fields, got {len(fields)}: {header[:100]}"
            )

        seq_id, marker, taxon, lineage = fields
        return cls(
            seq_id=seq_id,
            marker=marker,
            taxon=taxon.replace("_", " "),
            lineage=lineage,
        )

    def to_header(self) -> str:
        return f"{self.seq_id}|{self.marker}|{self.taxon.replace(' ', '_')}|{self.lineage.replace(' ', '_')}"

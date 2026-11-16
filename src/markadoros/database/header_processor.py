from markadoros.database.sequence_header import SequenceHeader


def process_bold_header(header: str) -> SequenceHeader:
    """Process a BOLD FASTA header into a SequenceHeader."""
    split_header = header.split("|")

    if len(split_header) != 4:
        raise ValueError(
            f"Invalid BOLD header: expected 4 pipe-separated fields, got {len(split_header)}. "
            f"Header: {header[:100]}..."
        )

    if "," not in split_header[3]:
        raise ValueError(
            f"Invalid BOLD header: lineage field (4th) should contain commas. "
            f"Header: {header[:100]}..."
        )

    seq_id = split_header[0]
    marker = split_header[1]
    lineage = split_header[3]

    lineage_split = lineage.split(",")
    filtered = [x for x in lineage_split if x != "None"]

    if lineage_split[-1] != "None":
        taxon = lineage_split[-2]
    else:
        taxon = filtered[-1] if filtered else "NA"

    return SequenceHeader(seq_id=seq_id, marker=marker, taxon=taxon, lineage=lineage)


def process_unite_header(header: str) -> SequenceHeader:
    """Process a UNITE FASTA header into a SequenceHeader."""
    split_header = header.split("|")

    if not len(split_header) == 5 or ";" not in split_header[4]:
        raise ValueError(f"Invalid UNITE header: {header}")

    return SequenceHeader(
        seq_id=split_header[1] + "/" + split_header[2],
        marker="ITS",
        taxon=split_header[0],
        lineage=split_header[4],
    )


def process_silva_lsu_header(header: str) -> SequenceHeader:
    """Process a SILVA LSU FASTA header into a SequenceHeader."""
    split_header = header.split(" ", 1)
    lineage = split_header[1]
    taxon = [x for x in lineage.split(";") if x is not None][-1]

    return SequenceHeader(
        seq_id=split_header[0], marker="LSU", taxon=taxon, lineage=lineage
    )


def process_silva_ssu_header(header: str) -> SequenceHeader:
    """Process a SILVA SSU FASTA header into a SequenceHeader."""
    split_header = header.split(" ", 1)
    lineage = split_header[1]
    taxon = lineage.split(";")[-1]

    return SequenceHeader(
        seq_id=split_header[0], marker="SSU", taxon=taxon, lineage=lineage
    )


def process_generic_header(header: str) -> SequenceHeader:
    """Process a generic FASTA header into a SequenceHeader."""
    if not len(header.split("|")) == 4:
        raise ValueError(
            f"Invalid header: {header}. Expected header with format '><seq_id>|<marker>|<taxon name>|<lineage>'"
        )

    return SequenceHeader.from_header(header)

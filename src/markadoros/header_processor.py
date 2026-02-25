def process_bold_header(header: str) -> tuple[str, str]:
    """Process a BOLD FASTA header and return the new FASTA header."""
    split_header = header.split("|")

    if not len(split_header) == 4 or "," not in split_header[3]:
        raise ValueError(f"Invalid BOLD header: {header}")

    seq_id = split_header[0]
    marker = split_header[1]
    lineage = split_header[3]
    taxon = [x for x in lineage.split(",") if x != "None"][-1]

    return marker, "|".join([seq_id, marker, taxon, lineage])


def process_unite_header(header: str) -> tuple[str, str]:
    """Process a UNITE FASTA header and return the new FASTA header."""
    split_header = header.split("|")

    if not len(split_header) == 5 or ";" not in split_header[4]:
        raise ValueError(f"Invalid UNITE header: {header}")

    seq_id = split_header[1] + ";" + split_header[2]
    marker = "ITS"
    lineage = split_header[4]
    taxon = split_header[0]

    return marker, "|".join([seq_id, marker, taxon, lineage])


def process_generic_header(header: str) -> tuple[str, str]:
    """Process a generic FASTA header and return the new FASTA header."""
    split_header = header.split("|")

    if not len(split_header) == 4:
        raise ValueError(
            f"Invalid header: {header}. Expected header with format '><seq_id>|<marker>|<taxon name>|<lineage>'"
        )

    marker = split_header[1]

    return marker, header

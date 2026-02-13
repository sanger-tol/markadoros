"""markadoros - A tool for barcode processing."""

from markadoros.cli import cli

VALID_MARKERS = {
    "COI": "COI-5P",
    "CYTB": "CYTB",
    "rbcL": "rbcL",
    "matK": "matK",
    "18S": "18S",
    "28S": "28S",
}

__all__ = ["cli", "VALID_MARKERS"]

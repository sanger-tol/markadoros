from pathlib import Path

import polars as pl
from isal import igzip
from loguru import logger


class BoldTSVProcessor:
    def __init__(self, tmpdir: Path):
        self.tmpdir = tmpdir
        if not self.tmpdir.exists():
            self.tmpdir.mkdir(parents=True)

        self.taxonomy_levels = [
            "subspecies",
            "species",
            "genus",
            "tribe",
            "subfamily",
            "family",
            "order",
            "class",
            "phylum",
            "kingdom",
        ]

    def select_best_row(self, group_df, taxonomy_levels):
        """For each group, find first row matching modal (most common) taxonomy at any level.

        Starts from finest taxonomic resolution and falls back to coarser levels if all are null.
        """
        for level in taxonomy_levels:
            non_null = group_df.filter(pl.col(level).is_not_null())

            if non_null.height == 0:
                continue

            # Count occurrences and sort to find modal value
            modal_df = (
                non_null.group_by(level)
                .agg(pl.len().alias("count"))
                .sort("count", descending=True)
                .head(1)
            )

            if modal_df.height > 0:
                modal_value = modal_df.select(level).item()
                match = non_null.filter(pl.col(level) == modal_value)
                if match.height > 0:
                    return match.head(1)

        # Fallback: return first row
        return group_df.head(1)

    def process_to_fasta(self, tsv: Path):
        selected_columns = [
            "processid",
            "bin_uri",
            "marker_code",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "subfamily",
            "tribe",
            "genus",
            "species",
            "subspecies",
            "nuc",
        ]

        schema = (
            pl.scan_csv(tsv, separator="\t", has_header=True, null_values=["None", ""])
            .select(selected_columns)
            .collect_schema()
        )

        logger.info(f"Processing BOLD TSV file with Polars: {tsv}")
        df = (
            pl.scan_csv(
                tsv,
                separator="\t",
                has_header=True,
                null_values=["None", ""],
                quote_char=None,
            )
            .select(selected_columns)
            .filter(
                pl.col("bin_uri").is_not_null()
                & pl.col("nuc").is_not_null()
                & (pl.col("marker_code") == "COI-5P")
            )
            .group_by("bin_uri")
            .map_groups(
                lambda group: self.select_best_row(group, self.taxonomy_levels),
                schema=schema,
            )
            .collect(engine="streaming")
        )
        ## Make the type-checker happy
        assert isinstance(df, pl.DataFrame)

        outfile = Path(self.tmpdir / "BOLD.BINS.fasta.gz")
        logger.info(f"Processing {df.height} records to FASTA: {outfile}")
        with igzip.open(outfile, "wb") as fout:
            collected = []
            for row in df.iter_rows(named=True):
                if row["nuc"] is None or len(row["nuc"]) == 0:
                    continue

                if row["processid"] is None or row["bin_uri"] is None:
                    continue

                # Reverse to go from kingdom (coarse) to subspecies (fine)
                lineage = [
                    row.get(level) or "None" for level in reversed(self.taxonomy_levels)
                ]
                lineage_str = ",".join(lineage) if lineage else "None"

                header = f">{row['processid']}/{row['bin_uri']}|{row['marker_code']}|Ignored|{lineage_str}\n"
                nuc = row.get("nuc") or ""

                collected.append(header.encode("utf-8") + nuc.encode("utf-8") + b"\n")

                if len(collected) >= 100_000:
                    fout.writelines(collected)
                    collected.clear()

            if collected:
                fout.writelines(collected)

        return outfile

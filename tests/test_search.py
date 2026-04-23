"""Integration tests for search command."""

from pathlib import Path

import pytest
from click.testing import CliRunner

from markadoros.cli import cli


@pytest.fixture
def cli_runner():
    """Return a Click CLI test runner."""
    return CliRunner()


class TestSearchCommand:
    """Tests for the search command with various input types."""

    def test_search_hic_short_reads_sr(
        self, cli_runner, test_data_dir, test_db_dir, tmp_path
    ):
        """Test search with Hi-C short reads using -x sr and -n 1000."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            db_index = test_db_dir / "db.json"

            result = cli_runner.invoke(
                cli,
                [
                    "search",
                    "-x",
                    "sr",
                    "-i",
                    str(db_index),
                    "-o",
                    str(outdir),
                    "-n",
                    "1000",
                    "--db",
                    "BOLD.COI",
                    str(Path(test_data_dir / "icHalSede1.hic.cram")),
                ],
            )

            summary_file = Path(outdir) / "icHalSede1.hic.COI.summary.json"

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert summary_file.exists()

    def test_search_hic_short_reads_sr_expected_taxon(
        self, cli_runner, test_data_dir, test_db_dir, tmp_path
    ):
        """Test search with Hi-C short reads using -x sr and -n 1000."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            db_index = test_db_dir / "db.json"

            result = cli_runner.invoke(
                cli,
                [
                    "search",
                    "-x",
                    "sr",
                    "-i",
                    str(db_index),
                    "-o",
                    str(outdir),
                    "-n",
                    "1000",
                    "--db",
                    "BOLD.COI",
                    "--expected-taxon",
                    "'Halyzia sedeimguttata'",
                    str(Path(test_data_dir / "icHalSede1.hic.cram")),
                ],
            )

            summary_file = Path(outdir) / "icHalSede1.hic.COI.summary.json"

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert summary_file.exists()

    def test_search_hic_short_reads_sr_expected_taxon_synonyms(
        self, cli_runner, test_data_dir, test_db_dir, tmp_path
    ):
        """Test search with Hi-C short reads using -x sr and -n 1000."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            db_index = test_db_dir / "db.json"

            result = cli_runner.invoke(
                cli,
                [
                    "search",
                    "-x",
                    "sr",
                    "-i",
                    str(db_index),
                    "-o",
                    str(outdir),
                    "-n",
                    "1000",
                    "--db",
                    "BOLD.COI",
                    "--expected-taxon",
                    "'Halyzia sedeimguttata'",
                    "--synonyms",
                    "'Melanostoma mellinum'",
                    str(Path(test_data_dir / "icHalSede1.hic.cram")),
                ],
            )

            summary_file = Path(outdir) / "icHalSede1.hic.COI.summary.json"

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert summary_file.exists()

    def test_search_pacbio_with_20_reads(
        self, cli_runner, test_data_dir, test_db_dir, tmp_path
    ):
        """Test search with PacBio reads using -n 20 reads."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            db_index = test_db_dir / "db.json"

            result = cli_runner.invoke(
                cli,
                [
                    "search",
                    "-x",
                    "sr",
                    "-i",
                    str(db_index),
                    "-o",
                    str(outdir),
                    "-n",
                    "20",
                    "--db",
                    "BOLD.COI",
                    str(Path(test_data_dir / "icHalSede1.pacbio.subset.fa.gz")),
                ],
            )

            summary_file = Path(outdir) / "icHalSede1.pacbio.subset.COI.summary.json"

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert summary_file.exists()

    def test_search_pacbio_with_1_read_fallback(
        self, cli_runner, test_data_dir, test_db_dir, tmp_path
    ):
        """Test search with PacBio reads using -n 1 read to test fallback mechanism."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            db_index = test_db_dir / "db.json"

            result = cli_runner.invoke(
                cli,
                [
                    "search",
                    "-x",
                    "sr",
                    "-i",
                    str(db_index),
                    "-o",
                    str(outdir),
                    "-n",
                    "1",
                    "--db",
                    "BOLD.COI",
                    str(Path(test_data_dir / "icHalSede1.pacbio.subset.fa.gz")),
                ],
            )

            summary_file = Path(outdir) / "icHalSede1.pacbio.subset.COI.summary.json"

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert summary_file.exists()

    def test_search_pacbio_no_hits(
        self, cli_runner, test_data_dir, test_db_dir, tmp_path
    ):
        """Test search with PacBio reads that have no hits in the database."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            db_index = test_db_dir / "db.json"

            result = cli_runner.invoke(
                cli,
                [
                    "search",
                    "-x",
                    "sr",
                    "-i",
                    str(db_index),
                    "-o",
                    str(outdir),
                    "-n",
                    "20",
                    "--db",
                    "BOLD.COI",
                    str(Path(test_data_dir / "icHalSede1.pacbio.nohit.fa.gz")),
                ],
            )

            summary_file = Path(outdir) / "icHalSede1.pacbio.nohit.COI.summary.json"

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert not summary_file.exists()

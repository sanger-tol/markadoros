"""Integration tests for database creation commands."""

from pathlib import Path

import pytest
from click.testing import CliRunner

from markadoros.cli import cli


@pytest.fixture
def cli_runner():
    """Return a Click CLI test runner."""
    return CliRunner()


def assert_db_files(outdir, prefix, marker):
    if not isinstance(outdir, Path):
        outdir = Path(outdir)

    subdir = f"{prefix}.{marker}"

    assert Path(outdir / "db.json").exists(), "Database index not created"
    assert Path(outdir / subdir / "db").exists(), "Database file not created"
    assert Path(outdir / subdir / "db.index").exists(), "Database file not created"
    assert Path(outdir / subdir / "db.lookup").exists(), "Database file not created"
    assert Path(outdir / subdir / "db.source").exists(), "Database file not created"
    assert Path(outdir / subdir / "db_h").exists(), "Database file not created"
    assert Path(outdir / subdir / "db_h.dbtype").exists(), "Database file not created"
    assert Path(outdir / subdir / "db_h.index").exists(), "Database file not created"
    assert Path(outdir / subdir / "db.dbtype").exists(), "Database file not created"
    assert Path(outdir / subdir / "taxon.json.gz").exists(), "Taxon JSON not created"


class TestBOLDFastaDatabase:
    """Tests for building databases from BOLD FASTA files."""

    def test_bold_fasta_database_creation(self, cli_runner, test_data_dir, tmp_path):
        """Test basic BOLD FASTA database creation with COI marker."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            result = cli_runner.invoke(
                cli,
                [
                    "database",
                    "-x",
                    "bold",
                    "--prefix",
                    "BOLD",
                    "--outdir",
                    str(outdir),
                    str(Path(test_data_dir / "BOLD.sub.fa.gz")),
                ],
            )

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert_db_files(outdir, "BOLD", "COI")

    def test_bold_fasta_with_exclude_file(self, cli_runner, test_data_dir, tmp_path):
        """Test BOLD FASTA database creation with exclude patterns."""
        exclude_file = test_data_dir / "wolbachia.exclude"

        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            result = cli_runner.invoke(
                cli,
                [
                    "database",
                    "-x",
                    "bold",
                    "--prefix",
                    "BOLD",
                    "--exclude-file",
                    str(exclude_file),
                    "--outdir",
                    str(outdir),
                    str(Path(test_data_dir / "BOLD.sub.fa.gz")),
                ],
            )

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert_db_files(outdir, "BOLD", "COI")

    def test_bold_fasta_with_min_length(self, cli_runner, test_data_dir, tmp_path):
        """Test BOLD FASTA database creation with minimum length filter."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            result = cli_runner.invoke(
                cli,
                [
                    "database",
                    "-x",
                    "bold",
                    "--prefix",
                    "BOLD",
                    "--min-length",
                    "500",
                    "--outdir",
                    str(outdir),
                    str(Path(test_data_dir / "BOLD.sub.fa.gz")),
                ],
            )

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert_db_files(outdir, "BOLD", "COI")

    def test_bold_fasta_with_create_index(self, cli_runner, test_data_dir, tmp_path):
        """Test BOLD FASTA database creation with index creation."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            result = cli_runner.invoke(
                cli,
                [
                    "database",
                    "-x",
                    "bold",
                    "--prefix",
                    "BOLD",
                    "--create-index",
                    "--outdir",
                    str(outdir),
                    str(Path(test_data_dir / "BOLD.sub.fa.gz")),
                ],
            )

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert_db_files(outdir, "BOLD", "COI")

    def test_bold_fasta_no_cleanup(self, cli_runner, test_data_dir, tmp_path):
        """Test BOLD FASTA database creation without cleanup of temporary files."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            result = cli_runner.invoke(
                cli,
                [
                    "database",
                    "-x",
                    "bold",
                    "--prefix",
                    "BOLD",
                    "--no-cleanup",
                    "--outdir",
                    str(outdir),
                    str(Path(test_data_dir / "BOLD.sub.fa.gz")),
                ],
            )

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert_db_files(outdir, "BOLD", "COI")

    def test_unite_fasta_database_creation(self, cli_runner, test_data_dir, tmp_path):
        """Test basic BOLD FASTA database creation with COI marker."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            result = cli_runner.invoke(
                cli,
                [
                    "database",
                    "-x",
                    "unite",
                    "--prefix",
                    "UNITE",
                    "--outdir",
                    str(outdir),
                    str(Path(test_data_dir / "UNITE.sub.fa.gz")),
                ],
            )

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert_db_files(outdir, "UNITE", "ITS")

    def test_silva_ssu_fasta_database_creation(
        self, cli_runner, test_data_dir, tmp_path
    ):
        """Test basic BOLD FASTA database creation with COI marker."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            result = cli_runner.invoke(
                cli,
                [
                    "database",
                    "-x",
                    "silva_ssu",
                    "--prefix",
                    "SILVA",
                    "--outdir",
                    str(outdir),
                    str(Path(test_data_dir / "SILVA.sub.fa.gz")),
                ],
            )

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert_db_files(outdir, "SILVA", "SSU")

    def test_silva_lsu_fasta_database_creation(
        self, cli_runner, test_data_dir, tmp_path
    ):
        """Test basic BOLD FASTA database creation with COI marker."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            result = cli_runner.invoke(
                cli,
                [
                    "database",
                    "-x",
                    "silva_lsu",
                    "--prefix",
                    "SILVA",
                    "--outdir",
                    str(outdir),
                    str(Path(test_data_dir / "SILVA.sub.fa.gz")),
                ],
            )

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert_db_files(outdir, "SILVA", "LSU")


class TestBOLDTSVDatabase:
    """Tests for building databases from BOLD TSV files."""

    def test_bold_tsv_database_creation(self, cli_runner, test_data_dir, tmp_path):
        """Test basic BOLD TSV database creation."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            result = cli_runner.invoke(
                cli,
                [
                    "bold-coi-from-tsv",
                    "--prefix",
                    "BOLD",
                    "--outdir",
                    str(outdir),
                    str(Path(test_data_dir / "BOLD.sub.tsv")),
                ],
            )

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert_db_files(outdir, "BOLD", "COI")

    def test_bold_tsv_with_exclude_file(self, cli_runner, test_data_dir, tmp_path):
        """Test BOLD TSV database creation with exclude patterns."""
        exclude_file = test_data_dir / "wolbachia.exclude"

        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            result = cli_runner.invoke(
                cli,
                [
                    "bold-coi-from-tsv",
                    "--prefix",
                    "BOLD",
                    "--exclude-file",
                    str(exclude_file),
                    "--outdir",
                    str(outdir),
                    str(Path(test_data_dir / "BOLD.sub.tsv")),
                ],
            )

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert_db_files(outdir, "BOLD", "COI")

    def test_bold_tsv_with_min_length(self, cli_runner, test_data_dir, tmp_path):
        """Test BOLD TSV database creation with minimum length filter."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            result = cli_runner.invoke(
                cli,
                [
                    "bold-coi-from-tsv",
                    "--prefix",
                    "BOLD",
                    "--min-length",
                    "500",
                    "--outdir",
                    str(outdir),
                    str(Path(test_data_dir / "BOLD.sub.tsv")),
                ],
            )

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert_db_files(outdir, "BOLD", "COI")

    def test_bold_tsv_with_create_index(self, cli_runner, test_data_dir, tmp_path):
        """Test BOLD TSV database creation with index creation."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            result = cli_runner.invoke(
                cli,
                [
                    "bold-coi-from-tsv",
                    "--prefix",
                    "BOLD",
                    "--create-index",
                    "--outdir",
                    str(outdir),
                    str(Path(test_data_dir / "BOLD.sub.tsv")),
                ],
            )

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert_db_files(outdir, "BOLD", "COI")

    def test_bold_tsv_no_cleanup(self, cli_runner, test_data_dir, tmp_path):
        """Test BOLD TSV database creation without cleanup of temporary files."""
        with cli_runner.isolated_filesystem(temp_dir=tmp_path) as outdir:
            result = cli_runner.invoke(
                cli,
                [
                    "bold-coi-from-tsv",
                    "--prefix",
                    "BOLD",
                    "--no-cleanup",
                    "--outdir",
                    str(outdir),
                    str(Path(test_data_dir / "BOLD.sub.tsv")),
                ],
            )

            assert result.exit_code == 0, f"Command failed with: {result.output}"
            assert_db_files(outdir, "BOLD", "COI")

import json
import shutil
import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def test_data_dir():
    """Return the path to the test data directory."""
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def session_tmp_dir():
    """Create a temporary directory for the entire test session."""
    tmpdir = tempfile.mkdtemp(prefix="markadoros_test_")
    yield Path(tmpdir)
    # Cleanup after session
    shutil.rmtree(tmpdir, ignore_errors=True)


@pytest.fixture(scope="session")
def generated_test_db_dir(session_tmp_dir):
    """Generate test database from BOLD.sub.fa.gz using the database command.

    This ensures the database is created fresh with correct paths for the
    current system, avoiding hardcoded absolute paths that fail on CI.
    """
    from click.testing import CliRunner

    from markadoros.cli import cli

    test_data_dir = Path(__file__).parent / "data"
    db_dir = session_tmp_dir / "test_db"
    db_dir.mkdir(parents=True, exist_ok=True)

    runner = CliRunner()

    # Create database from BOLD.sub.fa.gz
    result = runner.invoke(
        cli,
        [
            "database",
            "-x",
            "bold",
            "--prefix",
            "BOLD",
            "--outdir",
            str(db_dir),
            str(test_data_dir / "BOLD.sub.fa.gz"),
        ],
    )

    if result.exit_code != 0:
        raise RuntimeError(f"Failed to create test database: {result.output}")

    return db_dir


@pytest.fixture
def test_db_dir(generated_test_db_dir):
    """Return the path to the test database directory.

    This fixture now uses the dynamically generated database instead of
    a hardcoded one, ensuring it works on any system.
    """
    return generated_test_db_dir


@pytest.fixture
def db_index_path(temp_db_dir):
    """Return the path to the database index JSON file."""
    return temp_db_dir / "db.json"


@pytest.fixture
def cli_runner():
    """Return a Click CLI test runner."""
    from click.testing import CliRunner

    return CliRunner()

from pathlib import Path

import pytest


@pytest.fixture
def test_data_dir():
    """Return the path to the test data directory."""
    return Path(__file__).parent / "data"


@pytest.fixture
def test_db_dir():
    """Return the path to the test data directory."""
    return Path(__file__).parent / "db"


@pytest.fixture
def db_index_path(temp_db_dir):
    """Return the path to the database index JSON file."""
    return temp_db_dir / "db.json"


@pytest.fixture
def cli_runner():
    """Return a Click CLI test runner."""
    from click.testing import CliRunner

    return CliRunner()

"""Tests for ipral_file_infos.py."""
import subprocess
from pathlib import Path

import pytest

MAIN_DIR = Path(__file__).parent.parent
TEST_DIR = MAIN_DIR / "test"
EXE = MAIN_DIR / "ipral_file_infos.py"


@pytest.mark.parametrize("date_start, date_end", [("2020-12-01", "2020-12-31")])
def test_script(date_start, date_end):
    """Test if script is working."""
    ret = subprocess.call([EXE, date_start, date_end])

    assert ret == 0

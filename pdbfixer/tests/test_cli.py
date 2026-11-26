# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Test command-line interface.

"""

# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import os
import shlex
import subprocess

import pytest

# =============================================================================================
# UNIT TESTS
# =============================================================================================


def run_cli(arguments, expected_output=None):
    output = subprocess.check_output(["pdbfixer", *shlex.split(arguments)], text=True)
    if expected_output is not None:
        assert output == expected_output


def test_help():
    run_cli("--help")


def test_pdbid():
    run_cli("--pdbid 1LE1")


@pytest.mark.skipif(
    os.getenv("GITHUB_ACTION") is not None, reason="Cannot download during CI"
)
def test_url():
    run_cli(
        '--url "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=1LE1"'
    )

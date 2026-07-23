"""Keep the package version single-sourced and current.

The version appears in two files (``pyproject.toml`` and ``msep/__init__.py``) and is
also quoted in the manuscript. These tests fail if the two files disagree, so a bump
cannot be applied to one and forgotten in the other.
"""

import re
from pathlib import Path

import msep

REPO_ROOT = Path(__file__).resolve().parents[1]
EXPECTED_VERSION = "1.2.0"


def _pyproject_version() -> str:
    text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    assert match, "no version field in pyproject.toml"
    return match.group(1)


def test_dunder_version_matches_expected() -> None:
    assert msep.__version__ == EXPECTED_VERSION


def test_pyproject_matches_dunder() -> None:
    assert _pyproject_version() == msep.__version__

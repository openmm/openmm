"""Pytest hooks and defaults for UMA ice RPMD tests."""

from __future__ import annotations

from pathlib import Path

import pytest

# Executable-style modules named test_*.py that run on import (not pytest tests).
collect_ignore = [
    "test_all_models.py",
    "test_bead_scaling.py",
    "test_esen_large.py",
    "test_esen_openmm.py",
    "test_parallel_uma.py",
    "test_true_batch_rpmd.py",
]


def pytest_collection_modifyitems(config, items) -> None:
    """Tag all tests collected from this tree with ``ice_rpmd`` (unless already marked)."""
    ice_root = Path(__file__).resolve().parent
    for item in items:
        path = getattr(item, "path", None)
        if path is None:
            continue
        p = Path(path).resolve()
        if ice_root == p.parent or ice_root in p.parents:
            if item.get_closest_marker("ice_rpmd") is None:
                item.add_marker(pytest.mark.ice_rpmd)

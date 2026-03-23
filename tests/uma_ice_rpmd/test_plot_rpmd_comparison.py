"""Tests for plot_rpmd_comparison CSV parsing (kinetic temperature column)."""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

from plot_rpmd_comparison import _read_csv

_HEADER_FULL = (
    "step,time_ps,T_K,PE_kj_mol,"
    "q6_mean,q6_std,q6_p10,q6_p50,q6_p90,q6_min,q6_max,"
    "q_tet_mean,q_tet_std,q_tet_p10,q_tet_p50,q_tet_p90,q_tet_min,q_tet_max,"
    "n_q6_valid,n_q_tet_valid,n_oxygen\n"
)


def test_read_csv_parses_temperature_and_nan_for_empty_cells() -> None:
    content = (
        _HEADER_FULL
        + "0,0.0,,,0.1,0,0,0,0,0,0,0.2,0,0,0,0,0,0,1,1,1\n"
        + "1,0.001,240.0,-1.0,0.1,0,0,0,0,0,0,0.2,0,0,0,0,0,0,1,1,1\n"
    )
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, encoding="utf-8") as f:
        f.write(content)
        path = Path(f.name)
    try:
        t, q6, qt, tk = _read_csv(path)
        assert t.shape == (2,)
        assert np.isnan(tk[0])  # type: ignore[union-attr]
        assert tk is not None
        assert tk[1] == pytest.approx(240.0)
    finally:
        path.unlink(missing_ok=True)


def test_read_csv_returns_none_temperature_when_column_absent() -> None:
    content = (
        "step,time_ps,q6_mean,q_tet_mean\n"
        "0,0.0,0.1,0.2\n"
    )
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, encoding="utf-8") as f:
        f.write(content)
        path = Path(f.name)
    try:
        t, q6, qt, tk = _read_csv(path)
        assert tk is None
    finally:
        path.unlink(missing_ok=True)

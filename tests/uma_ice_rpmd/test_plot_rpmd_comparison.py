"""Tests for plot_rpmd_comparison CSV parsing (kinetic temperature column)."""
from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path

import numpy as np
import pytest

from scripts.plot_rpmd_comparison import (
    _enrich_title_with_system,
    _format_title_two_lines,
    _read_csv,
    _read_n_oxygen_from_csv,
    _save_cv_plane_figure,
)

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


def test_format_title_preserves_explicit_newline() -> None:
    s = "Line one\nLine two"
    assert _format_title_two_lines(s) == s


def test_format_title_splits_on_em_dash() -> None:
    s = "RPMD ice order: A vs B — 243 K, PILE-G"
    assert _format_title_two_lines(s) == "RPMD ice order: A vs B\n243 K, PILE-G"


def test_read_n_oxygen_from_csv() -> None:
    content = (
        _HEADER_FULL
        + "0,0.0,,,0.1,0,0,0,0,0,0,0.2,0,0,0,0,0,0,1,1,32\n"
    )
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, encoding="utf-8") as f:
        f.write(content)
        path = Path(f.name)
    try:
        assert _read_n_oxygen_from_csv(path) == 32
    finally:
        path.unlink(missing_ok=True)


def test_read_n_oxygen_missing_column() -> None:
    content = "step,time_ps,q6_mean,q_tet_mean\n0,0.0,0.1,0.2\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, encoding="utf-8") as f:
        f.write(content)
        path = Path(f.name)
    try:
        assert _read_n_oxygen_from_csv(path) is None
    finally:
        path.unlink(missing_ok=True)


def test_enrich_title_appends_molecules_and_beads() -> None:
    t = "A\nB"
    assert _enrich_title_with_system(t, 32, 32) == "A\nB · 32 molecules, 32 beads"


def test_enrich_title_partial() -> None:
    assert _enrich_title_with_system("Only mol", 16, None) == "Only mol · 16 molecules"
    assert _enrich_title_with_system("Only beads", None, 8) == "Only beads · 8 beads"


def test_enrich_title_noop_when_empty() -> None:
    assert _enrich_title_with_system("Plain", None, None) == "Plain"


@pytest.mark.skipif(
    importlib.util.find_spec("matplotlib") is None,
    reason="matplotlib not installed",
)
def test_save_cv_plane_writes_png(tmp_path: Path) -> None:
    """RPMD order-parameter CV plane figure is written without error."""
    out = tmp_path / "cv_plane.png"
    q6 = np.array([0.1, 0.15, 0.12])
    qt = np.array([0.2, 0.22, 0.21])
    _save_cv_plane_figure(
        out,
        has_openmm=True,
        q6_o=q6,
        qt_o=qt,
        reference_data=None,
        cuda_sequential_data=None,
        has_ipi=False,
        q6_i=np.array([]),
        qt_i=np.array([]),
        tip4p_data=None,
        title="Test title",
        n_molecules=4,
        n_beads_ann=8,
        dpi=72,
    )
    assert out.is_file()
    assert out.stat().st_size > 100

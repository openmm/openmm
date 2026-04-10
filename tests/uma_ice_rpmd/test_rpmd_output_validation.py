from __future__ import annotations

from pathlib import Path

import pytest

from ipi_thermo_utils import build_thermo_lookup
from rpmd_output_validation import summarize_order_csv, validate_order_csv_complete


def test_summarize_order_csv_reads_first_and_last_time(tmp_path: Path) -> None:
    csv_path = tmp_path / "order.csv"
    csv_path.write_text(
        "step,time_ps,q6_mean,q_tet_mean\n"
        "0,0.000000,0.1,0.2\n"
        "10,0.002500,0.11,0.21\n"
        "20,0.005000,0.12,0.22\n",
        encoding="utf-8",
    )
    summary = summarize_order_csv(csv_path)
    assert summary.row_count == 3
    assert summary.first_time_ps == pytest.approx(0.0)
    assert summary.last_time_ps == pytest.approx(0.005)


def test_validate_order_csv_complete_rejects_truncated_file(tmp_path: Path) -> None:
    csv_path = tmp_path / "truncated.csv"
    csv_path.write_text(
        "step,time_ps,q6_mean,q_tet_mean\n"
        "0,0.000000,0.1,0.2\n"
        "10,0.002500,0.11,0.21\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="incomplete"):
        validate_order_csv_complete(
            csv_path,
            expected_final_time_ps=0.01,
            label="OpenMM stamped order CSV",
        )


def test_build_thermo_lookup_prefers_centroid_temperature() -> None:
    by_time, by_step = build_thermo_lookup(
        {
            "step": [0.0, 10.0],
            "time": [0.0, 0.0025],
            "temperature": [400.0, 410.0],
            "temperature(nm=0)": [240.0, 241.5],
            "potential": [-1.0, -2.0],
        }
    )
    assert by_step[0][0] == "240.0000"
    assert by_step[10][0] == "241.5000"
    assert by_time[0.0][0] == "240.0000"
    assert by_time[0.0025][1] != ""


def test_build_thermo_lookup_prefers_ring_when_requested() -> None:
    by_time, by_step = build_thermo_lookup(
        {
            "step": [0.0, 10.0],
            "time": [0.0, 0.0025],
            "temperature": [400.0, 410.0],
            "temperature(nm=0)": [240.0, 241.5],
            "potential": [-1.0, -2.0],
        },
        prefer_ring_temperature=True,
    )
    assert by_step[0][0] == "400.0000"
    assert by_step[10][0] == "410.0000"
    assert by_time[0.0025][0] == "410.0000"

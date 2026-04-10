"""Tests for trajectory Step parsing / verify helpers in run_ipi_lammps_uma_rpmd."""
from __future__ import annotations

from pathlib import Path

import pytest

# Import after path setup in conftest / same dir
from run_ipi_lammps_uma_rpmd import (
    IPI_TRAJ_OUTPUT_STRIDE,
    _ipi_order_from_traj_argv,
    _last_step_bead0_traj,
    _write_ipi_input,
    _verify_traj_matches_requested_steps,
)


def test_last_step_bead0_traj_picks_last_step(tmp_path: Path) -> None:
    ipi = tmp_path / "ipi"
    ipi.mkdir()
    traj = ipi / "ice__i-pi.traj_00.xyz"
    traj.write_text(
        "3\n# CELL Step: 0 Bead: 0\nO 0 0 0\nH 0 0 0\nH 0 0 0\n"
        "3\n# CELL Step: 10 Bead: 0\nO 1 0 0\nH 0 0 0\nH 0 0 0\n"
        "3\n# CELL Step: 10000 Bead: 0\nO 2 0 0\nH 0 0 0\nH 0 0 0\n",
        encoding="utf-8",
    )
    assert _last_step_bead0_traj(ipi) == 10000


def test_last_step_bead0_traj_missing_returns_none(tmp_path: Path) -> None:
    assert _last_step_bead0_traj(tmp_path) is None


def test_verify_warns_when_traj_short(capsys, tmp_path: Path) -> None:
    ipi = tmp_path / "ipi"
    ipi.mkdir()
    traj = ipi / "ice__i-pi.traj_00.xyz"
    traj.write_text(
        "3\n# Step: 310\nO 0 0 0\nH 0 0 0\nH 0 0 0\n",
        encoding="utf-8",
    )
    _verify_traj_matches_requested_steps(ipi, 10000, IPI_TRAJ_OUTPUT_STRIDE)
    out = capsys.readouterr()
    assert "last Step in bead-0 traj = 310" in out.out
    assert "WARNING" in out.err
    assert "well short" in out.err


def test_verify_fatal_when_traj_short(tmp_path: Path) -> None:
    ipi = tmp_path / "ipi"
    ipi.mkdir()
    traj = ipi / "ice__i-pi.traj_00.xyz"
    traj.write_text(
        "3\n# Step: 310\nO 0 0 0\nH 0 0 0\nH 0 0 0\n",
        encoding="utf-8",
    )
    with pytest.raises(SystemExit) as exc:
        _verify_traj_matches_requested_steps(
            ipi, 10000, IPI_TRAJ_OUTPUT_STRIDE, fatal_if_short=True
        )
    assert exc.value.code == 1


def test_ipi_order_argv_default_matches_openmm_path_integral(tmp_path: Path) -> None:
    """LAMMPS/i-PI orchestrator must not pass --centroid unless explicitly requested."""
    script_dir = tmp_path / "uma"
    script_dir.mkdir()
    (script_dir / "ipi_order_from_traj.py").write_text("# stub\n", encoding="utf-8")
    traj = script_dir / "ipi" / "ice__i-pi.traj_00.xyz"
    traj.parent.mkdir()
    traj.touch()
    out_csv = script_dir / "pipeline_out" / "ice_order_ipi_rpmd.csv"
    argv = _ipi_order_from_traj_argv(
        script_dir=script_dir,
        traj=traj,
        beads=32,
        dt_fs=0.1,
        order_csv=out_csv,
        centroid=False,
    )
    assert "--centroid" not in argv
    assert argv.count("--beads") == 1
    i = argv.index("--beads")
    assert argv[i + 1] == "32"


def test_ipi_order_argv_centroid_opt_in(tmp_path: Path) -> None:
    script_dir = tmp_path / "uma"
    script_dir.mkdir()
    (script_dir / "ipi_order_from_traj.py").write_text("# stub\n", encoding="utf-8")
    traj = script_dir / "ipi" / "t.xyz"
    traj.parent.mkdir()
    traj.touch()
    argv = _ipi_order_from_traj_argv(
        script_dir=script_dir,
        traj=traj,
        beads=4,
        dt_fs=0.5,
        order_csv=script_dir / "o.csv",
        centroid=True,
    )
    assert "--centroid" in argv


def test_verify_ok_when_traj_complete(capsys, tmp_path: Path) -> None:
    ipi = tmp_path / "ipi"
    ipi.mkdir()
    traj = ipi / "ice__i-pi.traj_00.xyz"
    traj.write_text(
        "3\n# Step: 10000\nO 0 0 0\nH 0 0 0\nH 0 0 0\n",
        encoding="utf-8",
    )
    _verify_traj_matches_requested_steps(ipi, 10000, IPI_TRAJ_OUTPUT_STRIDE)
    out = capsys.readouterr()
    assert "10000" in out.out
    assert "WARNING" not in out.err


def test_write_ipi_input_thermalizes_initial_velocities(tmp_path: Path) -> None:
    input_xml = tmp_path / "input.xml"
    _write_ipi_input(
        input_xml,
        molecules=64,
        beads=32,
        total_steps=400,
        dt_fs=0.25,
        seed=284759,
        thermostat_mode="pile_g",
        tau_fs=1000.0,
    )
    text = input_xml.read_text(encoding="utf-8")
    assert '<velocities mode="thermal" units="kelvin">243</velocities>' in text

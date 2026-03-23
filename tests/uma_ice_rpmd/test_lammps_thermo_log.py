"""Tests for LAMMPS thermo log parsing (T_K / PE_kj_mol merge into order CSV)."""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from lammps_thermo_log import EV_TO_KJ_PER_EV, parse_lammps_thermo_log


def test_parse_lammps_thermo_log_sample(tmp_path: Path) -> None:
    log = tmp_path / "run.log"
    log.write_text(
        textwrap.dedent(
            """\
            some preamble
               Step          Time           Temp          PotEng         KinEng         TotEng
                    75   0              243           -133146.24      5.9993511     -133140.25
                   100   0.0025         172.68095     -133144.68      4.2632661     -133140.41
            """
        ),
        encoding="utf-8",
    )
    d = parse_lammps_thermo_log(log)
    assert d[75][0] == pytest.approx(243.0)
    assert d[75][1] == pytest.approx(-133146.24 * EV_TO_KJ_PER_EV)
    assert d[100][0] == pytest.approx(172.68095)
    assert d[100][1] == pytest.approx(-133144.68 * EV_TO_KJ_PER_EV)


def test_parse_missing_file() -> None:
    assert parse_lammps_thermo_log(Path("/nonexistent/no.log")) == {}

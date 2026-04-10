#!/usr/bin/env python3
"""
Parse LAMMPS text logs for ``thermo_style`` tables that include **Step**, **Temp**, and **PotEng**.

Metal units: ``PotEng`` is total potential energy in **eV**; we convert to **kJ/mol** using the same
factor as ``compare_forces_openmm_vs_lammps.EV_TO_KJ`` so CSV ``PE_kj_mol`` matches OpenMM exports.
"""
from __future__ import annotations

from pathlib import Path

# eV (total system energy in LAMMPS metal) -> kJ/mol scale used elsewhere in this test suite
EV_TO_KJ_PER_EV = 96.4853


def parse_lammps_thermo_log(log_path: Path) -> dict[int, tuple[float, float]]:
    """
    Parameters
    ----------
    log_path
        LAMMPS ``-log`` file (must exist on disk).

    Returns
    -------
    dict[int, tuple[float, float]]
        Mapping **MD step** -> ``(T_K, PE_kj_mol)``. Only rows where the active thermo header
        includes ``PotEng`` are stored (typical MD block). Minimization blocks with different
        columns are skipped. Later duplicate steps overwrite earlier ones.
    """
    if not log_path.is_file():
        return {}
    lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    idx_step: int | None = None
    idx_temp: int | None = None
    idx_pe: int | None = None
    thermo: dict[int, tuple[float, float]] = {}

    for line in lines:
        ls = line.strip()
        if not ls:
            continue
        parts = ls.split()
        if len(parts) < 2:
            continue
        # New thermo header (may appear after minimize, again before MD, etc.)
        if parts[0] == "Step" and "Temp" in ls:
            try:
                idx_step = parts.index("Step")
                idx_temp = parts.index("Temp")
                idx_pe = parts.index("PotEng") if "PotEng" in parts else None
            except ValueError:
                idx_step = idx_temp = idx_pe = None
            continue
        if idx_step is None or idx_temp is None:
            continue
        try:
            step = int(parts[0])
        except ValueError:
            continue
        if idx_pe is None:
            continue
        try:
            T = float(parts[idx_temp])
            pe_ev = float(parts[idx_pe])
        except (IndexError, ValueError):
            continue
        thermo[step] = (T, pe_ev * EV_TO_KJ_PER_EV)

    return thermo

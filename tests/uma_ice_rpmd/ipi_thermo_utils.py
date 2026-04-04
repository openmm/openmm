"""Read i-PI ``<properties>`` files (e.g. ``ice__i-pi.md``) for post-processing.

``ipi_order_from_traj.py`` prefers **centroid kinetic ``T_K``** computed from
bead velocity trajectories (``vtraj``) using the same formula as OpenMM
(``centroid_kinetic_energy_and_temperature``). When ``vtraj`` files are absent,
``T_K`` falls back to this merge: ``temperature(nm=0)`` then plain
``temperature`` (unless ``prefer_ring_temperature=True`` swaps that order when
nm=0 is pathological). Plain ``temperature`` is i-PI’s ring kinetic
temperature and is **not** the same as centroid kinetic T.

``potential`` is converted from electronvolt (i-PI default in our template) to
kJ/mol using 1 eV × N_A × e / (kJ/mol) ≈ 96.4853 × (value in eV) for the
**total** extensive potential — comparable in magnitude to OpenMM's total ``PE_kj_mol``.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np

try:
    from ipi.utils.parsing import read_output as ipi_read_output

    _HAS_IPI = True
except ImportError:
    ipi_read_output = None
    _HAS_IPI = False

EV_TO_KJ_MOL = 96.485332123310018  # CODATA eV to kJ/mol conversion factor

_HEADER_RE = re.compile(
    r"#\s*(column|cols\.)\s+(\d+)(?:-(\d+))?\s*-->\s*([^\s\{\(]+)(?:\{([^\}]+)\})?(?:\(([^\)]+)\))?(?:\{([^\}]+)\})?\s*:\s*(.*)"
)


def read_ipi_properties_file(path: Path) -> tuple[dict[str, np.ndarray], dict[str, tuple]]:
    """Load i-PI property output; uses ``ipi.utils.parsing.read_output`` if installed."""
    if _HAS_IPI:
        return ipi_read_output(str(path))
    return _read_ipi_md_fallback(path)


def _read_ipi_md_fallback(path: Path) -> tuple[dict[str, np.ndarray], dict[str, tuple]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    header_lines = [ln for ln in lines if ln.startswith("#")]
    data_lines = [ln for ln in lines if not ln.startswith("#") and ln.strip()]
    properties: dict[str, dict] = {}
    for line in header_lines:
        m = _HEADER_RE.match(line)
        if not m:
            continue
        col_type, start_col, end_col, prop_name, u1, args, u2, desc = m.groups()
        col_info = f"{start_col}-{end_col}" if end_col else start_col
        units = u2 or u1
        if args is not None:
            prop_name = f"{prop_name}({args})"
        properties[col_info] = {"name": prop_name, "units": units or "", "description": desc or ""}

    values_dict: dict[str, list] = {properties[c]["name"]: [] for c in properties}
    info_dict = {
        properties[c]["name"]: (properties[c]["units"], properties[c]["description"])
        for c in properties
    }
    for line in data_lines:
        parts = line.split()
        for col_info, pinfo in properties.items():
            if "-" in col_info:
                a, b = map(int, col_info.split("-"))
                chunk = parts[a - 1 : b]
            else:
                chunk = [parts[int(col_info) - 1]]
            values_dict[pinfo["name"]].append([float(x) for x in chunk])
    for k, v in values_dict.items():
        values_dict[k] = np.array(v).squeeze()
    return values_dict, info_dict


def build_thermo_lookup(
    vals: dict[str, Any],
    *,
    prefer_ring_temperature: bool = False,
) -> tuple[dict[float, tuple[str, str]], dict[int, tuple[str, str]]]:
    """Build (time_ps_rounded -> (T_K, PE)), (ipi_step -> (T_K, PE)).

    By default, T_K uses ``temperature(nm=0)`` when present, else ``temperature``
    (centroid NM T, closest to OpenMM PILE-G ``T_K``).  With
    *prefer_ring_temperature*, use ``temperature`` first, then fall back to
    ``temperature(nm=0)`` (avoids bad nm=0 columns in some i-PI outputs).

    PE_kj_mol string from ``potential`` (eV) converted to kJ/mol when present.
    """
    n = 0
    for k in ("step", "time", "temperature", "temperature(nm=0)", "potential"):
        if k in vals:
            n = max(n, len(np.asarray(vals[k]).ravel()))
    if n == 0:
        return {}, {}

    def col(name: str) -> np.ndarray | None:
        if name not in vals:
            return None
        a = np.asarray(vals[name], dtype=np.float64).ravel()
        if a.size == n:
            return a
        if a.size == 1 and n > 1:
            return np.full(n, a[0])
        return None

    steps = col("step")
    times = col("time")
    t_nm0 = col("temperature(nm=0)")
    t_ext = col("temperature")
    pot_ev = col("potential")

    by_time: dict[float, tuple[str, str]] = {}
    by_step: dict[int, tuple[str, str]] = {}

    for i in range(n):
        tk_f: float | None = None
        if prefer_ring_temperature:
            if t_ext is not None and np.isfinite(t_ext[i]):
                tk_f = float(t_ext[i])
            elif t_nm0 is not None and np.isfinite(t_nm0[i]):
                tk_f = float(t_nm0[i])
        else:
            if t_nm0 is not None and np.isfinite(t_nm0[i]):
                tk_f = float(t_nm0[i])
            elif t_ext is not None and np.isfinite(t_ext[i]):
                tk_f = float(t_ext[i])
        tk_s = f"{tk_f:.4f}" if tk_f is not None and np.isfinite(tk_f) else ""

        pe_s = ""
        if pot_ev is not None and np.isfinite(pot_ev[i]):
            pe_s = f"{float(pot_ev[i]) * EV_TO_KJ_MOL:.6f}"

        pair = (tk_s, pe_s)
        if steps is not None:
            by_step[int(round(float(steps[i])))] = pair
        if times is not None:
            key_t = round(float(times[i]), 9)
            by_time[key_t] = pair

    return by_time, by_step


def lookup_thermo(
    *,
    time_ps: float,
    step: int,
    by_time: dict[float, tuple[str, str]],
    by_step: dict[int, tuple[str, str]],
) -> tuple[str, str]:
    """Return (T_K, PE_kj) strings for a trajectory frame."""
    key_t = round(float(time_ps), 9)
    if key_t in by_time:
        return by_time[key_t]
    if step in by_step:
        return by_step[step]
    for delta in (1, -1):
        if step + delta in by_step:
            return by_step[step + delta]
    return ("", "")


def load_thermo_lookups(
    md_path: Path,
    *,
    prefer_ring_temperature: bool = False,
) -> tuple[dict[float, tuple[str, str]], dict[int, tuple[str, str]]]:
    vals, _info = read_ipi_properties_file(md_path)
    return build_thermo_lookup(vals, prefer_ring_temperature=prefer_ring_temperature)

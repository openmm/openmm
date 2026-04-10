#!/usr/bin/env python3
"""
Plot i-PI ``<properties>`` thermo output (kinetic temperatures + potential).

Default input: ``ipi/ice__i-pi.md`` (from ``run_ipi_lammps_uma_rpmd.py`` template:
``temperature`` = extended ring-polymer MD temperature; ``temperature(nm=0)`` =
normal-mode 0 (centroid) kinetic temperature — closest analogue to OpenMM
``T_K`` in ``ice_order_openmm_rpmd.csv`` (centroid velocity temperature).

Uses ``ipi.utils.parsing.read_output`` when available; falls back to a light
header parser for the same file format.

Usage (from ``tests/uma_ice_rpmd``)::

  python plot_ipi_thermo.py
  python plot_ipi_thermo.py --md ipi/ice__i-pi.md -o pipeline_out/ipi_kinetic_temperature.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

from ipi_thermo_utils import read_ipi_properties_file


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot i-PI properties MD file (temperatures, potential)")
    ap.add_argument(
        "--md",
        type=Path,
        default=_SCRIPT_DIR / "ipi" / "ice__i-pi.md",
        help="i-PI properties output (default: ipi/ice__i-pi.md)",
    )
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=_SCRIPT_DIR / "pipeline_out" / "ipi_kinetic_temperature.png",
        help="Output PNG path",
    )
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    if not _HAS_MPL:
        print("matplotlib required: pip install matplotlib", file=sys.stderr)
        sys.exit(1)

    if not args.md.is_file():
        print(
            f"Missing i-PI properties file: {args.md}\n"
            "  Re-run i-PI after regenerating input.xml (includes <properties>), "
            "e.g. python run_ipi_lammps_uma_rpmd.py ...",
            file=sys.stderr,
        )
        sys.exit(1)

    vals, info = read_ipi_properties_file(args.md)
    time_ps = None
    if "time" in vals:
        time_ps = np.asarray(vals["time"], dtype=np.float64).ravel()
    elif "step" in vals:
        # Fallback: cannot know dt without input; user should include time in properties
        print("Warning: no 'time' column; using step index as x (add time{picosecond} to properties).", file=sys.stderr)
        time_ps = np.asarray(vals["step"], dtype=np.float64).ravel()

    t_ext = vals.get("temperature")
    t_nm0 = vals.get("temperature(nm=0)")
    pot = vals.get("potential")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(3, 1, figsize=(6.5, 8.0), sharex=True, constrained_layout=True)

    ax0, ax1, ax2 = axes
    if t_ext is not None:
        lab = "temperature (extended RPMD)" + (f" [{info['temperature'][0]}]" if "temperature" in info else "")
        ax0.plot(time_ps, np.asarray(t_ext).ravel(), "C0-", lw=1.2, label=lab)
    ax0.axhline(243.0, color="k", ls=":", lw=1.0, alpha=0.6, label="Bath 243 K")
    ax0.set_ylabel("T (K)")
    ax0.set_title("i-PI kinetic temperatures (properties file)")
    ax0.legend(loc="best", fontsize=9)
    ax0.grid(True, alpha=0.3)

    if t_nm0 is not None:
        lab = "temperature(nm=0) centroid mode" + (
            f" [{info['temperature(nm=0)'][0]}]" if "temperature(nm=0)" in info else ""
        )
        ax1.plot(time_ps, np.asarray(t_nm0).ravel(), "C3-", lw=1.2, label=lab)
    ax1.axhline(243.0, color="k", ls=":", lw=1.0, alpha=0.6, label="Bath 243 K")
    ax1.set_ylabel("T (K)")
    ax1.legend(loc="best", fontsize=9)
    ax1.grid(True, alpha=0.3)

    if pot is not None:
        p = np.asarray(pot).ravel()
        ax2.plot(time_ps, p, "C2-", lw=1.0, label="potential" + (f" [{info['potential'][0]}]" if "potential" in info else ""))
        ax2.set_ylabel(f"PE ({info.get('potential', ('eV',))[0]})")
    ax2.set_xlabel("time (ps)")
    ax2.legend(loc="best", fontsize=9)
    ax2.grid(True, alpha=0.3)

    fig.savefig(args.output, dpi=args.dpi)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()

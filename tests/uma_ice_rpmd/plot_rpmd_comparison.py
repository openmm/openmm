#!/usr/bin/env python3
"""
Plot Q6 order parameter: OpenMM RPMD vs i-PI + LAMMPS UMA.

Reads CSV files produced by:
  - OpenMM RPMD: test_uma_ice_rpmd.py or similar (ice_order_*.csv)
  - i-PI+LAMMPS: ipi_order_from_traj.py (ice_order_ipi_rpmd.csv)

Usage:
  cd tests/uma_ice_rpmd
  python plot_rpmd_comparison.py --openmm pipeline_out/ice_order_openmm_rpmd.csv --ipi pipeline_out/ice_order_ipi_rpmd.csv -o rpmd_comparison.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

_SCRIPT_DIR = Path(__file__).resolve().parent


def _read_csv(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (time_ps, q6_mean, q_tet_mean) from ice order CSV.
    Handles empty T_K/PE columns (e.g. from ipi_order_from_traj).
    """
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        return np.array([]), np.array([]), np.array([])
    header = lines[0]
    cols = [c.strip() for c in header.split(",")]
    idx_t = cols.index("time_ps")
    idx_q6 = cols.index("q6_mean")
    idx_qt = cols.index("q_tet_mean")
    t_list, q6_list, qt_list = [], [], []
    for line in lines[1:]:
        parts = [p.strip() for p in line.split(",")]
        if len(parts) <= max(idx_t, idx_q6, idx_qt):
            continue
        try:
            t_list.append(float(parts[idx_t]))
            q6_list.append(float(parts[idx_q6]))
            qt_list.append(float(parts[idx_qt]))
        except ValueError:
            continue
    return np.array(t_list), np.array(q6_list), np.array(qt_list)


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot OpenMM vs i-PI+LAMMPS Q6 comparison")
    ap.add_argument("--openmm", type=Path, default=_SCRIPT_DIR / "pipeline_out" / "ice_order_openmm_rpmd.csv")
    ap.add_argument("--ipi", type=Path, default=_SCRIPT_DIR / "pipeline_out" / "ice_order_ipi_rpmd.csv")
    ap.add_argument("-o", "--output", type=Path, default=_SCRIPT_DIR / "pipeline_out" / "rpmd_comparison.png")
    args = ap.parse_args()

    if not _HAS_MPL:
        print("matplotlib required: pip install matplotlib", file=sys.stderr)
        sys.exit(1)
    if not args.openmm.is_file():
        print(f"Missing OpenMM CSV: {args.openmm}", file=sys.stderr)
        sys.exit(1)
    if not args.ipi.is_file():
        print(f"Missing i-PI CSV: {args.ipi}", file=sys.stderr)
        sys.exit(1)

    t_o, q6_o, qt_o = _read_csv(args.openmm)
    t_i, q6_i, qt_i = _read_csv(args.ipi)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax1.plot(t_o, q6_o, "b-", label="OpenMM RPMD", alpha=0.8)
    ax1.plot(t_i, q6_i, "r--", label="i-PI + LAMMPS UMA", alpha=0.8)
    ax1.set_ylabel(r"$\langle Q_6 \rangle$")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.plot(t_o, qt_o, "b-", label="OpenMM RPMD", alpha=0.8)
    ax2.plot(t_i, qt_i, "r--", label="i-PI + LAMMPS UMA", alpha=0.8)
    ax2.set_ylabel(r"$\langle q_{\mathrm{tet}} \rangle$")
    ax2.set_xlabel("Time (ps)")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    fig.suptitle("RPMD ice order: OpenMM vs i-PI+LAMMPS (32 mol, 4 beads, 243 K, PILE-L)")
    plt.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()

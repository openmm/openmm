#!/usr/bin/env python3
"""
Plot Q6 / q_tet / kinetic T: OpenMM UMA RPMD vs i-PI + LAMMPS UMA, optionally TIP4P RPMD.

Third panel: ``T_K`` (centroid kinetic temperature) vs time, with optional bath reference
line (default 243 K for PILE-G). Use this to check thermostat behaviour.

Reads CSV files produced by:
  - OpenMM UMA RPMD: test_uma_ice_rpmd.py (ice_order_openmm_rpmd.csv)
  - TIP4P RPMD: run_openmm_tip4p_rpmd.py (ice_order_tip4p_rpmd.csv)
  - i-PI+LAMMPS: ipi_order_from_traj.py (ice_order_ipi_rpmd.csv)

Usage:
  cd tests/uma_ice_rpmd
  python plot_rpmd_comparison.py --openmm pipeline_out/ice_order_openmm_rpmd.csv --ipi pipeline_out/ice_order_ipi_rpmd.csv --tip4p pipeline_out/ice_order_tip4p_rpmd.csv -o rpmd_comparison.png
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


def _read_csv(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    """Return (time_ps, q6_mean, q_tet_mean, T_K) from ice order CSV.

    ``T_K`` is the instantaneous kinetic temperature (PI estimator) when the
    column exists; missing or empty cells become ``nan``. If the
    ``T_K`` column is absent, the fourth return is ``None`` (e.g. some i-PI exports).
    """
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        return np.array([]), np.array([]), np.array([]), None
    header = lines[0]
    cols = [c.strip() for c in header.split(",")]
    idx_t = cols.index("time_ps")
    idx_q6 = cols.index("q6_mean")
    idx_qt = cols.index("q_tet_mean")
    has_tk = "T_K" in cols
    idx_tk = cols.index("T_K") if has_tk else -1
    t_list, q6_list, qt_list = [], [], []
    tk_list: list[float] = []
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
        if has_tk:
            if idx_tk < len(parts) and parts[idx_tk] != "":
                try:
                    tk_list.append(float(parts[idx_tk]))
                except ValueError:
                    tk_list.append(float("nan"))
            else:
                tk_list.append(float("nan"))
    t_arr = np.array(t_list)
    q6_arr = np.array(q6_list)
    qt_arr = np.array(qt_list)
    tk_arr: np.ndarray | None = np.array(tk_list, dtype=np.float64) if has_tk else None
    if tk_arr is not None and tk_arr.shape != t_arr.shape:
        # Should not happen if every parsed row appends tk
        tk_arr = None
    return t_arr, q6_arr, qt_arr, tk_arr


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot OpenMM vs i-PI+LAMMPS Q6 comparison")
    ap.add_argument(
        "--openmm",
        type=Path,
        default=_SCRIPT_DIR / "pipeline_out" / "ice_order_openmm_rpmd.csv",
        help="OpenMM UMA RPMD order CSV (test_uma_ice_rpmd / run_openmm_rpmd_reference).",
    )
    ap.add_argument(
        "--ipi",
        type=Path,
        default=None,
        help=(
            "i-PI+LAMMPS order CSV; omit to skip i-PI curves "
            f"(default pipeline path: {_SCRIPT_DIR / 'pipeline_out' / 'ice_order_ipi_rpmd.csv'})."
        ),
    )
    ap.add_argument(
        "--tip4p",
        type=Path,
        default=None,
        help="TIP4P/2005f OpenMM RPMD order CSV (run_openmm_tip4p_rpmd.py)",
    )
    ap.add_argument("-o", "--output", type=Path, default=_SCRIPT_DIR / "pipeline_out" / "rpmd_comparison.png")
    ap.add_argument(
        "--title",
        type=str,
        default="RPMD ice order: OpenMM UMA vs i-PI+LAMMPS UMA (243 K, PILE-G)",
        help="Figure suptitle (include bead/molecule counts as needed).",
    )
    ap.add_argument(
        "--bath-temperature-k",
        type=float,
        default=243.0,
        help="Target bath temperature (K) for horizontal reference line on the T_K panel (PILE-G).",
    )
    ap.add_argument(
        "--no-temperature-line",
        action="store_true",
        help="Do not draw the bath-temperature reference line on the kinetic T panel.",
    )
    args = ap.parse_args()

    if not _HAS_MPL:
        print("matplotlib required: pip install matplotlib", file=sys.stderr)
        sys.exit(1)

    has_openmm = args.openmm.is_file()
    ipi_path = args.ipi
    has_ipi = ipi_path is not None and Path(ipi_path).is_file()
    tip4p_data: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None] | None = None
    if args.tip4p is not None and args.tip4p.is_file():
        tip4p_data = _read_csv(args.tip4p)

    if not has_openmm and not has_ipi and tip4p_data is None:
        print(
            "Need at least one of: existing --openmm CSV, existing --ipi CSV, or --tip4p CSV.",
            file=sys.stderr,
        )
        sys.exit(1)

    t_o = q6_o = qt_o = np.array([])
    tk_o: np.ndarray | None = None
    if has_openmm:
        t_o, q6_o, qt_o, tk_o = _read_csv(args.openmm)
    t_i = q6_i = qt_i = np.array([])
    tk_i: np.ndarray | None = None
    if has_ipi:
        t_i, q6_i, qt_i, tk_i = _read_csv(Path(ipi_path))

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9, 8.5), sharex=True)
    if has_openmm and t_o.size:
        ax1.plot(t_o, q6_o, "b-", label="OpenMM UMA RPMD", alpha=0.8)
        ax2.plot(t_o, qt_o, "b-", label="OpenMM UMA RPMD", alpha=0.8)
    if has_ipi and t_i.size:
        ax1.plot(t_i, q6_i, "r--", label="i-PI + LAMMPS UMA", alpha=0.8)
        ax2.plot(t_i, qt_i, "r--", label="i-PI + LAMMPS UMA", alpha=0.8)
    if tip4p_data is not None:
        t_t, q6_t, qt_t, tk_t = tip4p_data
        ax1.plot(t_t, q6_t, "g-.", label="TIP4P/2005f RPMD", alpha=0.85, lw=1.8)
        ax2.plot(t_t, qt_t, "g-.", label="TIP4P/2005f RPMD", alpha=0.85, lw=1.8)
    ax1.set_ylabel(r"$\langle Q_6 \rangle$")
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    ax2.set_ylabel(r"$\langle q_{\mathrm{tet}} \rangle$ (Chau–Harding, PI)")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    def _plot_tk(
        ax: "plt.Axes",
        t_arr: np.ndarray,
        tk_arr: np.ndarray | None,
        label: str,
        style: str,
        lw: float | None = None,
    ) -> None:
        if tk_arr is None or t_arr.size == 0:
            return
        mask = np.isfinite(tk_arr)
        if not np.any(mask):
            return
        ax.plot(t_arr[mask], tk_arr[mask], style, label=label, alpha=0.85, lw=lw)

    if has_openmm and t_o.size:
        _plot_tk(ax3, t_o, tk_o, "OpenMM UMA RPMD", "b-")
    if has_ipi and t_i.size:
        _plot_tk(ax3, t_i, tk_i, "i-PI + LAMMPS UMA", "r--")
    if tip4p_data is not None:
        t_t, q6_t, qt_t, tk_t = tip4p_data
        _plot_tk(ax3, t_t, tk_t, "TIP4P/2005f RPMD", "g-.", lw=1.8)
    if not args.no_temperature_line:
        ax3.axhline(
            args.bath_temperature_k,
            color="k",
            ls=":",
            lw=1.0,
            alpha=0.6,
            label=f"Bath target ({args.bath_temperature_k:.0f} K)",
        )
    ax3.set_ylabel(r"$T_{\mathrm{kin}}$ (K)")
    ax3.set_xlabel("Time (ps)")
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)

    fig.suptitle(args.title)
    plt.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()

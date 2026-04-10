#!/usr/bin/env python3
"""
Plot Steinhardt Q6 and tetrahedral order (q_tet) vs time from ice_order*.csv files.

Typical use (RPMD classical FF vs OpenMM UMA), from tests/uma_ice_rpmd:

  python plot_order_vs_time.py \\
    --classical pipeline_out/ice_order_classical.csv \\
    --uma-rpmd pipeline_out/ice_order_openmm_rpmd.csv \\
    -o pipeline_out/order_vs_time.png
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _load(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def _float(row: dict, key: str) -> float:
    s = row.get(key, "").strip()
    if not s:
        return float("nan")
    return float(s)


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot Q6 and q_tet vs time from order CSVs")
    ap.add_argument(
        "--classical",
        type=Path,
        default=None,
        help="CSV from classical FF RPMD (e.g. ice_order_classical.csv)",
    )
    ap.add_argument(
        "--uma-rpmd",
        type=Path,
        default=None,
        help="CSV from OpenMM UMA RPMD (e.g. ice_order_openmm_rpmd.csv)",
    )
    ap.add_argument(
        "--uma-md",
        type=Path,
        default=None,
        help="Optional: classical MD UMA CSV for a third curve",
    )
    ap.add_argument("-o", "--output", type=Path, default=Path("order_vs_time.png"))
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, constrained_layout=True)

    def plot_series(ax_q6, ax_qt, path: Path, label: str, color: str, linestyle: str = "-") -> None:
        rows = _load(path)
        t = [_float(r, "time_ps") for r in rows]
        q6 = [_float(r, "q6_mean") for r in rows]
        qt = [_float(r, "q_tet_mean") for r in rows]
        ax_q6.plot(t, q6, color=color, ls=linestyle, lw=1.5, label=label)
        ax_qt.plot(t, qt, color=color, ls=linestyle, lw=1.5, label=label)

    if args.classical and args.classical.is_file():
        plot_series(axes[0], axes[1], args.classical, "Classical FF (RPMD)", "#1f77b4", "-")
    if args.uma_rpmd and args.uma_rpmd.is_file():
        plot_series(axes[0], axes[1], args.uma_rpmd, "OpenMM UMA (RPMD)", "#d62728", "-")
    if args.uma_md and args.uma_md.is_file():
        plot_series(axes[0], axes[1], args.uma_md, "OpenMM UMA (MD)", "#2ca02c", "--")

    axes[0].set_ylabel(r"Mean $Q_6$ (Steinhardt)")
    axes[0].set_title("Ice order parameters vs time")
    axes[0].legend(loc="best", fontsize=9)
    axes[0].grid(True, alpha=0.3)

    axes[1].set_ylabel(r"Mean $q_{\mathrm{tet}}$")
    axes[1].set_xlabel("Time (ps)")
    axes[1].legend(loc="best", fontsize=9)
    axes[1].grid(True, alpha=0.3)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi)
    print(f"Wrote {args.output.resolve()}")


if __name__ == "__main__":
    main()

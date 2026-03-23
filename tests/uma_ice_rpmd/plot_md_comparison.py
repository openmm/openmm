#!/usr/bin/env python3
"""
Plot **classical (single-bead) MD** ice order parameters only — no RPMD.

Three curves, same layout as ``plot_rpmd_comparison.py`` (⟨Q6⟩, ⟨q_tet⟩, T_kin vs time):

  - **LAMMPS UMA (MD)** — ``ice_order_uma_lammps.csv`` (``run_ice_pipeline`` / LAMMPS dump → order)
  - **OpenMM UMA (MD)** — ``ice_order_uma_openmm.csv`` (``run_openmm_ice_lammps_match``)
  - **TIP4P/2005f (MD)** — ``ice_order_classical.csv`` (``run_openmm_ice_classical_flex``)

This script does **not** read RPMD exports (``*_rpmd.csv``). Use ``plot_rpmd_comparison.py`` for path-integral runs.

Usage:
  cd tests/uma_ice_rpmd
  python plot_md_comparison.py -o pipeline_out/md_order_comparison.png
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
_OUT = _SCRIPT_DIR / "pipeline_out"

# Reuse CSV parsing and title helpers from the RPMD plotter (same column schema).
from plot_rpmd_comparison import (  # noqa: E402
    _DEFAULT_DPI,
    _FIGSIZE_IN,
    _PLOT_RC,
    _enrich_title_with_system,
    _format_title_two_lines,
    _read_csv,
    _read_n_oxygen_from_csv,
)

# Distinct colors — all solid MD trajectories (no RPMD curve styles here).
_STYLES: dict[str, tuple[str, str]] = {
    "lammps": ("#ff7f0e", "LAMMPS UMA (MD)"),
    "openmm": ("#1f77b4", "OpenMM UMA (MD)"),
    "tip4p": ("#2ca02c", "TIP4P/2005f (MD)"),
}


def _plot_tk(
    ax: "plt.Axes",
    t_arr: np.ndarray,
    tk_arr: np.ndarray | None,
    label: str,
    color: str,
    lw: float = 1.8,
) -> None:
    if tk_arr is None or t_arr.size == 0:
        return
    mask = np.isfinite(tk_arr)
    if not np.any(mask):
        return
    ax.plot(t_arr[mask], tk_arr[mask], "-", color=color, label=label, alpha=0.9, lw=lw)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Plot LAMMPS vs OpenMM vs TIP4P classical MD ice order (Q6, q_tet, T_kin)"
    )
    ap.add_argument(
        "--lammps",
        type=Path,
        default=_OUT / "ice_order_uma_lammps.csv",
        help="LAMMPS UMA MD order CSV (default: pipeline_out/ice_order_uma_lammps.csv).",
    )
    ap.add_argument(
        "--openmm-md",
        type=Path,
        dest="openmm_md",
        default=_OUT / "ice_order_uma_openmm.csv",
        help="OpenMM UMA MD order CSV (default: pipeline_out/ice_order_uma_openmm.csv).",
    )
    ap.add_argument(
        "--tip4p-md",
        type=Path,
        dest="tip4p_md",
        default=_OUT / "ice_order_classical.csv",
        help="TIP4P/2005f classical MD order CSV (default: pipeline_out/ice_order_classical.csv).",
    )
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=_OUT / "md_order_comparison.png",
        help="Output PNG path.",
    )
    ap.add_argument(
        "--title",
        type=str,
        default=(
            "Ice order (classical MD): LAMMPS UMA vs OpenMM UMA vs TIP4P/2005f\n"
            "(243 K, Langevin)"
        ),
        help="Figure suptitle; use \\n or split at em dash ' — '.",
    )
    ap.add_argument(
        "--bath-temperature-k",
        type=float,
        default=243.0,
        help="Horizontal reference on the T_kin panel (default 243 K).",
    )
    ap.add_argument(
        "--no-temperature-line",
        action="store_true",
        help="Do not draw the bath-temperature reference line.",
    )
    ap.add_argument("--dpi", type=int, default=_DEFAULT_DPI, help=f"PNG DPI (default {_DEFAULT_DPI}).")
    ap.add_argument(
        "--molecules",
        type=int,
        default=None,
        metavar="N",
        help="Molecule count in suptitle; default: infer n_oxygen from CSVs when present.",
    )
    args = ap.parse_args()

    if not _HAS_MPL:
        print("matplotlib required: pip install matplotlib", file=sys.stderr)
        sys.exit(1)

    paths = {
        "lammps": args.lammps,
        "openmm": args.openmm_md,
        "tip4p": args.tip4p_md,
    }
    data: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]] = {}
    for key, p in paths.items():
        if p.is_file():
            data[key] = _read_csv(p)
        else:
            print(f"  Skip (missing): {p}", file=sys.stderr)

    if not data:
        print(
            "Need at least one existing CSV among --lammps, --openmm-md, --tip4p-md.",
            file=sys.stderr,
        )
        sys.exit(1)

    with plt.rc_context(_PLOT_RC):
        fig, (ax1, ax2, ax3) = plt.subplots(
            3,
            1,
            figsize=_FIGSIZE_IN,
            sharex=True,
            constrained_layout=True,
        )
        for key in ("lammps", "openmm", "tip4p"):
            if key not in data:
                continue
            color, label = _STYLES[key]
            t_arr, q6, qt, tk = data[key]
            if t_arr.size == 0:
                continue
            ax1.plot(t_arr, q6, "-", color=color, label=label, alpha=0.9, lw=1.8)
            ax2.plot(t_arr, qt, "-", color=color, label=label, alpha=0.9, lw=1.8)
            _plot_tk(ax3, t_arr, tk, label, color)

        ax1.set_ylabel(r"$\langle Q_6 \rangle$")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        ax2.set_ylabel(r"$\langle q_{\mathrm{tet}} \rangle$")
        ax2.legend()
        ax2.grid(True, alpha=0.3)

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
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        title = _format_title_two_lines(args.title)
        n_mol = args.molecules
        if n_mol is None:
            for p in (args.lammps, args.openmm_md, args.tip4p_md):
                n_mol = _read_n_oxygen_from_csv(p)
                if n_mol is not None:
                    break
        # Classical MD: annotate molecules only (no bead count).
        fig.suptitle(_enrich_title_with_system(title, n_mol, None))
        fig.set_constrained_layout_pads(w_pad=0.02, h_pad=0.0, hspace=0.06)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(args.output, dpi=args.dpi, bbox_inches="tight", pad_inches=0.08)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Plot per-particle force component errors from debug_force_components.py.

Reads force_component_errors_per_particle.npz (or --input) and produces:
  1. err_bond, err_coulomb, err_lj, err_cavity vs particle_index (line or scatter).
  2. Distribution (histogram or CDF) per component so you can see which contributes most.

Usage:
  python plot_force_component_errors.py
  python plot_force_component_errors.py --input force_component_errors_per_particle.npz --out errors.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

_SCRIPT_DIR = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot per-particle force component errors (err vs particle index, and distributions)."
    )
    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help="NPZ or CSV from debug_force_components.py (default: force_component_errors_per_particle.npz in script dir)",
    )
    parser.add_argument("--out", type=str, default="force_component_errors.png", help="Output plot file")
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args()

    inp = Path(args.input) if args.input else _SCRIPT_DIR / "force_component_errors_per_particle.npz"
    if not inp.exists():
        csv_path = inp.with_suffix(".csv") if inp.suffix == ".npz" else _SCRIPT_DIR / "force_component_errors_per_particle.csv"
        if csv_path.exists():
            # Load from CSV
            data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
            particle_index = data[:, 0].astype(np.int32)
            err_bond = data[:, 1]
            err_coulomb = data[:, 2]
            err_lj = data[:, 3]
            err_cavity = data[:, 4]
        else:
            raise FileNotFoundError(
                f"Neither {inp} nor {csv_path} found. Run debug_force_components.py first."
            )
    else:
        if inp.suffix == ".npz":
            with np.load(inp, allow_pickle=False) as z:
                particle_index = np.asarray(z["particle_index"], dtype=np.int32)
                err_bond = np.asarray(z["err_bond"], dtype=np.float64)
                err_coulomb = np.asarray(z["err_coulomb"], dtype=np.float64)
                err_lj = np.asarray(z["err_lj"], dtype=np.float64)
                err_cavity = np.asarray(z["err_cavity"], dtype=np.float64)
        else:
            data = np.loadtxt(inp, delimiter=",", skiprows=1)
            particle_index = data[:, 0].astype(np.int32)
            err_bond, err_coulomb, err_lj, err_cavity = data[:, 1], data[:, 2], data[:, 3], data[:, 4]

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib required: pip install matplotlib", flush=True)
        return

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=False)

    # 1) err vs particle index (log scale for y so bond doesn't swamp others)
    ax = axes[0]
    x = np.arange(len(particle_index), dtype=np.int32) if particle_index.size and np.all(particle_index == np.arange(len(particle_index))) else particle_index
    for lab, arr, color in [
        ("bond", err_bond, "C0"),
        ("Coulomb", err_coulomb, "C1"),
        ("LJ", err_lj, "C2"),
        ("Cavity", err_cavity, "C3"),
    ]:
        valid = np.isfinite(arr)
        if np.any(valid):
            ax.scatter(x[valid], arr[valid], label=lab, s=4, alpha=0.6, c=color)
    ax.set_yscale("log")
    ax.set_xlabel("Particle index (OpenMM order)")
    ax.set_ylabel("|OpenMM − HOOMD| L2 (Ha/Bohr)")
    ax.set_title("Per-particle force component error")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    # 2) Distribution: CDF or histogram per component
    ax = axes[1]
    for lab, arr, color in [
        ("bond", err_bond, "C0"),
        ("Coulomb", err_coulomb, "C1"),
        ("LJ", err_lj, "C2"),
        ("Cavity", err_cavity, "C3"),
    ]:
        a = np.asarray(arr, dtype=np.float64)
        a = a[np.isfinite(a)]
        if a.size == 0:
            continue
        a = np.sort(a)
        cdf = np.arange(1, a.size + 1, dtype=np.float64) / a.size
        ax.plot(a, cdf, label=lab, color=color)
    ax.set_xscale("log")
    ax.set_xlabel("Error (Ha/Bohr)")
    ax.set_ylabel("CDF")
    ax.set_title("Distribution of per-particle error by component")
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_path = Path(args.out)
    fig.savefig(out_path, dpi=args.dpi, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()

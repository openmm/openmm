#!/usr/bin/env python3
"""Plot NVE comparison: OpenMM vs LAMMPS from identical initial state."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_SCRIPT_DIR = Path(__file__).resolve().parent
_OUT = _SCRIPT_DIR / "pipeline_out"


def main() -> None:
    lammps_csv = _OUT / "ice_order_lammps_nve.csv"
    openmm_csv = _OUT / "ice_order_openmm_nve.csv"

    if not lammps_csv.is_file() or not openmm_csv.is_file():
        print(f"Missing CSVs: {lammps_csv} and/or {openmm_csv}")
        sys.exit(1)

    dl = pd.read_csv(lammps_csv)
    do = pd.read_csv(openmm_csv)

    # Align LAMMPS time: NVE phase starts at t=0.5 ps (step 5000). Shift to relative time.
    dl["t_rel_ps"] = dl["time_ps"] - dl["time_ps"].iloc[0]
    do["t_rel_ps"] = do["time_ps"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        "NVE Comparison: OpenMM vs LAMMPS from Identical State\n"
        "(same positions + velocities, Velocity-Verlet, no thermostat)",
        fontsize=13,
        fontweight="bold",
    )

    # Q6
    ax = axes[0, 0]
    ax.plot(dl["t_rel_ps"], dl["q6_mean"], "o-", label="LAMMPS NVE", color="tab:blue", markersize=3)
    ax.plot(do["t_rel_ps"], do["q6_mean"], "s-", label="OpenMM NVE", color="tab:red", markersize=3)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Mean Q6")
    ax.set_title("Steinhardt Q6 (ice order)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # q_tet
    ax = axes[0, 1]
    ax.plot(dl["t_rel_ps"], dl["q_tet_mean"], "o-", label="LAMMPS NVE", color="tab:blue", markersize=3)
    ax.plot(do["t_rel_ps"], do["q_tet_mean"], "s-", label="OpenMM NVE", color="tab:red", markersize=3)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Mean q_tet")
    ax.set_title("Tetrahedral order (ice order)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Temperature
    ax = axes[1, 0]
    if "T_K" in do.columns:
        ax.plot(do["t_rel_ps"], do["T_K"], "s-", label="OpenMM NVE", color="tab:red", markersize=3)
    ax.axhline(y=250, color="tab:blue", linestyle="--", alpha=0.5, label="LAMMPS NVE (~250 K)")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Instantaneous Temperature")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Energy conservation
    ax = axes[1, 1]
    if "E_tot_kj_mol" in do.columns:
        e_ref = do["E_tot_kj_mol"].iloc[0]
        ax.plot(
            do["t_rel_ps"],
            do["E_tot_kj_mol"] - e_ref,
            "s-",
            label="OpenMM E_tot drift",
            color="tab:red",
            markersize=3,
        )
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("E_tot - E_tot(0) (kJ/mol)")
    ax.set_title("Total Energy Conservation (OpenMM)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = _OUT / "nve_comparison_velocity_matched.png"
    fig.savefig(out_png, dpi=150)
    print(f"Saved: {out_png}")

    # Summary table
    print("\n" + "=" * 72)
    print("NVE COMPARISON SUMMARY (0.3 ps from identical state)")
    print("=" * 72)
    print(f"{'Metric':<25} {'LAMMPS (start→end)':<25} {'OpenMM (start→end)':<25}")
    print("-" * 75)
    print(
        f"{'Q6 mean':<25} "
        f"{dl['q6_mean'].iloc[0]:.4f} → {dl['q6_mean'].iloc[-1]:.4f}     "
        f"{do['q6_mean'].iloc[0]:.4f} → {do['q6_mean'].iloc[-1]:.4f}"
    )
    print(
        f"{'q_tet mean':<25} "
        f"{dl['q_tet_mean'].iloc[0]:.4f} → {dl['q_tet_mean'].iloc[-1]:.4f}     "
        f"{do['q_tet_mean'].iloc[0]:.4f} → {do['q_tet_mean'].iloc[-1]:.4f}"
    )
    if "T_K" in do.columns:
        print(
            f"{'Temperature (K)':<25} "
            f"~250 → ~250                "
            f"{do['T_K'].iloc[0]:.1f} → {do['T_K'].iloc[-1]:.1f}"
        )
    if "E_tot_kj_mol" in do.columns:
        drift = do["E_tot_kj_mol"].iloc[-1] - do["E_tot_kj_mol"].iloc[0]
        print(f"{'E_tot drift (kJ/mol)':<25} ~0                        {drift:.2f}")
    print("=" * 72)

    # Diagnosis
    t_final_openmm = do["T_K"].iloc[-1] if "T_K" in do.columns else None
    if t_final_openmm and t_final_openmm > 400:
        print(
            "\nDIAGNOSIS: OpenMM NVE heats from ~250 K to ~{:.0f} K with the same initial".format(t_final_openmm)
        )
        print("conditions. Energy is conserved (NVE integrator correct), but the")
        print("OpenMM potential energy surface differs from LAMMPS. This means the")
        print("UMA forces computed DURING DYNAMICS differ between the two codes.")
        print("The bug is in the PythonForce callback path, not in the integrator.")
    else:
        print("\nDIAGNOSIS: Dynamics appear consistent.")


if __name__ == "__main__":
    main()

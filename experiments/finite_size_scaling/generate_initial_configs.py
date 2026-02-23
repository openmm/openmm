#!/usr/bin/env python3
"""
Generate initial GSD configurations for the finite-size scaling campaign.

For each system size N:
  1. Create a lattice config with equilibrium bond lengths via initlattice_equilibrium.py
  2. Print the command to run a long no-cavity equilibration on the cluster

The equilibration run (150 ns at T=100 K, no cavity) must be done separately
because it is too expensive to run inline. After equilibration, extract N_T
independent frames separated by 300 ps for production initial conditions.
"""

import argparse
import subprocess
import sys
from pathlib import Path

from config import (
    SYSTEM_SIZES, CONFIGS_DIR, TEMPERATURE_K, CAVITY_FREQ_CM,
    INITLATTICE_SCRIPT, RUN_SIMULATION_SCRIPT,
    get_campaign_table, EQUIL_SEPARATION_PS, EQUIL_TOTAL_NS,
)


def generate_lattice(n_molecules: int, seed: int = 42):
    """Generate a lattice GSD for the given system size."""
    CONFIGS_DIR.mkdir(parents=True, exist_ok=True)
    out_dir = CONFIGS_DIR / f"N{n_molecules}"
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable, str(INITLATTICE_SCRIPT),
        "--job-dir", str(out_dir),
        "--replica", "0",
        "--nmol", str(n_molecules),
        "--temperature", str(TEMPERATURE_K),
        "--seed", str(seed),
    ]
    print(f"[N={n_molecules}] Generating lattice config ...")
    print(f"  cmd: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr}")
        return False
    print(f"  OK: {out_dir / 'molecular-0.gsd'}")
    return True


def print_equilibration_commands(sizes=None):
    """Print commands needed to equilibrate each system size."""
    table = get_campaign_table()
    if sizes:
        table = [row for row in table if row["N"] in sizes]

    print("\n" + "=" * 70)
    print("EQUILIBRATION COMMANDS")
    print("Run these (no cavity, long NVT) to generate independent configs.")
    print("After completion, the GSD trajectory contains frames separated")
    print(f"by {EQUIL_SEPARATION_PS} ps that serve as initial conditions.")
    print("=" * 70)

    for row in table:
        n = row["N"]
        n_replicas = row["replicas"]
        total_ps = EQUIL_TOTAL_NS * 1000
        lattice_gsd = CONFIGS_DIR / f"N{n}" / "molecular-0.gsd"
        equil_gsd = CONFIGS_DIR / f"N{n}" / "equil.gsd"

        print(f"\n# N = {n} ({n_replicas} replicas needed)")
        print(f"# Requires {EQUIL_TOTAL_NS} ns equilibration, extracting {n_replicas} frames")
        min_runtime_ps = n_replicas * EQUIL_SEPARATION_PS
        actual_runtime_ps = max(total_ps, min_runtime_ps)
        print(f"python {RUN_SIMULATION_SCRIPT} \\")
        print(f"  --dimers {n} --constant-density \\")
        print(f"  --no-cavity --temp {TEMPERATURE_K} \\")
        print(f"  --dt 0.001 --equil 0 --prod {actual_runtime_ps:.0f} \\")
        print(f"  --no-dipole --no-minimize \\")
        print(f"  --bussi-tau 1.0 \\")
        print(f"  --initial-gsd {lattice_gsd} --initial-gsd-frame 0 \\")
        print(f"  --seed 42")
        print(f"# Then extract {n_replicas} frames separated by {EQUIL_SEPARATION_PS} ps")


def main():
    parser = argparse.ArgumentParser(description="Generate initial configs for finite-size scaling")
    parser.add_argument("--sizes", type=int, nargs="*", default=None,
                        help="System sizes to generate (default: all)")
    parser.add_argument("--lattice-only", action="store_true",
                        help="Only generate lattice configs, skip equilibration commands")
    args = parser.parse_args()

    sizes = args.sizes or SYSTEM_SIZES

    for n in sizes:
        generate_lattice(n)

    if not args.lattice_only:
        print_equilibration_commands(sizes)


if __name__ == "__main__":
    main()

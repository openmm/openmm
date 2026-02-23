#!/usr/bin/env python3
"""
Run a single production simulation for one (N, replica_id) pair.

This is the entry point called by both run_local.sh and SLURM array jobs.
It reads all parameters from config.py and invokes run_simulation.py with
the correct arguments.

Usage:
  python run_single.py --N 250 --replica 0
  python run_single.py --N 250 --replica 0 --no-cavity   # baseline run (lambda=0)
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

from config import (
    TEMPERATURE_K, CAVITY_FREQ_CM, DT_PS, SWITCH_TIME_PS, PRODUCTION_PS,
    BUSSI_TAU_PS, FKT_KMAG, FKT_NUM_WAVEVECTORS, FKT_REF_INTERVAL_PS,
    FKT_MAX_REFS, FKT_OUTPUT_PERIOD_PS, RUN_SIMULATION_SCRIPT,
    CONFIGS_DIR, RESULTS_DIR,
    lambda_scaled, num_replicas,
)


def run_single(n_molecules: int, replica_id: int, baseline: bool = False,
               extra_args: list = None):
    """Launch run_simulation.py for one (N, replica) pair."""

    lam = 0.0 if baseline else lambda_scaled(n_molecules)
    tag = "baseline" if baseline else "cavity"
    out_dir = RESULTS_DIR / f"N{n_molecules}" / tag / f"replica_{replica_id:04d}"
    out_dir.mkdir(parents=True, exist_ok=True)

    gsd_path = CONFIGS_DIR / f"N{n_molecules}" / "molecular-0.gsd"
    if not gsd_path.exists():
        print(f"ERROR: initial GSD not found: {gsd_path}")
        print("Run generate_initial_configs.py first.")
        sys.exit(1)

    seed = 42 + replica_id

    cmd = [
        sys.executable, str(RUN_SIMULATION_SCRIPT),
        "--dimers", str(n_molecules),
        "--constant-density",
        "--lambda", f"{lam:.6f}",
        "--cavity-freq", str(CAVITY_FREQ_CM),
        "--temp", str(TEMPERATURE_K),
        "--dt", str(DT_PS),
        "--equil", "0",
        "--prod", str(PRODUCTION_PS),
        "--bussi-tau", str(BUSSI_TAU_PS),
        "--no-dipole",
        "--no-minimize",
        "--initial-gsd", str(gsd_path),
        "--initial-gsd-frame", str(replica_id),
        "--seed", str(seed),
        "--enable-fkt",
        "--fkt-kmag", str(FKT_KMAG),
        "--fkt-num-wavevectors", str(FKT_NUM_WAVEVECTORS),
        "--fkt-ref-interval-ps", str(FKT_REF_INTERVAL_PS),
        "--fkt-max-refs", str(FKT_MAX_REFS),
        "--fkt-output-period-ps", str(FKT_OUTPUT_PERIOD_PS),
        "--report-interval", "10000",
    ]

    if baseline:
        cmd.append("--no-cavity")
    else:
        cmd.extend(["--switch-time", str(SWITCH_TIME_PS)])

    if extra_args:
        cmd.extend(extra_args)

    print(f"[N={n_molecules}, replica={replica_id}, {tag}]")
    print(f"  lambda = {lam:.6f}")
    print(f"  output dir: {out_dir}")
    print(f"  cmd: {' '.join(cmd[:8])} ...")

    result = subprocess.run(cmd, cwd=str(out_dir))
    return result.returncode


def main():
    parser = argparse.ArgumentParser(description="Run single finite-size scaling simulation")
    parser.add_argument("--N", type=int, required=True, help="Number of molecules")
    parser.add_argument("--replica", type=int, required=True, help="Replica index (0-based)")
    parser.add_argument("--no-cavity", action="store_true", dest="baseline",
                        help="Baseline run with lambda=0 (no cavity)")
    args, extra = parser.parse_known_args()

    rc = run_single(args.N, args.replica, baseline=args.baseline, extra_args=extra or None)
    sys.exit(rc)


if __name__ == "__main__":
    main()

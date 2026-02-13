#!/usr/bin/env python3
"""
Run OpenMM and cav-hoomd (HOOMD) cavity simulations side by side and average F(k,t) on the fly.

Runs run_simulation.py (OpenMM) and run_cav_hoomd_advanced.py (cav-hoomd) in two
independent threads so OpenMM does not wait for HOOMD (or vice versa). Each backend
runs N replicas, parses its F(k,t) files, updates running mean/variance via Welford's
algorithm, and writes its own averaged F(k,t) files (e.g. fkt_averaged_openmm_* and
fkt_averaged_hoomd_*).

By default both backends use the same initial conditions: --initial-gsd init-0.gsd and
--initial-gsd-frame 0 (one shared configuration for all replicas). Use --initial-gsd-frame -1
to use frame i for replica i instead.

Prerequisite for cav-hoomd: ensure init-0.gsd exists (500 molecules + cavity particle L;
same O-O/N-N composition as current runs; e.g. from a prior cavity run or create molecular-0.gsd then add cavity).

Usage (from examples/cavity/dimer_system/):
  python run_fkt_replica_average.py --replicas 500
  python run_fkt_replica_average.py --replicas 100 --no-hoomd
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Fixed run_simulation args: bare glassy system (no cavity), matches cav-hoomd --no-cavity
DEFAULT_RUN_SIMULATION_ARGS = [
    "--dimers", "250",
    "--g", "0.0",
    "--temp", "100.0",
    "--dt", "0.001",
    "--equil", "200",
    "--prod", "300",
    "--cavity-freq", "1555",
    "--report-interval", "10000",
    "--enable-fkt",
    "--fkt-output-period-ps", "0.1",
    "--fkt-ref-interval-ps", "200.0",
    "--fkt-max-refs", "5",
    "--constant-density",
    "--no-dipole",
    "--no-minimize",
    "--no-cavity",
]

# Cav-hoomd args: match OpenMM (runtime 200+300=500 ps, same F(k,t) params)
DEFAULT_RUN_CAV_HOOMD_ARGS = [
    "--no-cavity",
    "--temperature", "100",
    "--runtime", "500",
    "--fixed-timestep", "--timestep", "1.0",
    "--molecular-bath", "bussi", "--molecular-tau", "5.0",
    "--enable-fkt", "--fkt-kmag", "113.4", "--fkt-wavevectors", "50",
    "--fkt-ref-interval", "200.0", "--fkt-max-refs", "5",
    "--fkt-output-period-ps", "0.1",
    "--input-gsd", "init-0.gsd",
    "--device", "GPU",
    "--console-output-period-ps", "10",
]


def _lag_key(lag_ps: float, decimals: int = 3) -> float:
    """Round lag time for stable dict key across replicas."""
    return round(lag_ps, decimals)


def parse_fkt_file(path: Path) -> List[Tuple[float, float]]:
    """
    Parse one F(k,t) file: skip # lines, return list of (lag_time_ps, fkt_value).
    """
    out: List[Tuple[float, float]] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            try:
                lag_ps = float(parts[0])
                fkt = float(parts[1])
                out.append((lag_ps, fkt))
            except ValueError:
                continue
    return out


def parse_replica_fkt_files(
    directory: Path,
    replica_prefix: str,
    replica_id: int,
) -> Dict[int, List[Tuple[float, float]]]:
    """
    Find all {replica_prefix}_{id:04d}_fkt_ref_{ref_idx:03d}.txt and return
    ref_idx -> [(lag_ps, fkt), ...].
    """
    pattern = re.compile(
        re.escape(f"{replica_prefix}_{replica_id:04d}_fkt_ref_")
        + r"(\d+)\.txt"
    )
    result: Dict[int, List[Tuple[float, float]]] = {}
    for p in directory.iterdir():
        if not p.is_file():
            continue
        m = pattern.fullmatch(p.name)
        if m is None:
            continue
        ref_idx = int(m.group(1))
        result[ref_idx] = parse_fkt_file(p)
    return result


def parse_hoomd_replica_fkt_files(
    script_dir: Path,
    replica_id: int,
) -> Dict[int, List[Tuple[float, float]]]:
    """
    Find all no_cavity/prod-{replica_id}_fkt_ref_{ref_idx}.txt and return
    ref_idx -> [(lag_ps, fkt), ...].
    """
    no_cavity_dir = script_dir / "no_cavity"
    if not no_cavity_dir.is_dir():
        return {}
    pattern = re.compile(
        re.escape(f"prod-{replica_id}_fkt_ref_") + r"(\d+)\.txt"
    )
    result: Dict[int, List[Tuple[float, float]]] = {}
    for p in no_cavity_dir.iterdir():
        if not p.is_file():
            continue
        m = pattern.fullmatch(p.name)
        if m is None:
            continue
        ref_idx = int(m.group(1))
        result[ref_idx] = parse_fkt_file(p)
    return result


def run_single_hoomd_replica(
    script_dir: Path,
    replica_id: int,
    base_seed: int,
    initial_gsd: Optional[str] = None,
    initial_gsd_frame: int = 0,
) -> bool:
    """Run run_cav_hoomd_advanced.py once for a single replica; return True on success.
    Uses the same initial GSD and frame as OpenMM: initial_gsd overrides --input-gsd when set;
    if initial_gsd_frame >= 0 that frame is used for all replicas; if < 0, use replica_id.
    """
    frame = replica_id if initial_gsd_frame < 0 else initial_gsd_frame
    cav_args = list(DEFAULT_RUN_CAV_HOOMD_ARGS)
    if initial_gsd:
        for j in range(len(cav_args) - 1):
            if cav_args[j] == "--input-gsd":
                cav_args[j + 1] = initial_gsd
                break
    args = (
        [sys.executable, str(script_dir / "run_cav_hoomd_advanced.py")]
        + cav_args
        + ["--replicas", str(replica_id)]
        + ["--frame", str(frame)]
        + ["--seed", str(base_seed + replica_id)]
    )
    result = subprocess.run(args, cwd=script_dir)
    return result.returncode == 0


def welford_update(
    count: int, mean: float, m2: float, new_value: float
) -> Tuple[int, float, float]:
    """Single step of Welford's online algorithm; returns (new_count, new_mean, new_M2)."""
    count += 1
    delta = new_value - mean
    mean += delta / count
    delta2 = new_value - mean
    m2 += delta * delta2
    return (count, mean, m2)


def write_averaged_files(
    directory: Path,
    averaged_prefix: str,
    stats: Dict[Tuple[int, float], Tuple[int, float, float]],
    n_replicas: int,
) -> None:
    """
    stats: (ref_idx, lag_ps) -> (count, mean, M2).
    Write {averaged_prefix}_fkt_ref_{ref_idx:03d}.txt with lag_time_ps, mean, std.
    """
    # Group by ref_idx
    by_ref: Dict[int, List[Tuple[float, float, float]]] = {}  # ref -> [(lag_ps, mean, std)]
    for (ref_idx, lag_ps), (count, mean, m2) in stats.items():
        if count < 1:
            continue
        variance = m2 / (count - 1) if count > 1 else 0.0
        std = variance ** 0.5
        if ref_idx not in by_ref:
            by_ref[ref_idx] = []
        by_ref[ref_idx].append((lag_ps, mean, std))

    for ref_idx in sorted(by_ref.keys()):
        rows = by_ref[ref_idx]
        rows.sort(key=lambda r: r[0])
        path = directory / f"{averaged_prefix}_fkt_ref_{ref_idx:03d}.txt"
        with open(path, "w", encoding="utf-8") as f:
            f.write("# F(k,t) correlation function (averaged over replicas)\n")
            f.write(f"# N replicas: {n_replicas}\n")
            f.write("# lag_time_ps\tF(k,t)_mean\tstd\n")
            for lag_ps, mean, std in rows:
                f.write(f"{lag_ps:.6f}\t{mean:.8f}\t{std:.8f}\n")


def run_single_replica(
    script_dir: Path,
    replica_id: int,
    replica_prefix: str,
    base_seed: int,
    initial_gsd: Optional[str] = None,
    initial_gsd_frame: int = 0,
    gsd_in_nm: bool = False,
) -> bool:
    """Run run_simulation.py once with --fkt-output-prefix and --seed; return True on success.
    If initial_gsd is set, pass --initial-gsd and --initial-gsd-frame so OpenMM uses the same
    initial configuration as cav-hoomd: frame = initial_gsd_frame if >= 0, else replica_id.
    If gsd_in_nm is True, pass --gsd-in-nm so OpenMM does not convert from Bohr (GSD already in nm).
    """
    args = (
        [sys.executable, str(script_dir / "run_simulation.py")]
        + DEFAULT_RUN_SIMULATION_ARGS
        + ["--fkt-output-prefix", f"{replica_prefix}_{replica_id:04d}"]
        + ["--seed", str(base_seed + replica_id)]
    )
    if initial_gsd:
        frame = replica_id if initial_gsd_frame < 0 else initial_gsd_frame
        args += ["--initial-gsd", initial_gsd, "--initial-gsd-frame", str(frame)]
    if gsd_in_nm:
        args += ["--gsd-in-nm"]
    result = subprocess.run(args, cwd=script_dir)
    return result.returncode == 0


def _run_openmm_loop(
    script_dir: Path,
    args: Any,
) -> int:
    """Run OpenMM replicas in sequence; update stats_openmm and write averaged files. Returns n_success."""
    stats: Dict[Tuple[int, float], Tuple[int, float, float]] = {}
    n_success = 0
    for i in range(args.replicas):
        ok = run_single_replica(
            script_dir,
            i,
            args.replica_prefix,
            args.base_seed,
            initial_gsd=args.initial_gsd,
            gsd_in_nm=args.gsd_in_nm,
            initial_gsd_frame=args.initial_gsd_frame,
        )
        if not ok:
            print(f"[OpenMM] WARNING: Replica {i} failed; skipping F(k,t) merge.")
            continue
        n_success += 1
        data = parse_replica_fkt_files(script_dir, args.replica_prefix, i)
        for ref_idx, rows in data.items():
            for lag_ps, fkt in rows:
                key = (ref_idx, _lag_key(lag_ps))
                count, mean, m2 = stats.get(key, (0, 0.0, 0.0))
                count, mean, m2 = welford_update(count, mean, m2, fkt)
                stats[key] = (count, mean, m2)
        if n_success % args.write_every == 0:
            write_averaged_files(
                script_dir, args.averaged_openmm_prefix, stats, n_success
            )
        if (n_success > 0) and (n_success % args.inspect_every == 0):
            ref_indices = sorted({k[0] for k in stats})
            if ref_indices:
                ref0 = ref_indices[0]
                parts = [f"[OpenMM] Replica {n_success}/{args.replicas} done;"]
                for lag in (0.0, 100.0, 500.0):
                    key = (ref0, _lag_key(lag))
                    if key in stats:
                        count, mean, m2 = stats[key]
                        var = m2 / (count - 1) if count > 1 else 0.0
                        se = (var / count) ** 0.5 if count else 0.0
                        parts.append(
                            f" F(k,t)(lag={lag:.0f} ps)= {mean:.2f} ± {se:.2f} (n={count})"
                        )
                print(" ".join(parts))
    if n_success > 0 and n_success % args.write_every != 0:
        write_averaged_files(
            script_dir, args.averaged_openmm_prefix, stats, n_success
        )
    return n_success


def _run_hoomd_loop(
    script_dir: Path,
    args: Any,
) -> int:
    """Run cav-hoomd replicas in sequence; update stats_hoomd and write averaged files. Returns n_success."""
    stats: Dict[Tuple[int, float], Tuple[int, float, float]] = {}
    n_success = 0
    for i in range(args.replicas):
        ok = run_single_hoomd_replica(
            script_dir,
            i,
            args.base_seed,
            initial_gsd=args.initial_gsd or "init-0.gsd",
            initial_gsd_frame=args.initial_gsd_frame,
        )
        if not ok:
            print(f"[HOOMD] WARNING: Replica {i} failed; skipping F(k,t) merge.")
            continue
        n_success += 1
        data = parse_hoomd_replica_fkt_files(script_dir, i)
        for ref_idx, rows in data.items():
            for lag_ps, fkt in rows:
                key = (ref_idx, _lag_key(lag_ps))
                count, mean, m2 = stats.get(key, (0, 0.0, 0.0))
                count, mean, m2 = welford_update(count, mean, m2, fkt)
                stats[key] = (count, mean, m2)
        if n_success % args.write_every == 0:
            write_averaged_files(
                script_dir, args.averaged_hoomd_prefix, stats, n_success
            )
        if (n_success > 0) and (n_success % args.inspect_every == 0):
            ref_indices = sorted({k[0] for k in stats})
            if ref_indices:
                ref0 = ref_indices[0]
                parts = [f"[HOOMD] Replica {n_success}/{args.replicas} done;"]
                for lag in (0.0, 100.0, 500.0):
                    key = (ref0, _lag_key(lag))
                    if key in stats:
                        count, mean, m2 = stats[key]
                        var = m2 / (count - 1) if count > 1 else 0.0
                        se = (var / count) ** 0.5 if count else 0.0
                        parts.append(
                            f" F(k,t)(lag={lag:.0f} ps)= {mean:.2f} ± {se:.2f} (n={count})"
                        )
                print(" ".join(parts))
    if n_success > 0 and n_success % args.write_every != 0:
        write_averaged_files(
            script_dir, args.averaged_hoomd_prefix, stats, n_success
        )
    return n_success


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run OpenMM and cav-hoomd side by side and average F(k,t) on the fly.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--replicas",
        type=int,
        default=500,
        help="Number of replica runs per backend (default: 500)",
    )
    parser.add_argument(
        "--replica-prefix",
        type=str,
        default="fkt_replica",
        help="Prefix for per-replica OpenMM F(k,t) files (default: fkt_replica)",
    )
    parser.add_argument(
        "--averaged-openmm-prefix",
        type=str,
        default="fkt_averaged_openmm",
        help="Prefix for averaged OpenMM F(k,t) files (default: fkt_averaged_openmm)",
    )
    parser.add_argument(
        "--averaged-hoomd-prefix",
        type=str,
        default="fkt_averaged_hoomd",
        help="Prefix for averaged cav-hoomd F(k,t) files (default: fkt_averaged_hoomd)",
    )
    parser.add_argument(
        "--write-every",
        type=int,
        default=1,
        help="Write averaged files after every N replicas (default: 1)",
    )
    parser.add_argument(
        "--inspect-every",
        type=int,
        default=10,
        help="Print summary every N replicas (default: 10)",
    )
    parser.add_argument(
        "--base-seed",
        type=int,
        default=42,
        help="Random seed for replica 0; replica i uses base_seed + i (default: 42)",
    )
    parser.add_argument(
        "--no-openmm",
        action="store_true",
        help="Do not run OpenMM; only run cav-hoomd",
    )
    parser.add_argument(
        "--no-hoomd",
        action="store_true",
        help="Do not run cav-hoomd; only run OpenMM",
    )
    parser.add_argument(
        "--initial-gsd",
        type=str,
        default="init-0.gsd",
        dest="initial_gsd",
        help="GSD file for initial positions/box; used by both OpenMM and cav-hoomd (default: init-0.gsd). Set to empty to disable.",
    )
    parser.add_argument(
        "--initial-gsd-frame",
        type=int,
        default=0,
        dest="initial_gsd_frame",
        help="Frame index to use from --initial-gsd. Default 0 = same initial condition for all replicas. Use -1 to use replica_id as frame (different config per replica).",
    )
    parser.add_argument(
        "--gsd-in-nm",
        action="store_true",
        dest="gsd_in_nm",
        help="Pass --gsd-in-nm to OpenMM: GSD positions/box are already in nm (do not convert from Bohr). Use if init-0.gsd was written in nm.",
    )
    parser.add_argument("--no-initial-gsd", action="store_true", dest="no_initial_gsd", help="Do not load from GSD; OpenMM uses generated positions (Phase 3 debug).")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent

    # Normalize: empty string or --no-initial-gsd means no GSD
    if args.no_initial_gsd:
        args.initial_gsd = None
    elif args.initial_gsd and args.initial_gsd.strip() == "":
        args.initial_gsd = None
    if args.initial_gsd and not (script_dir / args.initial_gsd).exists() and not args.no_hoomd:
        print(f"WARNING: {args.initial_gsd} not found. Ensure it exists for cav-hoomd (e.g. from a prior cavity run or create molecular-0.gsd then add cavity).")

    print(f"Running {args.replicas} replicas per backend (OpenMM and cav-hoomd side by side).")
    print("OpenMM args: " + " ".join(DEFAULT_RUN_SIMULATION_ARGS))
    print("Cav-hoomd args: " + " ".join(DEFAULT_RUN_CAV_HOOMD_ARGS))

    openmm_result: List[int] = []
    hoomd_result: List[int] = []

    def run_openmm_and_store() -> None:
        openmm_result.append(_run_openmm_loop(script_dir, args))

    def run_hoomd_and_store() -> None:
        hoomd_result.append(_run_hoomd_loop(script_dir, args))

    t_openmm: Optional[threading.Thread] = None
    t_hoomd: Optional[threading.Thread] = None

    if not args.no_openmm:
        t_openmm = threading.Thread(target=run_openmm_and_store, daemon=False)
        t_openmm.start()
    if not args.no_hoomd:
        t_hoomd = threading.Thread(target=run_hoomd_and_store, daemon=False)
        t_hoomd.start()

    if t_openmm is not None:
        t_openmm.join()
    if t_hoomd is not None:
        t_hoomd.join()

    n_openmm = openmm_result[0] if openmm_result else 0
    n_hoomd = hoomd_result[0] if hoomd_result else 0
    print(f"Completed: OpenMM {n_openmm}/{args.replicas}, cav-hoomd {n_hoomd}/{args.replicas} successful replicas.")
    ok_openmm = (args.no_openmm or n_openmm == args.replicas)
    ok_hoomd = (args.no_hoomd or n_hoomd == args.replicas)
    return 0 if (ok_openmm and ok_hoomd) else 1


if __name__ == "__main__":
    sys.exit(main())

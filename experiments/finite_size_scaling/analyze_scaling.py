#!/usr/bin/env python3
"""
Post-processing for the finite-size scaling campaign.

Reads F(k,t) files from all replicas at each system size N, computes
ensemble-averaged F(k,t), extracts structural relaxation times tau_s,
and produces finite-size scaling plots.

Usage:
  python analyze_scaling.py                        # analyze all available sizes
  python analyze_scaling.py --sizes 250 500 1000   # specific sizes
  python analyze_scaling.py --plot-only             # skip re-averaging, just plot
"""

import argparse
import glob
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from config import (
    SYSTEM_SIZES, RESULTS_DIR, SWITCH_TIME_PS,
    FKT_REF_INTERVAL_PS, lambda_scaled, num_replicas,
)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_fkt_file(path: Path) -> Tuple[float, List[Tuple[float, float]]]:
    """
    Parse one F(k,t) file.

    Returns (reference_time_ps, [(lag_ps, fkt_value), ...]).
    """
    ref_time = None
    data = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                m = re.search(r"Reference time:\s+([\d.]+)", line)
                if m:
                    ref_time = float(m.group(1))
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                try:
                    data.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    continue
    return ref_time, data


def collect_fkt_files(n_molecules: int, tag: str = "cavity") -> Dict[int, List[Path]]:
    """
    Collect F(k,t) files for all replicas at system size N.

    Returns dict mapping ref_file_idx -> list of paths (one per replica).
    """
    base = RESULTS_DIR / f"N{n_molecules}" / tag
    if not base.exists():
        return {}

    by_ref: Dict[int, List[Path]] = {}
    for replica_dir in sorted(base.iterdir()):
        if not replica_dir.is_dir():
            continue
        for fkt_file in sorted(replica_dir.glob("*_fkt_ref_*.txt")):
            m = re.search(r"_fkt_ref_(\d+)\.txt$", fkt_file.name)
            if m:
                ref_idx = int(m.group(1))
                by_ref.setdefault(ref_idx, []).append(fkt_file)

    return by_ref


# ---------------------------------------------------------------------------
# Averaging (Welford's online algorithm)
# ---------------------------------------------------------------------------

def average_fkt_across_replicas(
    files: List[Path],
) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Compute mean and stderr of F(k,t) across replicas for one reference index.

    Returns (ref_time_ps, lag_times, mean_fkt, stderr_fkt, n_replicas).
    """
    ref_time = None
    lag_dict: Dict[float, List[float]] = {}

    for fpath in files:
        rt, data = parse_fkt_file(fpath)
        if rt is not None:
            ref_time = rt
        for lag, fkt in data:
            key = round(lag, 4)
            lag_dict.setdefault(key, []).append(fkt)

    if not lag_dict:
        return ref_time, np.array([]), np.array([]), np.array([]), 0

    lags_sorted = sorted(lag_dict.keys())
    lags = np.array(lags_sorted)
    means = np.array([np.mean(lag_dict[k]) for k in lags_sorted])
    stds = np.array([np.std(lag_dict[k], ddof=1) if len(lag_dict[k]) > 1 else 0.0
                     for k in lags_sorted])
    counts = np.array([len(lag_dict[k]) for k in lags_sorted])
    stderr = stds / np.sqrt(counts)

    return ref_time, lags, means, stderr, int(np.median(counts))


# ---------------------------------------------------------------------------
# Relaxation time extraction
# ---------------------------------------------------------------------------

def extract_tau_s(lag_times: np.ndarray, fkt_mean: np.ndarray,
                  threshold: float = 0.1) -> Optional[float]:
    """
    Extract relaxation time tau_s defined by F(k, tau_s) = threshold.

    Normalizes F(k,t) so that F(k,0) = 1 (ISF -> phi).
    Returns None if the ISF never drops below the threshold.
    """
    if len(lag_times) < 2 or fkt_mean[0] == 0:
        return None

    phi = fkt_mean / fkt_mean[0]

    below = np.where(phi <= threshold)[0]
    if len(below) == 0:
        return None

    idx = below[0]
    if idx == 0:
        return lag_times[0]

    # Linear interpolation between points straddling the threshold
    t0, t1 = lag_times[idx - 1], lag_times[idx]
    p0, p1 = phi[idx - 1], phi[idx]
    tau = t0 + (threshold - p0) * (t1 - t0) / (p1 - p0)
    return float(tau)


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def analyze_system_size(n_molecules: int, output_dir: Path) -> dict:
    """Analyze all F(k,t) data for one system size."""
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {
        "N": n_molecules,
        "lambda": lambda_scaled(n_molecules),
        "waiting_times": [],
        "tau_s": [],
        "n_replicas_actual": [],
    }

    # Cavity runs
    fkt_by_ref = collect_fkt_files(n_molecules, "cavity")
    if not fkt_by_ref:
        print(f"  [N={n_molecules}] No cavity F(k,t) data found")
        return results

    for ref_idx in sorted(fkt_by_ref.keys()):
        files = fkt_by_ref[ref_idx]
        ref_time, lags, mean_fkt, stderr_fkt, n_rep = average_fkt_across_replicas(files)

        if len(lags) == 0:
            continue

        # Waiting time relative to switch time
        if ref_time is not None:
            t_w = ref_time - SWITCH_TIME_PS
        else:
            t_w = ref_idx * FKT_REF_INTERVAL_PS

        # Save averaged F(k,t)
        out_file = output_dir / f"fkt_avg_N{n_molecules}_tw{t_w:.0f}.npz"
        np.savez(out_file, lag_ps=lags, fkt_mean=mean_fkt, fkt_stderr=stderr_fkt,
                 ref_time_ps=ref_time, t_w_ps=t_w, n_replicas=n_rep)

        tau = extract_tau_s(lags, mean_fkt)
        results["waiting_times"].append(t_w)
        results["tau_s"].append(tau)
        results["n_replicas_actual"].append(n_rep)

        tau_str = f"{tau:.2f} ps" if tau is not None else "N/A (not decayed)"
        print(f"  [N={n_molecules}] t_w={t_w:7.0f} ps  tau_s={tau_str:>16s}  ({n_rep} replicas)")

    # Baseline
    base_fkt = collect_fkt_files(n_molecules, "baseline")
    if base_fkt:
        tau_s_baseline = []
        for ref_idx in sorted(base_fkt.keys()):
            _, lags, mean_fkt, _, _ = average_fkt_across_replicas(base_fkt[ref_idx])
            tau = extract_tau_s(lags, mean_fkt)
            if tau is not None:
                tau_s_baseline.append(tau)
        if tau_s_baseline:
            results["tau_s_baseline"] = float(np.mean(tau_s_baseline))

    return results


def plot_scaling(all_results: List[dict], output_dir: Path):
    """Generate finite-size scaling plots."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plots")
        return

    # --- Plot 1: tau_s_tilde vs N at tw=0 ---
    fig, ax = plt.subplots(figsize=(8, 5))

    ns, taus_tw0, taus_baseline = [], [], []
    for r in all_results:
        if not r["tau_s"]:
            continue
        # First waiting time after switch
        idx_tw0 = None
        for i, tw in enumerate(r["waiting_times"]):
            if tw >= 0:
                idx_tw0 = i
                break
        if idx_tw0 is not None and r["tau_s"][idx_tw0] is not None:
            ns.append(r["N"])
            taus_tw0.append(r["tau_s"][idx_tw0])
            taus_baseline.append(r.get("tau_s_baseline", r["tau_s"][idx_tw0]))

    if ns:
        tau_tilde = [t / b if b else None for t, b in zip(taus_tw0, taus_baseline)]
        valid = [(n, tt) for n, tt in zip(ns, tau_tilde) if tt is not None]
        if valid:
            ns_v, tt_v = zip(*valid)
            ax.semilogx(ns_v, tt_v, "o-", color="C0", markersize=8)
            ax.set_xlabel("System size N (molecules)", fontsize=12)
            ax.set_ylabel(r"$\tilde{\tau}_s = \tau_s(\lambda) / \tau_s(0)$", fontsize=12)
            ax.set_title(r"Finite-size scaling at $t_w = 0$", fontsize=13)
            ax.axhline(1.0, color="gray", ls="--", alpha=0.5)
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            fig.savefig(output_dir / "tau_tilde_vs_N.png", dpi=200)
            print(f"  Saved {output_dir / 'tau_tilde_vs_N.png'}")
    plt.close(fig)

    # --- Plot 2: tau_s_tilde vs t_w for each N ---
    fig, ax = plt.subplots(figsize=(8, 5))
    cmap = plt.cm.viridis
    valid_results = [r for r in all_results if r["tau_s"]]
    for i, r in enumerate(valid_results):
        color = cmap(i / max(len(valid_results) - 1, 1))
        tws = np.array(r["waiting_times"])
        taus = np.array(r["tau_s"], dtype=float)
        baseline = r.get("tau_s_baseline", taus[0] if len(taus) > 0 else 1.0)
        mask = np.isfinite(taus) & (tws >= 0)
        if mask.any() and baseline:
            ax.plot(tws[mask], taus[mask] / baseline, "o-", color=color,
                    label=f"N={r['N']}", markersize=5)

    ax.set_xlabel(r"$t_w$ (ps)", fontsize=12)
    ax.set_ylabel(r"$\tilde{\tau}_s$", fontsize=12)
    ax.set_title("Non-thermal aging at constant effective coupling", fontsize=13)
    ax.axhline(1.0, color="gray", ls="--", alpha=0.5)
    ax.legend(fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_dir / "tau_tilde_vs_tw.png", dpi=200)
    print(f"  Saved {output_dir / 'tau_tilde_vs_tw.png'}")
    plt.close(fig)

    # --- Plot 3: raw tau_s vs N ---
    fig, ax = plt.subplots(figsize=(8, 5))
    if ns:
        ax.loglog(ns, taus_tw0, "s-", color="C1", markersize=8)
        ax.set_xlabel("System size N", fontsize=12)
        ax.set_ylabel(r"$\tau_s$ at $t_w = 0$ (ps)", fontsize=12)
        ax.set_title("Raw relaxation time vs system size", fontsize=13)
        ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(output_dir / "tau_s_vs_N.png", dpi=200)
    print(f"  Saved {output_dir / 'tau_s_vs_N.png'}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Analyze finite-size scaling results")
    parser.add_argument("--sizes", type=int, nargs="*", default=None,
                        help="System sizes to analyze (default: all available)")
    parser.add_argument("--plot-only", action="store_true",
                        help="Skip re-averaging, load existing NPZ and plot")
    args = parser.parse_args()

    sizes = args.sizes or SYSTEM_SIZES
    output_dir = RESULTS_DIR / "analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = []

    for n in sizes:
        print(f"\n--- N = {n} ---")
        sys_output = output_dir / f"N{n}"

        if args.plot_only:
            npz_files = sorted(sys_output.glob("fkt_avg_*.npz")) if sys_output.exists() else []
            if not npz_files:
                print(f"  No pre-computed data for N={n}")
                continue
            result = {"N": n, "lambda": lambda_scaled(n),
                      "waiting_times": [], "tau_s": [], "n_replicas_actual": []}
            for npz_file in npz_files:
                data = np.load(npz_file)
                tw = float(data["t_w_ps"])
                tau = extract_tau_s(data["lag_ps"], data["fkt_mean"])
                result["waiting_times"].append(tw)
                result["tau_s"].append(tau)
                result["n_replicas_actual"].append(int(data["n_replicas"]))
            all_results.append(result)
        else:
            result = analyze_system_size(n, sys_output)
            all_results.append(result)

    # Summary table
    print("\n" + "=" * 70)
    print(f"{'N':>7s}  {'lambda':>8s}  {'replicas':>8s}  {'tau_s(tw=0)':>14s}  {'tau_s(baseline)':>15s}")
    print("-" * 70)
    for r in all_results:
        tau_tw0 = None
        for i, tw in enumerate(r["waiting_times"]):
            if tw >= 0 and r["tau_s"][i] is not None:
                tau_tw0 = r["tau_s"][i]
                break
        n_rep = r["n_replicas_actual"][0] if r["n_replicas_actual"] else 0
        tau_str = f"{tau_tw0:.2f} ps" if tau_tw0 is not None else "N/A"
        base_str = f"{r.get('tau_s_baseline', 'N/A')}"
        if isinstance(r.get("tau_s_baseline"), float):
            base_str = f"{r['tau_s_baseline']:.2f} ps"
        print(f"{r['N']:7d}  {r['lambda']:8.4f}  {n_rep:8d}  {tau_str:>14s}  {base_str:>15s}")

    # Save summary
    summary_file = output_dir / "scaling_summary.npz"
    np.savez(summary_file,
             system_sizes=np.array([r["N"] for r in all_results]),
             lambdas=np.array([r["lambda"] for r in all_results]),
             data=all_results)
    print(f"\nSummary saved to {summary_file}")

    # Plots
    print("\nGenerating plots ...")
    plot_scaling(all_results, output_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()

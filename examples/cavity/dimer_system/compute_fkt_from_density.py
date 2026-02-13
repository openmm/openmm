#!/usr/bin/env python3
"""
Compute F(k,t) from bare density npz via equilibrium time-averaging and k-averaging.

Loads *_density.npz produced by run_simulation.py (when --enable-fkt and density npz
is written). Assumes the retained segment is in equilibrium and computes

  F(τ) = (1/n_k) Σ_k C_k(τ) / C_k(0)
  C_k(τ) = (1/(n-τ)) Σ_{t'} Re( ρ_k(t'+τ) ρ*_k(t') )

Output: npz and optionally text file with lag_time_ps, Fkt; optional plot.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def load_density_npz(path: Path) -> dict:
    """Load density npz; return dict with time_ps, rhok_real, rhok_imag, wavevectors, etc."""
    data = np.load(path, allow_pickle=False)
    required = ("time_ps", "rhok_real", "rhok_imag", "wavevectors")
    for k in required:
        if k not in data:
            raise ValueError(f"Missing '{k}' in {path}. Expected keys: {list(data.keys())}")
    return {k: data[k] for k in data.files}


def compute_fkt_from_density(
    time_ps: np.ndarray,
    rhok_real: np.ndarray,
    rhok_imag: np.ndarray,
    equil_ps: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute k-averaged, time-averaged F(k,t) normalized to F(0)=1.

    Parameters
    ----------
    time_ps : np.ndarray
        Shape (n_times,), simulation times in ps.
    rhok_real, rhok_imag : np.ndarray
        Shape (n_times, n_k), real and imaginary parts of ρ_k(t).
    equil_ps : float
        Discard frames with time_ps < equil_ps for time-averaging (default 0).

    Returns
    -------
    lag_time_ps : np.ndarray
        Shape (n_lag,), lag times in ps.
    Fkt : np.ndarray
        Shape (n_lag,), F(k,t) normalized so Fkt[0]=1.
    """
    rhok = rhok_real + 1j * rhok_imag  # (n_times, n_k)
    n_times, n_k = rhok.shape

    # Equilibration: keep frames with time_ps >= equil_ps
    if equil_ps > 0:
        mask = time_ps >= equil_ps
        if not np.any(mask):
            raise ValueError(
                f"equil_ps={equil_ps} ps leaves no frames (time_ps range [{time_ps.min():.2f}, {time_ps.max():.2f}])."
            )
        time_ps = time_ps[mask]
        rhok = rhok[mask]
        n_times = rhok.shape[0]

    # For each lag τ, C_k(τ) = (1/(n_times-τ)) * Σ_{t'} Re(rhok[t'+τ,k] * conj(rhok[t',k]))
    # F(τ) = (1/n_k) * Σ_k C_k(τ)
    F_raw = np.zeros(n_times, dtype=np.float64)
    for tau in range(n_times):
        n_sum = n_times - tau
        if n_sum <= 0:
            break
        # corr[t',k] = Re(rhok[t'+τ,k] * conj(rhok[t',k])) = Re(rhok[t'+τ,k])*Re(rhok[t',k]) + Im(rhok[t'+τ,k])*Im(rhok[t',k])
        # shape (n_sum,) for each k, then mean over k
        prod = np.real(rhok[tau:] * np.conj(rhok[:n_sum]))  # (n_sum, n_k)
        C_k = np.mean(prod, axis=0)  # (n_k,)
        F_raw[tau] = np.mean(C_k)

    # Normalize so F(0)=1
    if F_raw[0] <= 0:
        raise ValueError(
            f"F(0) = {F_raw[0]:.6e} <= 0; cannot normalize. Check density data."
        )
    Fkt = F_raw / F_raw[0]

    # Lag times: use mean dt from production segment
    if n_times > 1:
        dt_ps = np.mean(np.diff(time_ps))
    else:
        dt_ps = 0.0
    lag_time_ps = np.arange(n_times, dtype=np.float64) * dt_ps

    return lag_time_ps, Fkt


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute F(k,t) from bare density npz (equilibrium time- and k-averaging)."
    )
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        required=True,
        dest="input_path",
        help="Path to *_density.npz from run_simulation.py",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        dest="output_path",
        help="Output path for F(k,t) (default: {input_stem}_fkt_from_density.npz)",
    )
    parser.add_argument(
        "--equil-ps",
        type=float,
        default=0.0,
        dest="equil_ps",
        help="Discard this many ps from the start for time-averaging (default: 0)",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Write a simple plot lag_time_ps vs F(k,t)",
    )
    parser.add_argument(
        "--text",
        action="store_true",
        help="Also write two-column text file (lag_time_ps, Fkt)",
    )
    args = parser.parse_args()

    path = Path(args.input_path).resolve()
    if not path.exists():
        print(f"Error: file not found: {path}")
        return 1

    data = load_density_npz(path)
    time_ps = np.asarray(data["time_ps"], dtype=np.float64)
    rhok_real = np.asarray(data["rhok_real"], dtype=np.float64)
    rhok_imag = np.asarray(data["rhok_imag"], dtype=np.float64)

    n_times_orig = time_ps.shape[0]
    lag_time_ps, Fkt = compute_fkt_from_density(
        time_ps, rhok_real, rhok_imag, equil_ps=args.equil_ps
    )
    n_times_used = len(lag_time_ps)
    n_wavevectors = rhok_real.shape[1]

    out_stem = path.stem.replace("_density", "") if path.stem.endswith("_density") else path.stem
    out_path = args.output_path
    if out_path is None:
        out_path = path.parent / f"{out_stem}_fkt_from_density.npz"
    out_path = Path(out_path)

    np.savez(
        out_path,
        lag_time_ps=lag_time_ps,
        Fkt=Fkt,
        equil_ps_used=np.float64(args.equil_ps),
        n_times_used=np.int64(n_times_used),
        n_times_total=np.int64(n_times_orig),
        n_wavevectors=np.int64(n_wavevectors),
    )
    print(f"Wrote {out_path}")
    print(f"  equil_ps_used: {args.equil_ps}")
    print(f"  n_times_used: {n_times_used} (of {n_times_orig})")
    print(f"  n_wavevectors: {n_wavevectors}")
    print(f"  Fkt[0] = {Fkt[0]:.6f}")

    if args.text:
        txt_path = out_path.with_suffix(".txt")
        with open(txt_path, "w", encoding="utf-8") as f:
            f.write("# lag_time_ps\tF(k,t)\n")
            f.write(f"# equil_ps_used={args.equil_ps}, n_times_used={n_times_used}, n_wavevectors={n_wavevectors}\n")
            for t, fval in zip(lag_time_ps, Fkt):
                f.write(f"{t:.6f}\t{fval:.8f}\n")
        print(f"Wrote {txt_path}")

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("Warning: matplotlib not available, skipping --plot")
        else:
            fig_path = out_path.with_suffix(".png")
            plt.figure(figsize=(6, 4))
            plt.plot(lag_time_ps, Fkt, "b-", linewidth=1)
            plt.xlabel("Lag time (ps)")
            plt.ylabel("F(k,t)")
            plt.title(f"F(k,t) from density npz (equil_ps={args.equil_ps}, n={n_times_used})")
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(fig_path, dpi=150)
            plt.close()
            print(f"Wrote {fig_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

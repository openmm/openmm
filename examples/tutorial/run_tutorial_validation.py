#!/usr/bin/env python3
"""Validate mKA cavity MD tutorial physics (Section 2 + spectrum).

Run from repo root after building OpenMM:

    python examples/tutorial/run_tutorial_validation.py
    python examples/tutorial/run_tutorial_validation.py --steps 20000
"""

from __future__ import annotations

import argparse
import sys

from tutorial_common import OMEGA_C_CM1, run_nvt_single_dimer


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lambda-coupling", type=float, default=0.01)
    parser.add_argument("--temperature-K", type=float, default=100.0)
    parser.add_argument("--steps", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--platform", default="CPU")
    parser.add_argument(
        "--peak-tolerance-cm1",
        type=float,
        default=400.0,
        help="Allowed deviation of spectral peak from cavity frequency (cm^-1)",
    )
    parser.add_argument(
        "--temperature-tolerance-K",
        type=float,
        default=40.0,
        help="Allowed deviation of mean molecular T from bath (K)",
    )
    args = parser.parse_args()

    print("Running tutorial Section 2 validation...")
    print(
        f"  lambda={args.lambda_coupling}, T={args.temperature_K} K, "
        f"steps={args.steps}, platform={args.platform}"
    )

    result = run_nvt_single_dimer(
        lambda_coupling=args.lambda_coupling,
        temperature_K=args.temperature_K,
        n_steps=args.steps,
        seed=args.seed,
        platform_name=args.platform,
    )

    mean_t = result["mean_temperature_K"]
    mean_t_ph = result["mean_photon_temperature_K"]
    peak = result["peak_frequency_cm1"]
    omega_c = result["omega_c_cm1"]

    print(f"\nResults:")
    print(f"  Mean molecular T_kin  = {mean_t:.1f} K  (target {args.temperature_K} K)")
    print(f"  Mean photon T_kin     = {mean_t_ph:.1f} K")
    print(f"  Spectral peak         = {peak:.0f} cm^-1  (cavity omega_c = {omega_c:.0f} cm^-1)")

    failures = []
    if abs(mean_t - args.temperature_K) > args.temperature_tolerance_K:
        failures.append(
            f"molecular temperature {mean_t:.1f} K deviates from "
            f"{args.temperature_K} K by more than {args.temperature_tolerance_K} K"
        )
    if mean_t_ph > 5.0 * args.temperature_K:
        failures.append(
            f"photon temperature {mean_t_ph:.1f} K is much higher than bath "
            f"({args.temperature_K} K) — check displaceToEquilibrium"
        )
    if abs(peak - omega_c) > args.peak_tolerance_cm1:
        failures.append(
            f"spectral peak {peak:.0f} cm^-1 deviates from omega_c={omega_c:.0f} cm^-1 "
            f"by more than {args.peak_tolerance_cm1} cm^-1"
        )

    if failures:
        print("\nFAIL:")
        for msg in failures:
            print(f"  - {msg}")
        return 1

    print("\nPASS: tutorial physics checks succeeded.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

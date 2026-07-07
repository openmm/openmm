#!/usr/bin/env python3
"""Regression tests for mKA cavity MD tutorial physics."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

TUTORIAL_DIR = Path(__file__).resolve().parents[2] / "examples" / "tutorial"
sys.path.insert(0, str(TUTORIAL_DIR))

openmm = pytest.importorskip("openmm")

from tutorial_common import OMEGA_C_CM1, run_nvt_single_dimer, select_platform  # noqa: E402


@pytest.mark.timeout(300)
def test_nvt_temperature_and_spectrum_peak():
    """Section 2: molecular T near bath and spectral peak near omega_c."""
    result = run_nvt_single_dimer(
        lambda_coupling=0.01,
        temperature_K=100.0,
        n_steps=3000,
        seed=42,
        platform_name=select_platform().getName(),
    )

    assert abs(result["mean_system_temperature_K"] - 100.0) < 50.0, (
        f"Mean system T={result['mean_system_temperature_K']:.1f} K too far from 100 K"
    )
    assert abs(result["peak_frequency_cm1"] - OMEGA_C_CM1) < 500.0, (
        f"Peak {result['peak_frequency_cm1']:.0f} cm^-1 too far from "
        f"omega_c={OMEGA_C_CM1:.0f} cm^-1"
    )


@pytest.mark.timeout(120)
def test_displace_to_equilibrium_is_small():
    """Photon equilibrium displacement should be tiny for a single dimer."""
    from tutorial_common import (
        Units,
        build_single_aa_dimer_system,
        create_context,
    )

    omegac_au = Units.cm1_to_au(OMEGA_C_CM1)
    system, displacer, positions = build_single_aa_dimer_system(
        lambda_coupling=0.01, omegac_au=omegac_au
    )
    context = create_context(system, dt_fs=1.0, temperature_K=100.0, seed=42)
    context.setPositions(positions)
    displacer.displaceToEquilibrium(context, 0.01)

    pos = context.getState(getPositions=True).getPositions(asNumpy=True)
    q_xy = float((pos[2, 0] ** 2 + pos[2, 1] ** 2) ** 0.5)
    assert q_xy < 1e-3, f"Photon displacement {q_xy:.4f} nm is too large"

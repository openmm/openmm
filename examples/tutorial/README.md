# mKA Cavity MD Tutorial

Step-by-step notebook and validation scripts for cavity molecular dynamics with the modified Kob–Andersen (mKA) dimer model.

## Contents

| File | Purpose |
|------|---------|
| [`mka_cavity_md_tutorial.ipynb`](mka_cavity_md_tutorial.ipynb) | Interactive walkthrough (Sections 0–5) |
| [`tutorial_common.py`](tutorial_common.py) | Shared system builders and analysis helpers |
| [`run_tutorial_validation.py`](run_tutorial_validation.py) | Headless physics validation script |

## Prerequisites

Build OpenMM from this repository with Python wrappers enabled:

```bash
pixi install
pixi run smoke   # verify import
```

For the notebook, use the test environment with Jupyter:

```bash
pixi run -e test jupyter lab examples/tutorial/mka_cavity_md_tutorial.ipynb
```

## Run validation (no notebook required)

From the repository root:

```bash
python examples/tutorial/run_tutorial_validation.py
python examples/tutorial/run_tutorial_validation.py --steps 20000
python examples/tutorial/run_tutorial_validation.py --platform CPU
```

By default the validation script uses **CUDA** (mixed precision) when available, otherwise **CPU**.

This runs tutorial Section 2 (single A–A dimer + photon, NVT Bussi) and checks:

- Mean molecular kinetic temperature near the bath (100 K by default)
- Photon in-plane temperature stays bounded (no runaway heating from bad equilibrium displacement)
- Dipole spectrum peak near the cavity frequency (1560 cm⁻¹)

## Automated tests

```bash
pixi run -e test test-py
```

The tutorial regression lives in [`tests/tutorial/test_mka_tutorial_physics.py`](../tests/tutorial/test_mka_tutorial_physics.py).

## Physics notes

- **Finite-q displacement**: `displaceToEquilibrium()` sets `q_eq = -(λ/(m_ph·ω_c))·d_xy`, matching `CavityForce` with `K = m_ph·ω_c²`.
- **Bussi thermostat**: Applied to molecular atoms only; photon equilibrates via cavity coupling.
- **Photon temperature**: The photon is **not** Bussi-thermostatted. At weak coupling (λ ≈ 0.01) its in-plane kinetic temperature is **below** the molecular bath (often ~50 K when T_bath = 100 K). This is expected, not a DOF miscount. Report T using 2 DOF in the cavity plane (y, z when the dimer lies along x); 3-DOF Cartesian T is similar.
- **Dipole self-energy**: Always included in `CavityForce`; no toggle is required in Python setup code.
- **Spectrum**: Direct FFT of the dipole trace; peak should be near ω_c = 1560 cm⁻¹ for λ = 0.01.
- **Platform**: The notebook prefers CUDA; CPU/Reference is supported after the split-Verlet + Bussi integrator fix in `ReferenceVerletDynamics`.

## Troubleshooting

| Symptom | Likely cause |
|---------|----------------|
| `AttributeError: setIncludeDipoleSelfEnergy` | Remove that call; DSE is always on |
| Photon T >> bath T | Old `displaceToEquilibrium` bug; rebuild OpenMM from current branch |
| Peak ~1000 cm⁻¹ instead of ~1560 cm⁻¹ (CPU) | Old Reference split-Verlet bug; rebuild OpenMM from current branch |
| Photon T ≈ 50 K with T_bath = 100 K | Expected at weak λ; photon is not thermostatted (not a DOF bug) |

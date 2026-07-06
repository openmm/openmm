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
```

This runs tutorial Section 2 (single A–A dimer + photon, NVT Bussi) and checks:

- Mean molecular kinetic temperature near the bath (100 K by default)
- Photon temperature stays bounded (no runaway heating from bad equilibrium displacement)
- Dipole spectrum peak near the cavity frequency (1560 cm⁻¹)

## Automated tests

```bash
pixi run -e test test-py
```

The tutorial regression lives in [`tests/tutorial/test_mka_tutorial_physics.py`](../tests/tutorial/test_mka_tutorial_physics.py).

## Physics notes

- **Finite-q displacement**: `displaceToEquilibrium()` sets `q_eq = -(λ/(m_ph·ω_c))·d_xy`, matching `CavityForce` with `K = m_ph·ω_c²`.
- **Bussi thermostat**: Applied to molecular atoms only; photon equilibrates via cavity coupling.
- **Dipole self-energy**: Always included in `CavityForce`; no toggle is required in Python setup code.
- **Spectrum**: Section 5 uses Wiener–Khinchin autocorrelation; peak should be near ω_c = 1560 cm⁻¹ for λ = 0.01.

## Troubleshooting

| Symptom | Likely cause |
|---------|----------------|
| `AttributeError: setIncludeDipoleSelfEnergy` | Remove that call; DSE is always on |
| Photon T >> bath T | Old `displaceToEquilibrium` bug; rebuild OpenMM from current branch |
| Peak ~1000 cm⁻¹ instead of ~1560 cm⁻¹ | Same displacement bug injecting excess PE |
| Mean T_kin ≈ 50 K with T_bath = 100 K | Using total KE or wrong DOF in temperature estimate |

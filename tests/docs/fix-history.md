# Test Suite Fix History

Technical notes on fixes applied to the dimer and cavity simulation tests.

## Force Constants

HOOMD and OpenMM both use `E = 0.5 * k * (r - r0)²`; use the same k values. Do not double k.
- Correct: O-O k=0.73204, N-N k=1.4325 (a.u.) → ~1560, ~2325 cm⁻¹
- Wrong: 2× those values → frequencies √2 too high

## Cavity Coupling

CavityForce must receive `lambda_coupling` from t=0, not 0.0 with `setCouplingOnStep`. Otherwise coupling stays off during equilibration.

## LJ Parameters (Dimer)

Kob-Andersen model needs non-additive N-O cross-terms; Lorentz-Berthelot fails. Use explicit exceptions with sigma_NO, epsilon_NO from cav-hoomd.

## IR Spectrum (ML Water)

- Use unwrapped positions when computing dipoles: omit `enforcePeriodicBox` so molecules across box edges stay intact.
- ACF: subtract mean dipole, FFT-based method, ω² spectral weight.
- Real-time: `compute_ir_realtime.py` reads dipole file while simulation runs.

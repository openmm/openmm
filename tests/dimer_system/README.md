# Dimer System

Two-component O-O and N-N dimer toy model for cavity coupling tests.

## RPMD Cavity Simulation

`run_simulation_rpmd.py` runs Ring Polymer MD of 250 dimers coupled to a cavity photon. Uses PILE_G thermostat (Bussi centroid + Langevin internal modes).

```bash
python run_simulation_rpmd.py --dimers 250 --beads 8 --lambda 0.07
```

| Option | Default | Description |
|--------|---------|-------------|
| `--dimers` | 250 | Number of dimers |
| `--beads` | 8 | RPMD beads |
| `--lambda` | 0.0700 | Cavity coupling |
| `--temp` | 100 | Temperature (K) |
| `--dt` | 0.001 | Timestep (ps) |
| `--cavity-freq` | 1560 | Cavity frequency (cm⁻¹) |

Output: `cavity_diamer_lambda{lambda:.4f}.npz` with `time_ps`, `dipole_nm`, `metadata`.

## Performance

| Beads | dt | Speed (ns/day) |
|-------|----|----------------|
| 8 | 1 fs | ~120 |
| 4 | 3 fs | ~918 |
| 4 | 4 fs | ~1220 (recommended) |

Use `--beads 4 --dt 0.004` for fast production; `--beads 4 --dt 0.003` for higher accuracy.

## Classical MD

`run_simulation.py` runs classical cavity dimer. See `tests/docs/fix-history.md` for force constants and cavity coupling notes.

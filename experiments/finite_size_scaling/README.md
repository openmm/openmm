# Finite-Size Scaling of Non-Thermal Aging

Studies how non-thermal aging (cavity-induced structural relaxation slowdown) depends on system size N, from N=250 to N=128,000 molecules. The single-molecule coupling lambda is scaled as lambda(N) = lambda_ref * sqrt(N_ref / N) to keep the effective collective coupling constant.

## Quick start

```bash
# 1. Print the campaign table
python config.py

# 2. Generate lattice configs for all sizes
python generate_initial_configs.py

# 3. Run equilibration (long, do on cluster -- see printed commands)

# 4. Run production locally (small sizes)
./run_local.sh 250

# 5. Or submit to SLURM (NYU Torch)
./submit_all.sh --dry-run     # preview
./submit_all.sh               # submit all sizes

# 6. Analyze results
python analyze_scaling.py
```

## Files

| File | Purpose |
|------|---------|
| `config.py` | All parameters: sizes, lambdas, replicas, paths |
| `generate_initial_configs.py` | Create GSD lattice configs per N |
| `run_single.py` | Run one (N, replica) simulation |
| `run_local.sh` | Local execution (sequential or GNU parallel) |
| `submit_slurm.sh` | SLURM array job template for NYU Torch |
| `submit_all.sh` | Submit all sizes as SLURM array jobs |
| `analyze_scaling.py` | Ensemble averaging, tau_s extraction, plots |

## Protocol

Following the paper:
- T = 100 K, omega_c = 1560 cm-1, dt = 1 fs
- Coupling activates at t0 = 200 ps (`--switch-time 200`)
- Production: 2500 ps total
- F(k,t) at k = 6.0 a.u. (113.4 nm-1), 50 wavevectors, t_w every 200 ps
- Bussi thermostat (tau = 1.0 ps), Langevin cavity bath

## Replica scaling

F(k,t) self-averages as 1/N, so fewer replicas are needed for larger systems:
N_T(N) = ceil(500 * 250 / N), minimum 10.

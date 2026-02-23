"""
Central configuration for the finite-size scaling campaign.

All system sizes, coupling strengths, replica counts, and simulation parameters
for studying non-thermal aging as a function of N.

Physics: the effective collective coupling is lambda * sqrt(N). To maintain the
same Rabi splitting (same effective coupling) across system sizes, we scale
lambda(N) = lambda_ref * sqrt(N_ref / N).

Replica scaling: F(k,t) is self-averaging with variance ~ 1/N per trajectory.
To maintain constant statistical precision, N_T(N) = N_T_ref * N_ref / N,
with a minimum floor.
"""

import math
from pathlib import Path

# ---------------------------------------------------------------------------
# Reference system (paper: 250 dimers, 500 atoms, box = 40 Bohr)
# ---------------------------------------------------------------------------
N_REF = 250
L_REF_BOHR = 40.0
LAMBDA_REF = 0.098        # a.u., paper's intermediate coupling
N_REPLICAS_REF = 500
MIN_REPLICAS = 10

# ---------------------------------------------------------------------------
# System sizes (powers of 2 from 250)
# ---------------------------------------------------------------------------
SYSTEM_SIZES = [250, 500, 1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000]

# ---------------------------------------------------------------------------
# Physical parameters (matching paper exactly)
# ---------------------------------------------------------------------------
TEMPERATURE_K = 100.0
CAVITY_FREQ_CM = 1560.0       # resonant with A2 (O-O) stretch
DT_PS = 0.001                 # 1 fs
SWITCH_TIME_PS = 200.0        # coupling activates at t0 = 200 ps
PRODUCTION_PS = 2500.0        # total runtime (200 ps equilibrium + 2300 ps after switch)
BUSSI_TAU_PS = 1.0            # Bussi thermostat time constant
FRACTION_OO = 0.8             # 4:1 A2:B2 mixture
DENSITY_BOHR = N_REF / L_REF_BOHR**3   # 0.00390625 molecules/Bohr^3

# F(k,t) parameters
FKT_KMAG = 113.4              # nm^-1 (= 6.0 a.u., first peak of S(k))
FKT_NUM_WAVEVECTORS = 50
FKT_REF_INTERVAL_PS = 200.0   # new t_w reference every 200 ps
FKT_MAX_REFS = 15
FKT_OUTPUT_PERIOD_PS = 0.1    # sample F(k,t) every 0.1 ps

# Equilibration for initial configs (no cavity)
EQUIL_SEPARATION_PS = 300.0   # separation between independent configs
EQUIL_TOTAL_NS = 150.0        # total equilibration run

# Unit conversions
BOHR_TO_NM = 0.0529177


def box_size_bohr(n_molecules: int) -> float:
    """Box edge length in Bohr at constant density."""
    return L_REF_BOHR * (n_molecules / N_REF) ** (1.0 / 3.0)


def box_size_nm(n_molecules: int) -> float:
    return box_size_bohr(n_molecules) * BOHR_TO_NM


def lambda_scaled(n_molecules: int) -> float:
    """Single-molecule coupling scaled to maintain constant effective coupling."""
    return LAMBDA_REF * math.sqrt(N_REF / n_molecules)


def num_replicas(n_molecules: int) -> int:
    """Number of replicas to maintain constant statistical precision."""
    raw = math.ceil(N_REPLICAS_REF * N_REF / n_molecules)
    return max(raw, MIN_REPLICAS)


def get_campaign_table():
    """Return list of dicts with all per-size parameters."""
    rows = []
    for n in SYSTEM_SIZES:
        rows.append({
            "N": n,
            "lambda": lambda_scaled(n),
            "box_bohr": box_size_bohr(n),
            "box_nm": box_size_nm(n),
            "replicas": num_replicas(n),
            "atoms": 2 * n + 1,  # 2 atoms per dimer + 1 cavity particle
        })
    return rows


# ---------------------------------------------------------------------------
# Directory layout
# ---------------------------------------------------------------------------
CAMPAIGN_ROOT = Path(__file__).resolve().parent
RUN_SIMULATION_SCRIPT = (CAMPAIGN_ROOT / "../../examples/cavity/dimer_system/run_simulation.py").resolve()
INITLATTICE_SCRIPT = (CAMPAIGN_ROOT / "../../examples/cavity/dimer_system/initlattice_equilibrium.py").resolve()
CONFIGS_DIR = CAMPAIGN_ROOT / "configs"
RESULTS_DIR = CAMPAIGN_ROOT / "results"
LOGS_DIR = CAMPAIGN_ROOT / "logs"


if __name__ == "__main__":
    print(f"{'N':>7s}  {'lambda':>8s}  {'L (Bohr)':>9s}  {'L (nm)':>8s}  {'Replicas':>8s}  {'Atoms':>7s}")
    print("-" * 60)
    for row in get_campaign_table():
        print(f"{row['N']:7d}  {row['lambda']:8.4f}  {row['box_bohr']:9.1f}  {row['box_nm']:8.3f}  {row['replicas']:8d}  {row['atoms']:7d}")

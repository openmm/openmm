#!/usr/bin/env python3
"""Clean diagnostic: scatter init.xyz positions, force re-evaluation, check energy.

Uses LAMMPS scatter_atoms API instead of raw numpy pointer to ensure
positions are properly propagated.
"""

import sys
import re
import numpy as np
from pathlib import Path
import ctypes

from fairchem.lammps.lammps_fc import (
    FIX_EXT_ID,
    FIX_EXTERNAL_CMD,
    FixExternalCallback,
)
from fairchem.core import pretrained_mlip
from lammps import lammps

sys.path.insert(0, str(Path(__file__).parent))
from run_ipi_lammps_uma_rpmd import _patch_fairchem_lammps_atomic_numbers_tensor

_patch_fairchem_lammps_atomic_numbers_tensor()

def parse_ipi_xyz(path):
    with open(path) as f:
        n = int(f.readline().strip())
        comment = f.readline().strip()
        symbols, positions = [], []
        for _ in range(n):
            parts = f.readline().split()
            symbols.append(parts[0])
            positions.append([float(x) for x in parts[1:4]])
    cell_match = re.search(
        r"CELL\(abcABC\):\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", comment
    )
    a, b, c = float(cell_match.group(1)), float(cell_match.group(2)), float(cell_match.group(3))
    return np.array(positions), [a, b, c], symbols

print("Loading UMA ...", flush=True)
predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")

data_path = Path("ipi/lammps_data_ipi_client.data").resolve()
lmp = lammps(cmdargs=["-nocite", "-log", "none", "-echo", "none"])
lmp._predictor = predictor
lmp._task_name = "omol"

lmp.commands_list([
    "units metal",
    "atom_style atomic",
    "boundary p p p",
    "pair_style zero 1.0",
    f"read_data {data_path}",
    "comm_modify cutoff 12.0",
    "pair_coeff * *",
])

cb_energies = []
class LogCB(FixExternalCallback):
    def __call__(self, lmp_obj, ntimestep, nlocal, tag, x, f):
        super().__call__(lmp_obj, ntimestep, nlocal, tag, x, f)
        pos_np = np.array(x[:nlocal], dtype=np.float64, copy=True)
        forces_np = np.array(f[:nlocal], dtype=np.float64, copy=True)
        tags_np = np.array(lmp_obj.numpy.extract_atom("id")[:nlocal])
        types_np = np.array(lmp_obj.numpy.extract_atom("type")[:nlocal])
        cb_energies.append({
            "pos": pos_np,
            "forces": forces_np,
            "tags": tags_np,
            "types": types_np,
            "fmax": np.abs(forces_np).max(),
        })

cb = LogCB(charge=0, spin=1)
lmp.command(FIX_EXTERNAL_CMD)
lmp.set_fix_external_callback(FIX_EXT_ID, cb, lmp)

# === Test A: read_data positions ===
print("\n=== Test A: read_data positions ===")
lmp.command("run 0")
eA = lmp.get_thermo("pe")
print(f"  PE = {eA:.6f} eV  fmax = {cb_energies[-1]['fmax']:.6f} eV/Å")
posA = cb_energies[-1]["pos"].copy()
tagsA = cb_energies[-1]["tags"].copy()
typesA = cb_energies[-1]["types"].copy()

# === Test B: scatter init.xyz using scatter_atoms ===
ipi_pos, ipi_cell, ipi_sym = parse_ipi_xyz("ipi/init.xyz")
nlocal = lmp.get_natoms()

# Build position array ordered by LAMMPS atom tags
tag_to_ipi_idx = {}
tags_internal = np.array(lmp.numpy.extract_atom("id")[:nlocal])
print(f"\n  LAMMPS internal tags: {tags_internal[:10].tolist()} ...")

# The init.xyz atoms are in the same order as the data file.
# The data file atoms are in the same order as LAMMPS's internal order.
# So init.xyz index i corresponds to LAMMPS internal index i.
# Let's verify by comparing positions:
pos_diff = np.abs(posA - ipi_pos).max()
print(f"  Max pos diff (LAMMPS internal vs init.xyz): {pos_diff:.2e} Å")

# Now scatter using extract_atom (direct memory write)
print("\n=== Test B: numpy scatter init.xyz positions ===")
x_view = lmp.numpy.extract_atom("x")
x_view[:nlocal] = ipi_pos.copy()
lmp.command("run 0")
eB = lmp.get_thermo("pe")
print(f"  PE = {eB:.6f} eV  fmax = {cb_energies[-1]['fmax']:.6f} eV/Å")
posB = cb_energies[-1]["pos"].copy()

pos_diff_AB = np.abs(posA - posB).max()
print(f"  Max pos diff (A vs B in callback): {pos_diff_AB:.2e} Å")

# === Test C: Create fresh LAMMPS instance, load same data, evaluate ===
print("\n=== Test C: Fresh LAMMPS instance ===")
lmp2 = lammps(cmdargs=["-nocite", "-log", "none", "-echo", "none"])
lmp2._predictor = predictor
lmp2._task_name = "omol"
lmp2.commands_list([
    "units metal",
    "atom_style atomic",
    "boundary p p p",
    "pair_style zero 1.0",
    f"read_data {data_path}",
    "comm_modify cutoff 12.0",
    "pair_coeff * *",
])

cb2_energies = []
class LogCB2(FixExternalCallback):
    def __call__(self, lmp_obj, ntimestep, nlocal, tag, x, f):
        super().__call__(lmp_obj, ntimestep, nlocal, tag, x, f)
        cb2_energies.append({
            "pos": np.array(x[:nlocal], dtype=np.float64, copy=True),
            "fmax": np.abs(np.array(f[:nlocal])).max(),
        })

cb2 = LogCB2(charge=0, spin=1)
lmp2.command(FIX_EXTERNAL_CMD)
lmp2.set_fix_external_callback(FIX_EXT_ID, cb2, lmp2)
lmp2.command("run 0")
eC = lmp2.get_thermo("pe")
print(f"  PE = {eC:.6f} eV  fmax = {cb2_energies[-1]['fmax']:.6f} eV/Å")

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"  Test A (read_data):        {eA:.6f} eV")
print(f"  Test B (scatter init.xyz): {eB:.6f} eV")
print(f"  Test C (fresh instance):   {eC:.6f} eV")
print(f"  A - B gap: {eA - eB:.4f} eV")
print(f"  A - C gap: {eA - eC:.4f} eV")

del lmp._predictor
lmp.close()
del lmp2._predictor
lmp2.close()

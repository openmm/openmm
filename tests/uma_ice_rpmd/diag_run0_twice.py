#!/usr/bin/env python3
"""Minimal test: run 0 twice on the same LAMMPS state, check if energy changes.

Also tests: explicit position comparison in the callback between successive calls.
"""

import sys
import numpy as np
from pathlib import Path

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

call_log = []

class DetailedCallback(FixExternalCallback):
    def __call__(self, lmp_obj, ntimestep, nlocal, tag, x, f):
        pos_before = np.array(x[:nlocal], dtype=np.float64, copy=True)
        super().__call__(lmp_obj, ntimestep, nlocal, tag, x, f)
        pos_after = np.array(x[:nlocal], dtype=np.float64, copy=True)
        forces = np.array(f[:nlocal], dtype=np.float64, copy=True)

        n_rows_x = len(x)
        types = np.array(lmp_obj.numpy.extract_atom("type")[:nlocal])
        tags = np.array(lmp_obj.numpy.extract_atom("id")[:nlocal])

        call_log.append({
            "ntimestep": ntimestep,
            "nlocal": nlocal,
            "n_rows_x": n_rows_x,
            "pos": pos_before.copy(),
            "pos_after": pos_after.copy(),
            "forces": forces.copy(),
            "fmax": np.abs(forces).max(),
            "types": types.copy(),
            "tags": tags.copy(),
            "pos_x_range": (pos_before[:, 0].min(), pos_before[:, 0].max()),
            "pos_changed": not np.allclose(pos_before, pos_after),
        })

cb = DetailedCallback(charge=0, spin=1)
lmp.command(FIX_EXTERNAL_CMD)
lmp.set_fix_external_callback(FIX_EXT_ID, cb, lmp)

print("\n=== run 0 (first) ===")
lmp.command("run 0")
e1 = lmp.get_thermo("pe")
print(f"  PE = {e1:.6f} eV")

print("\n=== run 0 (second, no position change) ===")
lmp.command("run 0")
e2 = lmp.get_thermo("pe")
print(f"  PE = {e2:.6f} eV")
print(f"  ΔE = {e2 - e1:.6f} eV")

print(f"\nCallback calls: {len(call_log)}")
for i, cl in enumerate(call_log):
    print(f"  cb#{i}: step={cl['ntimestep']}  nlocal={cl['nlocal']}  "
          f"n_rows_x={cl['n_rows_x']}  "
          f"fmax={cl['fmax']:.6f}  "
          f"pos_x=[{cl['pos_x_range'][0]:.4f}, {cl['pos_x_range'][1]:.4f}]  "
          f"pos_changed_by_callback={cl['pos_changed']}")

if len(call_log) >= 2:
    pos_diff = np.abs(call_log[0]["pos"] - call_log[1]["pos"]).max()
    type_match = np.array_equal(call_log[0]["types"], call_log[1]["types"])
    tag_match = np.array_equal(call_log[0]["tags"], call_log[1]["tags"])
    print(f"\n  Position diff between calls: {pos_diff:.2e} Å")
    print(f"  Types match: {type_match}")
    print(f"  Tags match: {tag_match}")
    print(f"\n  cb#0 first 5 tags: {call_log[0]['tags'][:5].tolist()}")
    print(f"  cb#0 first 5 types: {call_log[0]['types'][:5].tolist()}")
    print(f"  cb#0 first 5 pos X: {call_log[0]['pos'][:5, 0].tolist()}")
    print(f"  cb#1 first 5 pos X: {call_log[1]['pos'][:5, 0].tolist()}")

# Now test: wrap_positions effect
from ase.geometry import wrap_positions as ase_wrap
boxlo, boxhi, xy, yz, xz, periodicity, _ = lmp.extract_box()
from fairchem.lammps.lammps_fc import restricted_cell_from_lammps_box
cell = restricted_cell_from_lammps_box(boxlo, boxhi, xy, yz, xz)
cell_np = cell.squeeze().numpy()

if call_log:
    pos0 = call_log[0]["pos"]
    wrapped0 = ase_wrap(pos0, cell=cell_np, pbc=periodicity, eps=0)
    wrap_diff = np.abs(pos0 - wrapped0).max()
    print(f"\n  Wrap diff (cb#0): {wrap_diff:.4f} Å")
    print(f"  Wrapped X range: [{wrapped0[:, 0].min():.4f}, {wrapped0[:, 0].max():.4f}]")
    print(f"  Raw X range:     [{pos0[:, 0].min():.4f}, {pos0[:, 0].max():.4f}]")
    print(f"  Cell: {cell_np}")
    print(f"  Box: lo={boxlo}, hi={boxhi}")

del lmp._predictor
lmp.close()

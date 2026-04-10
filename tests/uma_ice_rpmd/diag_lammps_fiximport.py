#!/usr/bin/env python3
"""Diagnose the 323 eV PE gap: simulate what fix_ipi does to LAMMPS.

Steps:
1. Load the minimized LAMMPS data file (same as i-PI client)
2. Set up fix_external with UMA callback (same as i-PI client)
3. Evaluate energy on the initial positions (from read_data) — should match minimize
4. Manually set positions to the init.xyz coordinates (simulating fix_ipi)
5. Evaluate energy — check if it matches standalone ASE or shows the 323 eV gap
6. Also test with wrapped positions to check PBC sensitivity
"""

import sys
import re
import numpy as np
from pathlib import Path


def parse_ipi_xyz(path: str):
    """Parse i-PI extended XYZ, return positions (Å), cell [a,b,c], symbols."""
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


def main():
    device = sys.argv[1] if len(sys.argv) > 1 else "cuda"

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

    print("Loading UMA predictor ...", flush=True)
    predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device=device)

    data_path = Path("ipi/lammps_data_ipi_client.data").resolve()
    print(f"Using LAMMPS data: {data_path}", flush=True)

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

    energies_captured = []

    class LoggingCallback(FixExternalCallback):
        def __call__(self, lmp_obj, ntimestep, nlocal, tag, x, f):
            super().__call__(lmp_obj, ntimestep, nlocal, tag, x, f)
            boxlo, boxhi, *_ = lmp_obj.extract_box()
            pos = np.array(x[:nlocal], dtype=np.float64)
            forces = np.array(f[:nlocal], dtype=np.float64)
            e = lmp_obj.get_thermo("pe")
            energies_captured.append({
                "step": ntimestep,
                "pe": e,
                "pos_x_range": (pos[:, 0].min(), pos[:, 0].max()),
                "pos_y_range": (pos[:, 1].min(), pos[:, 1].max()),
                "pos_z_range": (pos[:, 2].min(), pos[:, 2].max()),
                "box": (boxlo, boxhi),
                "fmax": np.abs(forces).max(),
            })

    cb = LoggingCallback(charge=0, spin=1)
    lmp.command(FIX_EXTERNAL_CMD)
    lmp.set_fix_external_callback(FIX_EXT_ID, cb, lmp)

    # --- Test 1: Energy on the read_data positions (minimized structure) ---
    print("\n=== Test 1: read_data positions (minimized, symmetric box) ===")
    lmp.command("run 0")
    if energies_captured:
        e1 = energies_captured[-1]
        print(f"  PE = {e1['pe']:.4f} eV")
        print(f"  Pos X: [{e1['pos_x_range'][0]:.4f}, {e1['pos_x_range'][1]:.4f}]")
        print(f"  Pos Y: [{e1['pos_y_range'][0]:.4f}, {e1['pos_y_range'][1]:.4f}]")
        print(f"  Pos Z: [{e1['pos_z_range'][0]:.4f}, {e1['pos_z_range'][1]:.4f}]")
        print(f"  Box:   [{e1['box'][0]}, {e1['box'][1]}]")
        print(f"  Fmax:  {e1['fmax']:.6f} eV/Å")

    # --- Test 2: Manually set positions to init.xyz coords (simulating fix_ipi) ---
    print("\n=== Test 2: Scatter init.xyz positions (simulating fix_ipi) ===")
    ipi_pos, ipi_cell, ipi_symbols = parse_ipi_xyz("ipi/init.xyz")
    nlocal = lmp.get_natoms()
    assert nlocal == len(ipi_pos), f"Atom count mismatch: {nlocal} vs {len(ipi_pos)}"

    x_ptr = lmp.numpy.extract_atom("x")
    x_ptr[:nlocal] = ipi_pos
    lmp.command("run 0")
    if len(energies_captured) > 1:
        e2 = energies_captured[-1]
        print(f"  PE = {e2['pe']:.4f} eV")
        print(f"  Pos X: [{e2['pos_x_range'][0]:.4f}, {e2['pos_x_range'][1]:.4f}]")
        print(f"  Pos Y: [{e2['pos_y_range'][0]:.4f}, {e2['pos_y_range'][1]:.4f}]")
        print(f"  Pos Z: [{e2['pos_z_range'][0]:.4f}, {e2['pos_z_range'][1]:.4f}]")
        print(f"  Box:   [{e2['box'][0]}, {e2['box'][1]}]")
        print(f"  Fmax:  {e2['fmax']:.6f} eV/Å")

    # --- Test 3: Same positions but reset box to fix_ipi convention (centered) ---
    print("\n=== Test 3: init.xyz positions + reset box to [-L/2, L/2] ===")
    a, b, c = ipi_cell
    lmp.command(
        f"change_box all x final {-0.5*a:.8f} {0.5*a:.8f} "
        f"y final {-0.5*b:.8f} {0.5*b:.8f} "
        f"z final {-0.5*c:.8f} {0.5*c:.8f}"
    )
    x_ptr[:nlocal] = ipi_pos
    lmp.command("run 0")
    if len(energies_captured) > 2:
        e3 = energies_captured[-1]
        print(f"  PE = {e3['pe']:.4f} eV")
        print(f"  Pos X: [{e3['pos_x_range'][0]:.4f}, {e3['pos_x_range'][1]:.4f}]")
        print(f"  Pos Y: [{e3['pos_y_range'][0]:.4f}, {e3['pos_y_range'][1]:.4f}]")
        print(f"  Pos Z: [{e3['pos_z_range'][0]:.4f}, {e3['pos_z_range'][1]:.4f}]")
        print(f"  Box:   [{e3['box'][0]}, {e3['box'][1]}]")
        print(f"  Fmax:  {e3['fmax']:.6f} eV/Å")

    # --- Test 4: Bohr round-trip (Å -> Bohr -> Å, same as fix_ipi) ---
    print("\n=== Test 4: Bohr round-trip positions ===")
    BOHR_TO_ANG = 0.52917721067
    pos_bohr = ipi_pos / BOHR_TO_ANG
    pos_roundtrip = pos_bohr * BOHR_TO_ANG
    max_diff = np.abs(ipi_pos - pos_roundtrip).max()
    print(f"  Max position diff after Å->Bohr->Å: {max_diff:.2e} Å")
    x_ptr[:nlocal] = pos_roundtrip
    lmp.command("run 0")
    if len(energies_captured) > 3:
        e4 = energies_captured[-1]
        print(f"  PE = {e4['pe']:.4f} eV")
        print(f"  Fmax:  {e4['fmax']:.6f} eV/Å")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for i, ec in enumerate(energies_captured):
        print(f"  Test {i+1}: PE = {ec['pe']:.4f} eV  fmax = {ec['fmax']:.6f} eV/Å")
    print()
    print("Reference: standalone ASE = -133146.2510 eV")
    print("Reference: i-PI step-0    = -132823 eV (the mystery)")
    print("If all tests show ~-133146 eV, the PE gap is from i-PI's internal handling, not LAMMPS.")

    del lmp._predictor
    lmp.close()


if __name__ == "__main__":
    main()

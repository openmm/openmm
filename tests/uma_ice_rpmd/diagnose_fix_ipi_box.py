#!/usr/bin/env python3
"""
Diagnose whether LAMMPS ``fix ipi``-style box re-centering changes UMA energy.

``fix_ipi`` maps the i-PI cell to a symmetric box with ``boxlo = -L/2``, ``boxhi = +L/2``
and PBC-wraps Cartesian positions. FairChem's ``FixExternalCallback`` uses
``ase.geometry.wrap_positions`` with the LAMMPS cell matrix (lengths only). This script:

1. Loads ``data.ice_uma_64``, registers UMA ``fix external`` (same pattern as
   ``run_ipi_lammps_uma_rpmd._lammps_uma_minimize_to_ipi_init``).
2. Runs ``minimize`` to the same tolerances as the i-PI driver (optional quick mode).
3. Records ``etotal`` (metal units, eV) and the potential from ``thermo_pe``.
4. Applies ``displace_atoms`` + ``change_box`` to mimic ``fix_ipi`` symmetric box.
5. Runs ``run 0`` again and compares energies.

Empirically (UMA ``uma-s-1p1``, ice 64 mol), |ΔE| is ~1e-5 eV: **box re-centering is
not** what caused the ~323 eV i-PI vs LAMMPS-minimize gap. The driver still symmetrizes
the minimized ``write_data`` and centers ``init.xyz`` so the client ``read_data``
geometry matches i-PI before the first socket exchange (see ``run_ipi_lammps_uma_rpmd``).

Usage (from ``tests/uma_ice_rpmd``, conda env with fairchem + lammps + CUDA):

  python diagnose_fix_ipi_box.py
  python diagnose_fix_ipi_box.py --quick
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_ROOT = _SCRIPT_DIR

if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))


def main() -> int:
    ap = argparse.ArgumentParser(description="Diagnose UMA energy vs LAMMPS box centering")
    ap.add_argument(
        "--data",
        type=Path,
        default=_ROOT / "lammps" / "data.ice_uma_64",
        help="LAMMPS data file",
    )
    ap.add_argument(
        "--quick",
        action="store_true",
        help="Short minimize (100 iter) instead of production caps",
    )
    ap.add_argument(
        "--device",
        default=os.environ.get("UMA_DEVICE", "cuda"),
        help="Torch device for UMA (default cuda or UMA_DEVICE)",
    )
    ap.add_argument(
        "--task-name",
        default=os.environ.get("UMA_TASK_NAME", "omol"),
        help="FairChem dataset tag for UMA (must match model mapping, default omol)",
    )
    args = ap.parse_args()

    if not args.data.is_file():
        print(f"Missing data file: {args.data}", file=sys.stderr)
        return 2

    try:
        from fairchem.lammps.lammps_fc import (
            FIX_EXT_ID,
            FIX_EXTERNAL_CMD,
            FixExternalCallback,
        )
        from lammps import lammps as Lammps
    except ImportError as exc:
        print(f"Requires fairchem-lammps + lammps: {exc}", file=sys.stderr)
        return 3

    try:
        from fairchem.core import pretrained_mlip
    except ImportError:
        from fairchem.core.calculate import pretrained_mlip

    from run_ipi_lammps_uma_rpmd import _patch_fairchem_lammps_atomic_numbers_tensor

    _patch_fairchem_lammps_atomic_numbers_tensor()

    model = os.environ.get("UMA_MODEL", "uma-s-1p1")
    task_name = args.task_name
    print(f"Loading UMA ({model}) on {args.device} ...", flush=True)
    predictor = pretrained_mlip.get_predict_unit(model, device=args.device)

    log_path = _ROOT / "pipeline_out" / "diagnose_fix_ipi_box.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    machine = os.environ.get("LAMMPS_MACHINE_NAME")

    lmp = Lammps(name=machine, cmdargs=["-nocite", "-log", str(log_path), "-echo", "screen"])
    lmp._predictor = predictor
    lmp._task_name = task_name

    try:
        lmp.commands_list(
            [
                "units metal",
                "atom_style atomic",
                "boundary p p p",
                f"read_data {args.data.resolve()}",
                "comm_modify cutoff 12.0",
                "pair_style zero 1.0",
                "pair_coeff * *",
            ]
        )
        lmp.command(FIX_EXTERNAL_CMD)
        lmp.set_fix_external_callback(FIX_EXT_ID, FixExternalCallback(charge=0, spin=1), lmp)

        if args.quick:
            lmp.command("minimize 0.0 1.0e-4 100 300")
        else:
            lmp.command("minimize 0.0 1.0e-6 2000 50000")

        lmp.command("run 0")
        e0 = lmp.get_thermo("etotal")
        pe0 = lmp.get_thermo("pe")

        boxlo, boxhi, xy, yz, xz, *_ = lmp.extract_box()
        lx = boxhi[0] - boxlo[0]
        ly = boxhi[1] - boxlo[1]
        lz = boxhi[2] - boxlo[2]
        cx = 0.5 * (boxlo[0] + boxhi[0])
        cy = 0.5 * (boxlo[1] + boxhi[1])
        cz = 0.5 * (boxlo[2] + boxhi[2])

        print("\n--- After minimize (original LAMMPS box) ---", flush=True)
        print(f"  boxlo = [{boxlo[0]:.6f}, {boxlo[1]:.6f}, {boxlo[2]:.6f}]", flush=True)
        print(f"  boxhi = [{boxhi[0]:.6f}, {boxhi[1]:.6f}, {boxhi[2]:.6f}]", flush=True)
        print(f"  etotal = {e0:.8f} eV", flush=True)
        print(f"  pe     = {pe0:.8f} eV", flush=True)

        # Mimic fix_ipi: shift atoms so origin-centered box, then reset box limits
        lmp.command(f"displace_atoms all move {-cx:.16f} {-cy:.16f} {-cz:.16f}")
        lmp.command(
            f"change_box all x final {-0.5 * lx:.16f} {0.5 * lx:.16f} "
            f"y final {-0.5 * ly:.16f} {0.5 * ly:.16f} "
            f"z final {-0.5 * lz:.16f} {0.5 * lz:.16f}"
        )
        lmp.command("run 0")
        e1 = lmp.get_thermo("etotal")
        pe1 = lmp.get_thermo("pe")

        boxlo2, boxhi2, *_ = lmp.extract_box()
        print("\n--- After displace_atoms + change_box (fix_ipi-style) ---", flush=True)
        print(f"  boxlo = [{boxlo2[0]:.6f}, {boxlo2[1]:.6f}, {boxlo2[2]:.6f}]", flush=True)
        print(f"  boxhi = [{boxhi2[0]:.6f}, {boxhi2[1]:.6f}, {boxhi2[2]:.6f}]", flush=True)
        print(f"  etotal = {e1:.8f} eV", flush=True)
        print(f"  pe     = {pe1:.8f} eV", flush=True)

        de = e1 - e0
        dpe = pe1 - pe0
        print("\n=== Summary ===", flush=True)
        print(f"  Δ etotal = {de:.6e} eV", flush=True)
        print(f"  Δ pe     = {dpe:.6e} eV", flush=True)
        print(f"  ASCII_DELTA_etotal_eV = {de:.6e}", flush=True)
        if abs(de) < 1e-3 and abs(dpe) < 1e-3:
            print(
                "  Interpretation: energy invariant to box centering (UMA+ASE consistent).",
                flush=True,
            )
        elif abs(de) < 0.1:
            print(
                "  Interpretation: small numerical drift only; unlikely to explain ~323 eV i-PI gap.",
                flush=True,
            )
        else:
            print(
                "  Interpretation: significant sensitivity to box origin / wrap path — "
                "pre-center init.xyz for i-PI (see convert_lammps_to_ipi_xyz).",
                flush=True,
            )
    finally:
        del lmp._predictor
        lmp.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

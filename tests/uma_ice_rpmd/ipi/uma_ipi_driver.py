"""i-PI PES driver for FairChem UMA models via ASE calculator interface.

Follows the same pattern as ``ipi.pes.mace``: subclass ``ASEDriver`` and
set ``self.ase_calculator`` to a ``FAIRChemCalculator``.

Command-line (via ``i-pi-py_driver -m custom``):

    i-pi-py_driver -P uma_ipi_driver.py -m custom \
        -a 127.0.0.1 -p 65535 \
        -o template=init.xyz,model=uma-s-1p1,device=cuda,task_name=omol

Parameters
----------
template : str
    Path to an ASE-readable structure file (xyz with i-PI CELL comment).
model : str, default ``"uma-s-1p1"``
    FairChem pretrained model name passed to ``get_predict_unit``.
device : str, default ``"cuda"``
    PyTorch device for inference.
task_name : str, default ``"omol"``
    UMA task name (``omol``, ``omat``, ``oc20``, …).

Notes
-----
i-PI writes xyz files with a ``# CELL(abcABC): a b c α β γ`` comment line.
ASE's generic xyz reader does not parse this format and leaves the cell as
zeros, which FAIRChem rejects with ``AllZeroUnitCellError``.  We parse the
cell explicitly in ``check_parameters`` and set it on the template atoms.
"""

import re
import numpy as np
from ipi.pes.ase import ASEDriver

__DRIVER_NAME__ = "uma"
__DRIVER_CLASS__ = "UMADriver"

# Pattern: # CELL(abcABC):  a  b  c  alpha  beta  gamma
_CELL_RE = re.compile(
    r"CELL\(abcABC\):\s+"
    r"([\d.eE+\-]+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)"
    r"\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)"
)


def _read_ipi_cell_angstrom(xyz_path: str) -> np.ndarray | None:
    """Return a (3,3) cell matrix in Angstrom from an i-PI xyz comment, or None."""
    with open(xyz_path) as fh:
        fh.readline()  # natoms line
        comment = fh.readline()
    m = _CELL_RE.search(comment)
    if m is None:
        return None
    a, b, c = float(m.group(1)), float(m.group(2)), float(m.group(3))
    alpha, beta, gamma = (
        np.deg2rad(float(m.group(4))),
        np.deg2rad(float(m.group(5))),
        np.deg2rad(float(m.group(6))),
    )
    # Build cell from lattice parameters (assumes orthorhombic when α=β=γ=90°)
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = np.sqrt(max(c**2 - cx**2 - cy**2, 0.0))
    cell = np.array([
        [a, 0.0, 0.0],
        [b * np.cos(gamma), b * np.sin(gamma), 0.0],
        [cx, cy, cz],
    ])
    return cell


class UMADriver(ASEDriver):
    """i-PI driver that evaluates forces with a FairChem UMA calculator.

    Overrides ``compute_structure`` to fix two issues in ``ASEDriver``:
    1. ``has_stress=False`` passes ``np.zeros(9)`` to the stress-shape check,
       which raises ValueError because (9,) is neither (6,) nor (3,3).
    2. i-PI's CELL(abcABC) comment is not parsed by ASE's xyz reader, so the
       cell must be read manually and set on the template atoms.
    """

    def __init__(
        self,
        template,
        model="uma-s-1p1",
        device="cuda",
        task_name="omol",
        *args,
        **kwargs,
    ):
        self.model = model
        self.device = device
        self.task_name = task_name
        # Request energy+forces only; we supply a zero virial ourselves.
        super().__init__(template, has_stress=False, *args, **kwargs)

    def check_parameters(self):
        """Load template structure, fix cell, and attach FAIRChemCalculator."""
        try:
            from fairchem.core import FAIRChemCalculator
            from fairchem.core.calculate import pretrained_mlip
        except ImportError as err:
            msg = str(err).lower()
            hint = (
                "FairChem loads Ray Serve, which needs a working pydantic / pydantic-core pair. "
                "If you see pydantic_core or validate_core_schema errors, repair the stack, e.g.\n"
                "  pip install --ignore-installed 'pydantic-core==2.41.5' 'pydantic==2.12.5'\n"
                "(adjust versions to match `pip show pydantic` / your Ray install), or recreate the conda env."
            )
            if "pydantic" in msg or "validate_core" in msg:
                raise RuntimeError(f"{err}\n{hint}") from err
            raise

        super().check_parameters()

        # ASE's xyz reader does not parse i-PI's CELL(abcABC) comment; fix it.
        cell = _read_ipi_cell_angstrom(self.template)
        if cell is not None:
            self.template_ase.set_cell(cell)
            print(f"[UMADriver] Cell set from i-PI comment: {np.diag(cell)} Å", flush=True)
        else:
            print("[UMADriver] WARNING: could not parse cell from template xyz", flush=True)

        self.template_ase.pbc = [True, True, True]
        self.template_ase.info["charge"] = 0
        self.template_ase.info["spin"] = 1

        print(f"[UMADriver] Loading model={self.model} on {self.device} ...", flush=True)
        predict_unit = pretrained_mlip.get_predict_unit(
            self.model, device=self.device
        )
        self.ase_calculator = FAIRChemCalculator(
            predict_unit, task_name=self.task_name
        )
        try:
            _df = bool(predict_unit.direct_forces)
        except Exception as exc:
            print(
                f"[UMADriver] Could not read predict_unit.direct_forces: {exc}",
                flush=True,
            )
        else:
            _desc = (
                "direct force head (predict uses torch.no_grad)"
                if _df
                else "autograd from model energy"
            )
            print(
                f"[UMADriver] predict_unit.direct_forces={_df} — {_desc}",
                flush=True,
            )
        print(f"[UMADriver] FAIRChemCalculator ready (task={self.task_name})", flush=True)

    def compute_structure(self, cell, pos):
        """Evaluate energy and forces; return zero virial (no stress requested).

        Overrides ``ASEDriver.compute_structure`` to avoid the ``ValueError``
        raised when ``has_stress=False`` leaves ``stress = np.zeros(9)``
        (shape ``(9,)``) which fails the ``(3,3)`` symmetry check.
        """
        from ipi.utils.units import unit_to_internal, unit_to_user

        # i-PI uses atomic units; convert to Angstrom for ASE
        pos_ang = unit_to_user("length", "angstrom", pos)
        cell_ang = unit_to_user("length", "angstrom", cell.T)

        structure = self.template_ase.copy()
        structure.positions[:] = pos_ang
        structure.cell[:] = cell_ang
        structure.pbc = [True, True, True]
        structure.info["charge"] = 0
        structure.info["spin"] = 1
        structure.calc = self.ase_calculator

        properties = structure.get_properties(["energy", "forces"])

        pot_ipi = np.asarray(
            unit_to_internal("energy", "electronvolt", properties["energy"]),
            np.float64,
        )
        force_ipi = np.asarray(
            unit_to_internal("force", "ev/ang", properties["forces"]),
            np.float64,
        )
        vir_ipi = np.zeros((3, 3), dtype=np.float64)

        return pot_ipi, force_ipi, vir_ipi, ""

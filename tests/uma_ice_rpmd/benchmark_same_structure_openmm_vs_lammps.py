"""Legacy import path; implementation in ``ice_rpmd_workflow.benchmark_same_structure_openmm_vs_lammps``."""
from __future__ import annotations

import importlib
import sys

_sub = importlib.import_module(
    "ice_rpmd_workflow.benchmark_same_structure_openmm_vs_lammps"
)
sys.modules[__name__] = _sub

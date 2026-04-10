"""Legacy import path; implementation in ``ice_rpmd_workflow.rpmd_thermo_utils``."""
from __future__ import annotations

import importlib
import sys

_sub = importlib.import_module("ice_rpmd_workflow.rpmd_thermo_utils")
sys.modules[__name__] = _sub

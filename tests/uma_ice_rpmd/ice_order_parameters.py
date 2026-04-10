"""Legacy import path; implementation in ``ice_rpmd_workflow.ice_order_parameters``."""
from __future__ import annotations

import importlib
import sys

_sub = importlib.import_module("ice_rpmd_workflow.ice_order_parameters")
sys.modules[__name__] = _sub

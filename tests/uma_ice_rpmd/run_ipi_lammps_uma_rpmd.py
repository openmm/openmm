#!/usr/bin/env python3
"""Shim: implementation in ``scripts/run_ipi_lammps_uma_rpmd`` (same module API + CLI)."""
from __future__ import annotations

import importlib
import sys

_mod = importlib.import_module("scripts.run_ipi_lammps_uma_rpmd")
sys.modules[__name__] = _mod

if __name__ == "__main__":
    _mod.main()

"""CUDA tensor bridge for OpenMM contexts.

This package exposes a small wrapper around a C++ extension that copies
positions and forces between PyTorch CUDA tensors and OpenMM CUDA buffers
without staging through host memory.
"""

from __future__ import annotations

import ctypes
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import torch  # noqa: E402, F401  # preload libtorch shared libraries


def _candidate_openmm_libs() -> list[Path]:
    candidates: list[Path] = []
    for env_name in ("OPENMM_LIB_DIR", "OPENMM_CUDA_LIB_DIR"):
        value = os.environ.get(env_name)
        if value:
            candidates.append(Path(value) / "libOpenMM.so")
    openmm_dir = os.environ.get("OPENMM_DIR")
    if openmm_dir:
        candidates.append(Path(openmm_dir) / "build" / "libOpenMM.so")
        candidates.append(Path(openmm_dir) / "lib" / "libOpenMM.so")
    # Editable source-tree fallback: .../openmm/python/openmm_cuda_bridge/openmm_cuda_bridge/__init__.py
    repo_openmm = Path(__file__).resolve().parents[3]
    candidates.append(repo_openmm / "build" / "libOpenMM.so")
    candidates.append(repo_openmm / "lib" / "libOpenMM.so")
    return candidates


def _preload_matching_openmm() -> None:
    for candidate in _candidate_openmm_libs():
        if candidate.exists():
            ctypes.CDLL(os.fspath(candidate), mode=ctypes.RTLD_GLOBAL)
            return


_preload_matching_openmm()

from ._openmm_cuda_bridge import CudaBridge as _CudaBridge


def _context_pointer(context: Any) -> int:
    """Return the raw C++ pointer from an OpenMM SWIG-wrapped Context."""
    swig_ptr = getattr(context, "this", None)
    if swig_ptr is None:
        raise TypeError("Expected an OpenMM SWIG Context object with a .this pointer")
    return int(swig_ptr)


@dataclass
class CudaBridge:
    """Python wrapper that accepts a normal OpenMM ``Context`` object."""

    context: Any
    groups: int = 0xFFFFFFFF

    def __post_init__(self) -> None:
        self._bridge = _CudaBridge(_context_pointer(self.context), int(self.groups))

    @property
    def num_atoms(self) -> int:
        return int(self._bridge.num_atoms())

    @property
    def padded_num_atoms(self) -> int:
        return int(self._bridge.padded_num_atoms())

    @property
    def device_index(self) -> int:
        return int(self._bridge.device_index())

    def set_positions(self, positions_angstrom):
        self._bridge.set_positions(positions_angstrom)

    def compute_forces(self):
        return self._bridge.compute_forces()

    def evaluate(self, positions_angstrom):
        return self._bridge.evaluate(positions_angstrom)


__all__ = ["CudaBridge"]

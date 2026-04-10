#!/usr/bin/env python3
"""Minimal OpenMM CUDA Context (no ML). Diagnoses driver vs OpenMM CUDA plugin issues.

OpenMM raises a generic "No compatible CUDA device is available" when *every*
cuCtxCreate attempt fails (or no GPUs are enumerated). That is not always PTX 222.

Usage:
  OPENMM_PLUGIN_DIR=$CONDA_PREFIX/lib/plugins python scripts/smoke_openmm_cuda.py

Exit 0 on success, 1 on failure (prints exception).
"""
from __future__ import annotations

import ctypes
import glob
import os
import subprocess
import sys

# https://docs.nvidia.com/cuda/cuda-driver-api/group__CUDA__TYPES.html
_CU_ERRORS: dict[int, str] = {
    0: "CUDA_SUCCESS",
    1: "CUDA_ERROR_INVALID_VALUE",
    2: "CUDA_ERROR_OUT_OF_MEMORY",
    3: "CUDA_ERROR_NOT_INITIALIZED",
    100: "CUDA_ERROR_NO_DEVICE",
    101: "CUDA_ERROR_INVALID_DEVICE",
    200: "CUDA_ERROR_INVALID_IMAGE",
    201: "CUDA_ERROR_INVALID_CONTEXT",
    208: "CUDA_ERROR_CONTEXT_ALREADY_CURRENT",
    209: "CUDA_ERROR_MAP_FAILED",
    212: "CUDA_ERROR_LAUNCH_FAILED",
    214: "CUDA_ERROR_PRIMARY_CONTEXT_ACTIVE",
    219: "CUDA_ERROR_UNKNOWN",
    222: "CUDA_ERROR_UNSUPPORTED_PTX_VERSION",
    999: "CUDA_ERROR_COMPAT_NOT_SUPPORTED_ON_DEVICE",
}


def _cu_str(code: int) -> str:
    return _CU_ERRORS.get(code, f"CUresult {code}")


def _probe_driver_api() -> None:
    """Print libcuda device count and first-device cuCtxCreate result (OpenMM swallows these)."""
    vis = os.environ.get("CUDA_VISIBLE_DEVICES", "<unset>")
    print(f"smoke_openmm_cuda: CUDA_VISIBLE_DEVICES={vis}", file=sys.stderr)
    try:
        lib = ctypes.CDLL("libcuda.so.1")
    except OSError as e:
        print(f"smoke_openmm_cuda: libcuda.so.1 not loaded: {e}", file=sys.stderr)
        return

    cu_init = lib.cuInit
    cu_init.argtypes = [ctypes.c_uint]
    cu_init.restype = ctypes.c_int

    n = ctypes.c_int(0)
    cu_device_get_count = lib.cuDeviceGetCount
    cu_device_get_count.argtypes = [ctypes.POINTER(ctypes.c_int)]
    cu_device_get_count.restype = ctypes.c_int

    r = cu_init(0)
    if r != 0:
        print(f"smoke_openmm_cuda: cuInit failed: {_cu_str(r)}", file=sys.stderr)
        return

    r = cu_device_get_count(ctypes.byref(n))
    if r != 0:
        print(f"smoke_openmm_cuda: cuDeviceGetCount failed: {_cu_str(r)}", file=sys.stderr)
        return

    print(f"smoke_openmm_cuda: libcuda reports {n.value} device(s)", file=sys.stderr)
    if n.value < 1:
        print(
            "  → No GPUs visible to the driver (wrong CUDA_VISIBLE_DEVICES, no GPU, "
            "or driver/module issue).",
            file=sys.stderr,
        )
        return

    cu_device_get = lib.cuDeviceGet
    cu_device_get.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.c_int]
    cu_device_get.restype = ctypes.c_int

    dev = ctypes.c_int(0)
    r = cu_device_get(ctypes.byref(dev), 0)
    if r != 0:
        print(f"smoke_openmm_cuda: cuDeviceGet(0) failed: {_cu_str(r)}", file=sys.stderr)
        return

    ctx = ctypes.c_void_p(0)
    # Match OpenMM: CU_CTX_MAP_HOST | CU_CTX_SCHED_SPIN (cuda.h).
    flags = 0x08 | 0x01
    # Prefer legacy 3-arg API (matches OpenMM when CUDA_VERSION < 13000).
    for sym in ("cuCtxCreate_v2", "cuCtxCreate"):
        if not hasattr(lib, sym):
            continue
        fn = getattr(lib, sym)
        fn.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.c_uint, ctypes.c_int]
        fn.restype = ctypes.c_int
        r = fn(ctypes.byref(ctx), flags, dev.value)
        print(
            f"smoke_openmm_cuda: {sym}(device 0) → {_cu_str(r)}",
            file=sys.stderr,
        )
        if r == 2:
            print(
                "  → GPU memory exhausted: quit other CUDA jobs (nvidia-smi), "
                "or free VRAM; OpenMM then reports a generic 'No compatible CUDA device'.",
                file=sys.stderr,
            )
        if r == 0 and ctx.value:
            cu_ctx_destroy = lib.cuCtxDestroy_v2
            if not hasattr(lib, "cuCtxDestroy_v2"):
                cu_ctx_destroy = lib.cuCtxDestroy
            cu_ctx_destroy.argtypes = [ctypes.c_void_p]
            cu_ctx_destroy.restype = ctypes.c_int
            cu_ctx_destroy(ctx)
        break
    else:
        print("smoke_openmm_cuda: no cuCtxCreate symbol in libcuda", file=sys.stderr)


def _probe_torch_cuda() -> None:
    try:
        import torch
    except ImportError:
        print("smoke_openmm_cuda: torch not installed (skip torch probe)", file=sys.stderr)
        return
    ok = bool(torch.cuda.is_available())
    print(f"smoke_openmm_cuda: torch.cuda.is_available()={ok}", file=sys.stderr)
    if ok:
        try:
            torch.zeros(1, device="cuda")
            print("smoke_openmm_cuda: torch allocated tensor on cuda:0 OK", file=sys.stderr)
        except Exception as e:
            print(f"smoke_openmm_cuda: torch cuda tensor failed: {e!r}", file=sys.stderr)


def _nvidia_driver_major() -> int | None:
    try:
        out = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=driver_version",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
        if out.returncode != 0 or not out.stdout.strip():
            return None
        line = out.stdout.strip().splitlines()[0].strip()
        if not line or not line[0].isdigit():
            return None
        return int(line.split(".", 1)[0])
    except (ValueError, IndexError, FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return None


def _legacy_cuda_workarounds() -> bool:
    """True → LD_PRELOAD system libnvrtc (driver <550 or OPENMM_FORCE_LEGACY_CUDA=1)."""
    if os.environ.get("OPENMM_SKIP_NVRTC_LD_PREPEND") == "1":
        return False
    if os.environ.get("OPENMM_FORCE_LEGACY_CUDA") == "1":
        return True
    if os.environ.get("OPENMM_FORCE_LEGACY_CUDA") == "0":
        return False
    thr = int(os.environ.get("OPENMM_LEGACY_DRIVER_MAJOR_MIN", "550"))
    maj = _nvidia_driver_major()
    if maj is None:
        return False
    return maj < thr


def _ensure_nvrtc_ld_preload() -> None:
    """On legacy drivers, LD_PRELOAD toolkit libnvrtc (conda python RPATH beats LD_LIBRARY_PATH)."""
    if not _legacy_cuda_workarounds():
        return
    candidates = [
        "/usr/local/cuda-12.0/lib64",
        "/usr/local/cuda-12.1/lib64",
        "/usr/local/cuda-12.2/lib64",
    ]
    env_root = os.environ.get("CUDAToolkit_ROOT", "").strip()
    if env_root:
        candidates.append(os.path.join(env_root, "lib64"))
    for lib64 in candidates:
        nvrtc = ""
        for name in ("libnvrtc.so.12", "libnvrtc.so.12.0", "libnvrtc.so.12.1", "libnvrtc.so.12.2"):
            p = os.path.join(lib64, name)
            if os.path.isfile(p):
                nvrtc = p
                break
        if not nvrtc:
            hits = glob.glob(os.path.join(lib64, "libnvrtc.so.12*"))
            if not hits:
                continue
            nvrtc = hits[0]
        prev = os.environ.get("LD_PRELOAD", "")
        parts = [p for p in prev.split(":") if p]
        if nvrtc not in parts:
            os.environ["LD_PRELOAD"] = nvrtc + (":" + prev if prev else "")
            print(f"smoke_openmm_cuda: set LD_PRELOAD={nvrtc} (NVRTC; beats conda python RPATH)", file=sys.stderr)
        break


def _ensure_plugin_dir() -> None:
    prefix = os.environ.get("CONDA_PREFIX") or os.environ.get("OPENMM_LIB_PATH", "")
    if prefix and not os.environ.get("OPENMM_PLUGIN_DIR"):
        if os.path.isdir(os.path.join(prefix, "lib", "plugins")):
            os.environ["OPENMM_PLUGIN_DIR"] = os.path.join(prefix, "lib", "plugins")
        elif os.path.basename(prefix) == "lib" and os.path.isdir(os.path.join(prefix, "plugins")):
            os.environ["OPENMM_PLUGIN_DIR"] = os.path.join(prefix, "plugins")


def main() -> int:
    _ensure_nvrtc_ld_preload()
    _ensure_plugin_dir()
    try:
        from openmm import Context, Platform, System, VerletIntegrator, Vec3, unit
    except ImportError as e:
        print(f"smoke_openmm_cuda: import openmm failed: {e}", file=sys.stderr)
        return 1

    try:
        platform = Platform.getPlatformByName("CUDA")
    except Exception as e:
        print(f"smoke_openmm_cuda: CUDA platform unavailable: {e}", file=sys.stderr)
        return 1

    system = System()
    system.addParticle(1.0)  # mass in amu
    integrator = VerletIntegrator(1.0 * unit.femtoseconds)
    try:
        context = Context(system, integrator, platform)
    except Exception as e:
        print(f"smoke_openmm_cuda: Context(CUDA) failed: {e!r}", file=sys.stderr)
        print("smoke_openmm_cuda: --- driver API probe (OpenMM hides cuCtxCreate errors) ---", file=sys.stderr)
        _probe_driver_api()
        _probe_torch_cuda()
        msg = str(e)
        if "UNSUPPORTED_PTX_VERSION" in msg or "222" in msg:
            print(
                "  → PTX 222: upgrade NVIDIA driver to 550+ (e.g. 580), or set "
                "OPENMM_FORCE_LEGACY_CUDA=1 and ensure /usr/local/cuda-12.0/lib64/libnvrtc.so.12 exists "
                "(this script LD_PRELOADs it only when driver major < 550 or force is set).",
                file=sys.stderr,
            )
        else:
            print(
                "  Hints: If libcuda shows 0 devices, fix visibility/driver. "
                "If cuCtxCreate fails with 222 → PTX newer than driver can JIT. "
                "If torch works but OpenMM fails with a different error → rebuild OpenMM CUDA "
                "(see install_openmm_fairchem_base.sh).",
                file=sys.stderr,
            )
        return 1

    context.setPositions([Vec3(0, 0, 0) * unit.nanometer])
    state = context.getState(getEnergy=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"smoke_openmm_cuda: OK (CUDA Context created, PE={pe:.6g} kJ/mol)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

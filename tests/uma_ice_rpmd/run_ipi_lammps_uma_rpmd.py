#!/usr/bin/env python3
"""
Orchestrate i-PI RPMD with a force-evaluation client (LAMMPS or Python driver).

Two-process architecture that mirrors i-PI's own design:

  1. **i-PI server** — integrates ring-polymer dynamics, owns positions/momenta,
     listens on a TCP socket (``ffsocket``).
  2. **Force client** — connects to that socket and returns energy, forces, and
     virial for each bead configuration i-PI sends.

The ``--client`` flag selects which force backend evaluates UMA:

  * ``lammps`` (default): LAMMPS Python API with ``fix ipi`` (socket client) and
    ``fix external pf/callback`` (FairChem UMA forces).  This is the standard
    i-PI + LAMMPS workflow.
  * ``python``: ``i-pi-py_driver`` with a custom ASE-based PES
    (``ipi/uma_ipi_driver.py``).  Fewer dependencies (no LAMMPS build) but
    requires ``pip install ipi fairchem-core``.

Prerequisites:
  LAMMPS client:  pip install fairchem-lammps fairchem-core lammps
  Python client:  pip install ipi fairchem-core

Usage:
  cd tests/uma_ice_rpmd
  python run_ipi_lammps_uma_rpmd.py --steps 1000 --device cuda
  python run_ipi_lammps_uma_rpmd.py --molecules 64 --beads 32 --dt-fs 0.1 --steps 10000
  python run_ipi_lammps_uma_rpmd.py --client python --steps 100   # ASE driver instead

i-PI progress is written to ``--ipi-log`` (default ``ipi/i-pi_run.log``). Each launch
backs up the previous log to ``i-pi_run_backup_YYYYMMDD_HHMMSS.log`` if it exists,
then truncates the main log. After a successful client run, the driver prints the
last ``Step:`` from ``ice__i-pi.traj_00.xyz`` so you can confirm length vs ``--prod``.

Optional ``--order-csv PATH`` runs ``ipi_order_from_traj.py`` after the client exits
(LAMMPS may raise on i-PI shutdown; the driver catches that and still post-processes).

**OpenMM parity:** by default (``--ipi-minimize``, on) the driver runs **LAMMPS + UMA
energy minimization** on the same ``read_data`` geometry as OpenMM’s
``LocalEnergyMinimizer`` leg, then overwrites ``ipi/init.xyz`` from the minimized
structure (same ``convert_lammps_to_ipi_xyz.py`` path). Use ``--no-ipi-minimize`` to
restore the previous raw-data ``init.xyz`` behaviour.
By default this uses the **path-integral** estimator (``ice_order_metrics_path_integral``)
and **per-atom** PBC wrap, matching OpenMM / LAMMPS parity. Pass ``--order-centroid``
only if you want Q6/q_tet on a bead-averaged structure (not comparable to OpenMM's RPMD CSV).
Use ``--order-molecular-wrap`` for legacy molecular wrap.
"""
from __future__ import annotations

import argparse
import fcntl
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent

# ── Exclusive process lock (acquired in main() after argparse, so ``--help`` works) ──
_LOCK_PATH = Path("/tmp/run_ipi_uma_rpmd.lock")


def _acquire_run_lock() -> None:
    fh = open(_LOCK_PATH, "w")  # noqa: SIM115 – process lifetime
    try:
        fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print(
            f"[run_ipi_lammps_uma_rpmd] Another instance is already running "
            f"(lock: {_LOCK_PATH}).  Exiting.",
            flush=True,
        )
        sys.exit(0)
    # Keep handle alive for lock scope; assign on module for clarity (optional)
    _acquire_run_lock._fh = fh  # type: ignore[attr-defined]


if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

IPI_DIR = _SCRIPT_DIR / "ipi"
LAMMPS_DIR = _SCRIPT_DIR / "lammps"
IPI_PORT_DEFAULT = 65535
UMA_DRIVER_DEFAULT = IPI_DIR / "uma_ipi_driver.py"
# Must match ``<trajectory ... stride=...>`` and ``<properties ... stride=...>`` in ``_write_ipi_input``.
IPI_TRAJ_OUTPUT_STRIDE = 10
# Base name in XML; i-PI builds ``<outtemplate>.<filename>`` then ``<system.prefix>_<that>``
# (same pattern as ``ice__i-pi.traj`` → default ``ice__i-pi.md``).
IPI_PROPERTIES_FILENAME = "md"

# LAMMPS ``run`` budget as i-PI **client**: i-PI ``<total_steps>`` is authoritative for physics length.
# The client should stop when the server closes the driver connection, not because LAMMPS exhausted
# the same integer as ``total_steps`` first.  Value is below ``2**63`` for signed 64-bit parsers.
LAMMPS_IPI_CLIENT_RUN_CAP_DEFAULT = 9_000_000


# ── Utilities ────────────────────────────────────────────────────────────────

def _wait_for_ipi_ready(log_path: Path, timeout_s: float = 60.0, interval_s: float = 0.2) -> bool:
    """Return True when i-PI has finished initialising and opened its listen socket.

    Watches the i-PI log file for the message that i-PI prints just after opening the
    server socket and starting the ForceField polling thread.  We deliberately avoid
    doing a TCP test-connect here: i-PI has a known bug (bitwise ``|`` instead of
    ``&`` in the handshake condition) that causes it to treat any connection attempt
    — including an immediately-closed probe — as a successful client handshake, which
    leaves a dead ghost client in its list and can stall the ForceField dispatch loop.
    """
    ready_markers = (
        "Starting the polling thread main loop",
        "interfacesocket.open: Created inet socket",
    )
    t0 = time.time()
    while time.time() - t0 < timeout_s:
        try:
            text = log_path.read_text(encoding="utf-8", errors="replace")
            if any(m in text for m in ready_markers):
                return True
        except OSError:
            pass
        time.sleep(interval_s)
    return False


def _tcp_port_in_use(host: str, port: int) -> bool:
    """True if something is already listening on host:port (e.g. stale i-PI).

    Uses /proc/net/tcp instead of a real TCP connect to avoid triggering i-PI's
    ghost-client bug.  Falls back to a short-lived connect only on non-Linux systems.
    """
    try:
        hex_port = f"{port:04X}"
        for path in ("/proc/net/tcp", "/proc/net/tcp6"):
            try:
                for line in Path(path).read_text().splitlines()[1:]:
                    parts = line.split()
                    if len(parts) >= 4 and parts[3] == "0A":
                        local = parts[1]
                        if local.split(":")[-1].upper() == hex_port:
                            return True
            except OSError:
                pass
        return False
    except Exception:
        import socket as _socket
        try:
            with _socket.socket(_socket.AF_INET, _socket.SOCK_STREAM) as s:
                s.settimeout(0.2)
                return s.connect_ex((host, port)) == 0
        except OSError:
            return False


def _backup_ipi_log_if_present(log_path: Path) -> None:
    """Copy existing log to ``stem_backup_YYYYMMDD_HHMMSS.suffix`` before truncating."""
    try:
        if not log_path.is_file() or log_path.stat().st_size == 0:
            return
    except OSError:
        return
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = log_path.parent / f"{log_path.stem}_backup_{ts}{log_path.suffix}"
    shutil.copy2(log_path, backup)
    print(f"Previous i-PI log saved as {backup}", flush=True)


_STEP_IN_COMMENT_RE = re.compile(r"Step:\s*(\d+)")


def _bead0_traj_path(ipi_dir: Path) -> Path | None:
    for name in ("ice__i-pi.traj_00.xyz", "ice__i-pi.traj_0.xyz"):
        p = ipi_dir / name
        if p.is_file():
            return p
    return None


def _ipi_order_from_traj_argv(
    *,
    script_dir: Path,
    traj: Path,
    beads: int,
    dt_fs: float,
    order_csv: Path,
    centroid: bool = False,
    pbc_wrap: str = "atom",
    thermo_prefer_ring_temperature: bool = False,
) -> list[str]:
    """Command argv for ``ipi_order_from_traj.py`` (path-integral order by default).

    *pbc_wrap* is ``atom`` (default, LAMMPS/Fairchem parity), ``molecular``, or ``none``.
    """
    try:
        traj_arg = str(traj.resolve().relative_to(script_dir.resolve()))
    except ValueError:
        traj_arg = str(traj.resolve())
    order_script = script_dir / "ipi_order_from_traj.py"
    cmd: list[str] = [
        sys.executable,
        str(order_script),
        "--traj",
        traj_arg,
        "--beads",
        str(beads),
        "--dt-fs",
        str(dt_fs),
        "-o",
        str(order_csv),
    ]
    if centroid:
        cmd.append("--centroid")
    if pbc_wrap == "none":
        cmd.append("--no-wrap")
    elif pbc_wrap == "molecular":
        cmd.append("--molecular-wrap")
    if thermo_prefer_ring_temperature:
        cmd.append("--thermo-prefer-ring-temperature")
    return cmd


_LAMMPS_IPI_EXIT_MARKER = "Got EXIT message from i-PI"


def _is_lammps_ipi_normal_shutdown(exc: BaseException) -> bool:
    """True when LAMMPS ``fix ipi`` raises after i-PI ends the run and closes the socket."""
    return _LAMMPS_IPI_EXIT_MARKER in str(exc)


def _run_ipi_order_csv_postprocess(
    *,
    script_dir: Path,
    ipi_dir: Path,
    beads: int,
    dt_fs: float,
    order_csv: Path,
    centroid: bool = False,
    pbc_wrap: str = "atom",
    thermo_prefer_ring_temperature: bool = False,
) -> None:
    """Run ``ipi_order_from_traj``; abort if bead-0 trajectory is missing."""
    traj = _bead0_traj_path(ipi_dir)
    if traj is None:
        print(
            "ERROR: --order-csv requested but no bead-0 trajectory found "
            f"({ipi_dir}/ice__i-pi.traj_00.xyz or traj_0.xyz).",
            file=sys.stderr,
            flush=True,
        )
        sys.exit(1)
    order_csv = order_csv.expanduser().resolve()
    order_csv.parent.mkdir(parents=True, exist_ok=True)
    thermo_path = script_dir / "ipi" / "ice__i-pi.md"
    if not thermo_path.is_file():
        print(
            "ERROR: i-PI thermo file missing; cannot merge T_K / PE into order CSV "
            f"({thermo_path}).",
            file=sys.stderr,
            flush=True,
        )
        sys.exit(1)
    cmd = _ipi_order_from_traj_argv(
        script_dir=script_dir,
        traj=traj,
        beads=beads,
        dt_fs=dt_fs,
        order_csv=order_csv,
        centroid=centroid,
        pbc_wrap=pbc_wrap,
        thermo_prefer_ring_temperature=thermo_prefer_ring_temperature,
    )
    cmd.extend(["--thermo", str(thermo_path)])
    est = "centroid structure" if centroid else "path-integral (OpenMM-parity)"
    print(f"Order CSV post-process ({est}): {' '.join(cmd)}", flush=True)
    subprocess.check_call(cmd, cwd=str(script_dir))


def _last_step_bead0_traj(ipi_dir: Path) -> int | None:
    """Scan bead-0 xyz for the last ``Step:`` in a comment line (streaming, O(n) lines)."""
    traj = _bead0_traj_path(ipi_dir)
    if traj is None:
        return None
    last: int | None = None
    try:
        with traj.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                m = _STEP_IN_COMMENT_RE.search(line)
                if m:
                    last = int(m.group(1))
    except OSError:
        return None
    return last


def _verify_traj_matches_requested_steps(
    ipi_dir: Path,
    expected_steps: int,
    stride: int = IPI_TRAJ_OUTPUT_STRIDE,
    *,
    fatal_if_short: bool = False,
) -> None:
    """Print last traj Step vs ``expected_steps``; warn or abort if trajectory is too short."""
    last = _last_step_bead0_traj(ipi_dir)
    if last is None:
        print(
            "Verify: no bead-0 trajectory (ice__i-pi.traj_00.xyz) or no Step: lines.",
            flush=True,
        )
        return
    print(
        f"Verify: last Step in bead-0 traj = {last} "
        f"(requested {expected_steps} i-PI steps, traj output stride {stride}).",
        flush=True,
    )
    if last + stride < expected_steps:
        msg = (
            f"Trajectory ends at Step {last}, well short of {expected_steps} i-PI steps. "
            "Causes: interrupted run (Ctrl+C → SOFTEXIT), OOM kill of i-PI, stale ipi/ files, "
            "or i-PI exiting while LAMMPS continued — not a --prod cap."
        )
        if fatal_if_short:
            print(f"ERROR: {msg}", file=sys.stderr, flush=True)
            sys.exit(1)
        print(f"WARNING: {msg}", file=sys.stderr, flush=True)


def _make_ipi_env() -> dict[str, str]:
    """Build an environment dict with unbuffered output, PyTorch CUDA tweaks,
    and LD_LIBRARY_PATH pointing at the conda/pip NVIDIA libs."""
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    env.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")

    import glob as _glob
    sp = (
        Path(sys.executable).parent.parent
        / "lib"
        / f"python{sys.version_info.major}.{sys.version_info.minor}"
        / "site-packages"
    )
    nvidia_libs = sorted(_glob.glob(str(sp / "nvidia" / "*" / "lib")))
    if nvidia_libs:
        existing = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = ":".join(nvidia_libs) + (":" + existing if existing else "")
    return env


def _preload_mpich_lib_for_pip_lammps() -> None:
    """Pip ``liblammps.so`` links **MPICH** ``libmpi.so.12`` (not Open MPI ``libmpi.so.40``).

    ``conda install openmpi`` does not provide ``.so.12``.  ``liblammps.so`` also sets an
    ``RPATH`` so ``LD_LIBRARY_PATH`` alone is unreliable; preloading via absolute
    ``ctypes.CDLL(..., RTLD_GLOBAL)`` matches what the dynamic linker needs before ``lammps()``.
    """
    import ctypes

    prefix = Path(sys.executable).resolve().parent.parent
    candidates: list[Path] = []
    env_mpi = prefix / "lib" / "libmpi.so.12"
    if env_mpi.is_file():
        candidates.append(env_mpi)
    if prefix.parent.name == "envs":
        base_mpi = prefix.parent.parent / "lib" / "libmpi.so.12"
        if base_mpi.is_file() and base_mpi.resolve() != env_mpi.resolve():
            candidates.append(base_mpi)

    for p in candidates:
        try:
            ctypes.CDLL(str(p), mode=ctypes.RTLD_GLOBAL)
            if p == env_mpi:
                return
            print(
                f"[run_ipi_lammps_uma_rpmd] Pre-loaded MPICH for pip LAMMPS: {p}\n"
                "  (install into this env: conda install mpich — avoids using base Miniconda lib)",
                flush=True,
            )
            return
        except OSError:
            continue


def _pyfftw_alloc_broken() -> bool:
    """True if aligned allocation fails (some conda/pyfftw/numpy builds).

    Tests in the CURRENT env first. If pyfftw is absent here, also probes the
    base-env Python that ``i-pi`` actually uses (``i-pi`` lives in the base
    Miniconda env, not necessarily in the active conda env).
    """
    try:
        import pyfftw

        if hasattr(pyfftw, "empty_aligned"):
            pyfftw.empty_aligned((4, 9), dtype="float64", n=16)
        else:
            pyfftw.n_byte_align_empty((4, 9), 16, "float64")
        return False
    except ImportError:
        pass
    except (ValueError, OSError, RuntimeError):
        return True

    ipi_bin = shutil.which("i-pi")
    if ipi_bin is None:
        return False
    ipi_python = Path(ipi_bin).resolve().parent / "python"
    if not ipi_python.is_file():
        ipi_python = Path(ipi_bin).resolve().parent / "python3"
    if not ipi_python.is_file():
        return False
    try:
        result = subprocess.run(
            [str(ipi_python), "-c",
             "import pyfftw; pyfftw.empty_aligned((4,9), dtype='float64', n=16)"],
            capture_output=True, timeout=10,
        )
        return result.returncode != 0
    except Exception:
        return True


def _ipi_subprocess_env(base: dict[str, str]) -> dict[str, str]:
    """Environment for the ``i-pi`` subprocess; prepends a stub so broken PyFFTW is skipped.

    i-PI only catches ``ImportError`` around ``import pyfftw``; a broken install can
    raise ``ValueError`` during allocation instead.  A stub ``pyfftw`` package that
    raises ``ImportError`` forces the NumPy FFT path.
    """
    env = dict(base)
    stub = _SCRIPT_DIR / "_ipi_pyfftw_stub"
    if stub.is_dir() and _pyfftw_alloc_broken():
        prepend = str(stub.resolve())
        old_pp = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = prepend + (":" + old_pp if old_pp else "")
        print(
            "[run_ipi_lammps_uma_rpmd] Broken PyFFTW detected; i-PI subprocess uses NumPy FFT "
            f"(stub on PYTHONPATH: {prepend})",
            flush=True,
        )
    return env


def _write_ipi_input(
    output_path: Path,
    molecules: int,
    beads: int,
    total_steps: int,
    *,
    dt_fs: float = 0.5,
    seed: int = 284759,
    port: int = IPI_PORT_DEFAULT,
    thermostat_mode: str = "pile_g",
    tau_fs: float = 1000.0,
) -> None:
    """Write i-PI input.xml with natoms=3*molecules, nbeads=beads.

    thermostat_mode: i-PI ``pile_g`` (SVR on centroid + PILE on internals, matches
    OpenMM ``RPMDIntegrator.PileG`` intent) or ``pile_l`` (Langevin on centroid).
    """
    xml = f"""\
<?xml version="1.0" encoding="UTF-8"?>
<simulation verbosity="medium">
  <!-- RPMD: {molecules} mol, {beads} beads, 243 K, thermostat={thermostat_mode}, {dt_fs} fs -->
  <total_steps>{total_steps}</total_steps>
  <prng>
    <seed>{seed}</seed>
  </prng>
  <output>
    <trajectory filename="traj" stride="{IPI_TRAJ_OUTPUT_STRIDE}" format="xyz">positions</trajectory>
    <trajectory filename="vtraj" stride="{IPI_TRAJ_OUTPUT_STRIDE}" format="xyz">velocities</trajectory>
    <properties stride="{IPI_TRAJ_OUTPUT_STRIDE}" filename="{IPI_PROPERTIES_FILENAME}">
      [ step, time{{picosecond}}, temperature{{kelvin}}, temperature(nm=0){{kelvin}}, potential{{electronvolt}} ]
    </properties>
    <checkpoint filename="RESTART" stride="{total_steps}"/>
  </output>

  <ffsocket mode="inet" name="lammps">
    <address>127.0.0.1</address>
    <port>{port}</port>
    <timeout>3600</timeout>
  </ffsocket>

  <system prefix="ice_">
    <forces>
      <force forcefield="lammps"/>
    </forces>
    <ensemble>
      <temperature units="kelvin">243</temperature>
    </ensemble>
    <motion mode="dynamics">
      <dynamics mode="nvt">
        <timestep units="femtosecond">{dt_fs}</timestep>
        <thermostat mode="{thermostat_mode}">
          <tau units="femtosecond">{tau_fs}</tau>
        </thermostat>
      </dynamics>
    </motion>
    <initialize nbeads="{beads}">
      <file mode="xyz" units="angstrom">init.xyz</file>
      <velocities mode="thermal" units="kelvin">243</velocities>
    </initialize>
  </system>
</simulation>
"""
    output_path.write_text(xml)


# ── Client launchers ─────────────────────────────────────────────────────────

def _patch_fairchem_lammps_atomic_numbers_tensor() -> None:
    """FairChem ``atomic_data_from_lammps_data`` uses ``torch.tensor(atomic_numbers)`` without dtype.

    On NumPy 2 + some PyTorch builds this raises ``RuntimeError: Could not infer dtype of numpy.int64``.
    The error occurs inside a LAMMPS ctypes callback and is swallowed (``Exception ignored``), so
    UMA forces are never applied and i-PI may stop far short of ``total_steps``.
    """
    import fairchem.lammps.lammps_fc as lfc

    if getattr(lfc, "_uma_ipi_atomic_numbers_dtype_patch", False):
        return

    import numpy as np
    import torch
    from fairchem.core.datasets.atomic_data import AtomicData

    def atomic_data_from_lammps_data(
        x,
        atomic_numbers,
        nlocal,
        cell,
        periodicity,
        task_name,
        charge: int = 0,
        spin: int = 0,
    ):
        pos = torch.as_tensor(x, dtype=torch.float32)
        pbc = torch.tensor(periodicity, dtype=torch.bool).unsqueeze(0)
        edge_index = torch.empty((2, 0), dtype=torch.long)
        cell_offsets = torch.empty((0, 3), dtype=torch.float32)
        nedges = torch.tensor([0], dtype=torch.long)
        tags = torch.zeros(nlocal, dtype=torch.long)
        fixed = torch.zeros(nlocal, dtype=torch.long)
        batch = torch.zeros(nlocal, dtype=torch.long)
        z = torch.as_tensor(np.asarray(atomic_numbers, dtype=np.int64), dtype=torch.long)
        return AtomicData(
            pos=pos,
            atomic_numbers=z,
            cell=cell,
            pbc=pbc,
            natoms=torch.tensor([nlocal], dtype=torch.long),
            edge_index=edge_index,
            cell_offsets=cell_offsets,
            nedges=nedges,
            charge=torch.LongTensor([charge]),
            spin=torch.LongTensor([spin]),
            fixed=fixed,
            tags=tags,
            batch=batch,
            dataset=[task_name],
        )

    lfc.atomic_data_from_lammps_data = atomic_data_from_lammps_data
    lfc._uma_ipi_atomic_numbers_dtype_patch = True


def _patch_fairchem_escn_moe_numpy_bool_indexing() -> None:
    """FairChem UMA ``escn_moe`` uses NumPy boolean masks to index PyTorch tensors.

    With NumPy 2 + recent PyTorch, ``tensor[numpy_bool_mask]`` raises
    ``RuntimeError: Could not infer dtype of numpy.bool_`` inside LAMMPS ``fix external``
    callbacks (errors are swallowed, so forces are never applied).
    """
    try:
        from fairchem.core.common.utils import conditional_grad
        from fairchem.core.models.uma import escn_moe
    except ImportError:
        return
    if getattr(escn_moe, "_uma_escn_moe_bool_index_patch", False):
        return

    import numpy as np
    import torch

    def _bool_index(mask_np: np.ndarray, device: torch.device) -> torch.Tensor:
        return torch.as_tensor(
            np.asarray(mask_np, dtype=np.bool_),
            device=device,
            dtype=torch.bool,
        )

    @conditional_grad(torch.enable_grad())
    def _moe_forward(self, data, emb):  # noqa: ANN001
        self.global_mole_tensors.mole_sizes = torch.zeros(
            data.natoms.shape[0], dtype=torch.int, device=emb["batch"].device
        ).scatter(0, emb["batch"], 1, reduce="add")
        self.global_mole_tensors.natoms = emb["batch"].shape[0]
        data_batch_full = data.batch_full.cpu()
        self.global_mole_tensors.expert_mixing_coefficients = (
            torch.zeros(
                data.natoms.shape[0],
                len(self.dataset_name_to_exp),
                dtype=data.pos.dtype,
            )
            .scatter(
                1,
                torch.tensor(
                    [
                        self.dataset_name_to_exp[dataset_name]
                        for dataset_name in data.dataset
                    ],
                ).unsqueeze(1),
                1.0,
            )
            .to(data.pos.device)
        )
        head_output = self.head(data, emb)
        np_dataset_names = np.array(data.dataset)
        full_output = {}
        for dataset_name in self.dataset_names:
            dataset_mask = np_dataset_names == dataset_name
            for key, mole_output_tensor in head_output.items():
                output_tensor = mole_output_tensor.new_zeros(mole_output_tensor.shape)
                if dataset_mask.any():
                    mask_t = _bool_index(dataset_mask, output_tensor.device)
                    if output_tensor.shape[0] == dataset_mask.shape[0]:
                        output_tensor[mask_t] = mole_output_tensor[mask_t]
                    else:
                        idx_cpu = torch.where(
                            torch.as_tensor(np.asarray(dataset_mask, dtype=np.bool_))
                        )[0]
                        atoms_mask = torch.isin(data_batch_full, idx_cpu)
                        output_tensor[atoms_mask] = mole_output_tensor[atoms_mask]
                full_output[f"{dataset_name}_{key}"] = (
                    {key: output_tensor} if self.wrap_property else output_tensor
                )
        return full_output

    @conditional_grad(torch.enable_grad())
    def _single_forward(self, data, emb):  # noqa: ANN001
        data_batch_full = data.batch_full.cpu()
        head_output = self.head(data, emb)
        if self.merged_on_dataset is not None:
            full_output = {}
            for key in head_output:
                full_output[f"{self.merged_on_dataset}_{key}"] = (
                    {key: head_output[key]} if self.wrap_property else head_output[key]
                )
                nan_tensor = head_output[key].new_full(
                    head_output[key].shape, float("nan")
                )
                for dataset in self.non_merged_dataset_names:
                    full_output[f"{dataset}_{key}"] = (
                        {key: nan_tensor} if self.wrap_property else nan_tensor
                    )
            return full_output
        assert set(data.dataset) <= set(self.dataset_names), (
            f"Input dataset names: {set(data.dataset)} must be a strict subset of "
            f"model's valid datset names: {set(self.dataset_names)}"
        )
        np_dataset_names = np.array(data.dataset)
        full_output = {}
        for dataset_name in self.dataset_names:
            dataset_mask = np_dataset_names == dataset_name
            for key, head_output_tensor in head_output.items():
                output_tensor = head_output_tensor.new_zeros(head_output_tensor.shape)
                if dataset_mask.any():
                    mask_t = _bool_index(dataset_mask, output_tensor.device)
                    if output_tensor.shape[0] == dataset_mask.shape[0]:
                        output_tensor[mask_t] = head_output_tensor[mask_t]
                    else:
                        idx_cpu = torch.where(
                            torch.as_tensor(np.asarray(dataset_mask, dtype=np.bool_))
                        )[0]
                        atoms_mask = torch.isin(data_batch_full, idx_cpu)
                        output_tensor[atoms_mask] = head_output_tensor[atoms_mask]
                full_output[f"{dataset_name}_{key}"] = (
                    {key: output_tensor} if self.wrap_property else output_tensor
                )
        return full_output

    escn_moe.DatasetSpecificMoEWrapper.forward = _moe_forward
    escn_moe.DatasetSpecificSingleHeadWrapper.forward = _single_forward
    escn_moe._uma_escn_moe_bool_index_patch = True


# Defaults aligned with ``lammps/in.ice_uma.lmp`` and OpenMM ``maxIterations=2000`` in
# ``test_uma_ice_rpmd.py`` (LAMMPS ``maxiter`` / ``maxeval`` are separate caps).
_IPI_MIN_ETOL_DEFAULT = 0.0
_IPI_MIN_FTOL_DEFAULT = 1.0e-6
_IPI_MIN_MAXITER_DEFAULT = 2000
_IPI_MIN_MAXEVAL_DEFAULT = 50000


def _lammps_uma_minimize_to_ipi_init(
    *,
    data_path: Path,
    init_xyz: Path,
    model: str,
    device: str,
    task_name: str,
    ipi_log: Path,
    etol: float,
    ftol: float,
    maxiter: int,
    maxeval: int,
) -> None:
    """LAMMPS + fix external UMA minimize ``data_path``, then write ``init_xyz`` for i-PI.

    Mirrors the FairChem pattern used in ``lammps/run_lammps_uma_ice.py``: external
    fix must be active before ``minimize``.  Does **not** start ``fix ipi``.
    """
    try:
        from fairchem.lammps.lammps_fc import (
            FIX_EXT_ID,
            FIX_EXTERNAL_CMD,
            FixExternalCallback,
        )
        from lammps import lammps
    except ImportError as exc:
        print(
            f"UMA minimization requires fairchem-lammps + lammps: {exc}\n"
            "  Install or use --no-ipi-minimize.",
            file=sys.stderr,
            flush=True,
        )
        raise
    _patch_fairchem_lammps_atomic_numbers_tensor()
    # escn_moe numpy.bool_ indexing is fixed in the installed escn_moe.py directly;
    # the runtime monkeypatch is no longer needed and was overriding the fix.
    try:
        from fairchem.core import pretrained_mlip
    except ImportError:
        from fairchem.core.calculate import pretrained_mlip

    print(
        f"UMA minimization (LAMMPS+fix external) before i-PI: "
        f"etol={etol} ftol={ftol} maxiter={maxiter} maxeval={maxeval}",
        flush=True,
    )
    print(f"Loading UMA predictor ({model}) on {device} ...", flush=True)
    predictor = pretrained_mlip.get_predict_unit(model, device=device)

    log_file = str(ipi_log.parent.parent / "pipeline_out" / "lammps_uma_minimize.log")
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    machine = os.environ.get("LAMMPS_MACHINE_NAME")
    _preload_mpich_lib_for_pip_lammps()
    try:
        lmp = lammps(name=machine, cmdargs=["-nocite", "-log", log_file, "-echo", "screen"])
    except OSError as exc:
        if "libmpi.so.12" in str(exc):
            print(
                "LAMMPS (pip) needs MPICH libmpi.so.12 (Open MPI only ships libmpi.so.40).\n"
                "  Fix: conda install -n <env> mpich\n"
                f"  Underlying error: {exc}",
                file=sys.stderr,
                flush=True,
            )
        raise
    lmp._predictor = predictor
    lmp._task_name = task_name
    try:
        lmp.commands_list([
            "units metal",
            "atom_style atomic",
            "atom_modify sort 0 0.0",
            "boundary p p p",
            f"read_data {data_path.resolve()}",
            "comm_modify cutoff 12.0",
            "pair_style zero 1.0",
            "pair_coeff * *",
        ])
        lmp.command(FIX_EXTERNAL_CMD)
        lmp.set_fix_external_callback(FIX_EXT_ID, FixExternalCallback(charge=0, spin=1), lmp)
        print("UMA fix external registered → minimize", flush=True)
        lmp.command(f"minimize {etol} {ftol} {maxiter} {maxeval}")

        # Match ``fix ipi`` / i-PI convention: symmetric orthorhombic box [-L/2, L/2].
        # ``convert_lammps_to_ipi_xyz`` centers Cartesian coords the same way; keeping
        # the written data file aligned avoids a mismatched internal state before the
        # first socket exchange.
        boxlo, boxhi, xy, yz, xz, *_rest = lmp.extract_box()
        lx = boxhi[0] - boxlo[0]
        ly = boxhi[1] - boxlo[1]
        lz = boxhi[2] - boxlo[2]
        cx = 0.5 * (boxlo[0] + boxhi[0])
        cy = 0.5 * (boxlo[1] + boxhi[1])
        cz = 0.5 * (boxlo[2] + boxhi[2])
        if abs(xy) + abs(xz) + abs(yz) < 1e-12:
            lmp.command(f"displace_atoms all move {-cx:.16f} {-cy:.16f} {-cz:.16f}")
            lmp.command(
                f"change_box all x final {-0.5 * lx:.16f} {0.5 * lx:.16f} "
                f"y final {-0.5 * ly:.16f} {0.5 * ly:.16f} "
                f"z final {-0.5 * lz:.16f} {0.5 * lz:.16f}"
            )

        # Persist for LAMMPS ``read_data`` in the i-PI client so the initial topology
        # matches ``init.xyz`` (same minimized geometry as the socket-driven steps).
        min_data = IPI_DIR / "lammps_data_ipi_client.data"
        lmp.command(f"write_data {min_data.resolve()}")
        print(f"Wrote minimized LAMMPS data {min_data}", flush=True)
    finally:
        del lmp._predictor
        lmp.close()

    converter = IPI_DIR / "convert_lammps_to_ipi_xyz.py"
    subprocess.check_call(
        [
            sys.executable,
            str(converter),
            "--data",
            str(min_data.resolve()),
            "-o",
            str(init_xyz),
        ],
        cwd=str(_SCRIPT_DIR),
    )
    print(f"Wrote minimized i-PI template {init_xyz}", flush=True)


def _run_client_lammps(
    *,
    data_path: Path,
    port: int,
    lammps_run_steps: int,
    model: str,
    device: str,
    task_name: str,
    ipi_env: dict[str, str],
    ipi_log: Path,
) -> None:
    """Start LAMMPS as an i-PI force client with UMA via fix external.

    ``lammps_run_steps`` is the LAMMPS ``run`` budget (typically huge); i-PI ``<total_steps>``
    decides when the server stops and closes the socket.

    Uses the LAMMPS Python API so that fix external pf/callback can be
    registered in-process (no separate LAMMPS binary needed).
    """
    try:
        from fairchem.lammps.lammps_fc import (
            FIX_EXT_ID,
            FIX_EXTERNAL_CMD,
            FixExternalCallback,
        )
        from lammps import lammps
    except ImportError as exc:
        print(
            f"LAMMPS client requires: pip install fairchem-lammps fairchem-core lammps\n  {exc}",
            file=sys.stderr,
        )
        sys.exit(1)
    _patch_fairchem_lammps_atomic_numbers_tensor()
    # escn_moe numpy.bool_ indexing is fixed in the installed escn_moe.py directly.
    try:
        from fairchem.core import pretrained_mlip
    except ImportError:
        from fairchem.core.calculate import pretrained_mlip

    print(f"Loading UMA predictor ({model}) on {device} ...", flush=True)
    predictor = pretrained_mlip.get_predict_unit(model, device=device)

    log_file = str(ipi_log.parent.parent / "pipeline_out" / "lammps_ipi_client.log")
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    machine = os.environ.get("LAMMPS_MACHINE_NAME")
    _preload_mpich_lib_for_pip_lammps()
    try:
        lmp = lammps(name=machine, cmdargs=["-nocite", "-log", log_file, "-echo", "screen"])
    except OSError as exc:
        if "libmpi.so.12" in str(exc):
            print(
                "LAMMPS (pip) needs MPICH libmpi.so.12 (Open MPI only ships libmpi.so.40).\n"
                "  Fix: conda install -n <env> mpich\n"
                f"  Underlying error: {exc}",
                file=sys.stderr,
                flush=True,
            )
        raise
    lmp._predictor = predictor
    lmp._task_name = task_name

    # ``write_data`` from the i-PI minimize leg can emit a ``Pair Coeffs`` section;
    # LAMMPS requires ``pair_style`` before ``read_data`` in that case.  Do not call
    # ``pair_coeff`` before ``read_data`` — the simulation box is not defined yet (err0033).
    lmp.commands_list([
        "units metal",
        "atom_style atomic",
        "atom_modify sort 0 0.0",
        "boundary p p p",
        "pair_style zero 1.0",
        f"read_data {data_path.resolve()}",
        "comm_modify cutoff 12.0",
        "pair_coeff * *",
    ])

    lmp.command(FIX_EXTERNAL_CMD)
    base_cb = FixExternalCallback(charge=0, spin=1)

    class _DiagCallback:
        """Thin wrapper that logs energy/positions for the first few callbacks."""

        def __init__(self, inner, max_log: int = 3):
            self._inner = inner
            self._max_log = max_log
            self._n = 0

        def __call__(self, lmp_obj, ntimestep, nlocal, tag, x, f):
            self._inner(lmp_obj, ntimestep, nlocal, tag, x, f)
            if self._n < self._max_log:
                boxlo, boxhi, *_ = lmp_obj.extract_box()
                import numpy as _np
                pos_arr = _np.array(x[:nlocal], dtype=_np.float64)
                e_global = lmp_obj.get_thermo("pe")
                print(
                    f"[DIAG cb#{self._n}] step={ntimestep}  nlocal={nlocal}  "
                    f"fix_ext_energy={lmp_obj.get_thermo('pe'):.6f} eV  "
                    f"box=[{boxlo[0]:.4f}..{boxhi[0]:.4f}, "
                    f"{boxlo[1]:.4f}..{boxhi[1]:.4f}, "
                    f"{boxlo[2]:.4f}..{boxhi[2]:.4f}]  "
                    f"pos_range_x=[{pos_arr[:,0].min():.4f}, {pos_arr[:,0].max():.4f}]  "
                    f"pos_range_y=[{pos_arr[:,1].min():.4f}, {pos_arr[:,1].max():.4f}]  "
                    f"pos_range_z=[{pos_arr[:,2].min():.4f}, {pos_arr[:,2].max():.4f}]  "
                    f"fmax={_np.abs(_np.array(f[:nlocal])).max():.6f} eV/A",
                    flush=True,
                )
            self._n += 1

    diag_cb = _DiagCallback(base_cb, max_log=5)
    lmp.set_fix_external_callback(FIX_EXT_ID, diag_cb, lmp)
    print("UMA fix external registered (with diagnostics)", flush=True)

    lmp.command(f"fix ipi_client all ipi 127.0.0.1 {port}")
    lmp.command(f"run {lammps_run_steps}")

    del lmp._predictor
    lmp.close()
    print("LAMMPS client finished", flush=True)


def _run_client_python(
    *,
    init_xyz: Path,
    uma_driver_path: Path,
    port: int,
    model: str,
    device: str,
    task_name: str,
    ipi_env: dict[str, str],
    ipi_log: Path,
) -> None:
    """Start the i-pi-py_driver subprocess with the custom UMA PES driver."""
    driver_opts = (
        f"template={init_xyz.resolve()},"
        f"model={model},"
        f"device={device},"
        f"task_name={task_name}"
    )
    driver_cmd = [
        "i-pi-py_driver",
        "-P", str(uma_driver_path),
        "-m", "custom",
        "-a", "127.0.0.1",
        "-p", str(port),
        "-o", driver_opts,
    ]
    print(f"Launching UMA driver: {' '.join(driver_cmd)}", flush=True)
    driver_log_path = ipi_log.parent.parent / "pipeline_out" / "uma_driver.log"
    driver_log_path.parent.mkdir(parents=True, exist_ok=True)
    driver_log_f = open(driver_log_path, "w", encoding="utf-8", buffering=1)
    try:
        proc = subprocess.Popen(
            driver_cmd,
            stdout=driver_log_f,
            stderr=subprocess.STDOUT,
            env=ipi_env,
        )
        print(f"UMA driver started (pid {proc.pid}); log: {driver_log_path}", flush=True)
        proc.wait()
        rc = proc.returncode
        if rc != 0:
            print(f"UMA driver exited with code {rc}", file=sys.stderr, flush=True)
        else:
            print("UMA driver finished successfully", flush=True)
    finally:
        if proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=5)
        driver_log_f.close()


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Run i-PI RPMD server + force client (LAMMPS or Python UMA driver)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Client selection
    ap.add_argument(
        "--client",
        choices=["lammps", "python"],
        default="lammps",
        help="Force client: lammps (default, LAMMPS Python API + fix ipi + fix external UMA) "
        "or python (i-pi-py_driver + ASE FAIRChemCalculator).",
    )

    # Simulation parameters
    ap.add_argument("--steps", type=int, default=None,
                    help="Override: MD steps (default: from --prod and --dt-fs)")
    ap.add_argument("--prod", type=float, default=None,
                    help="Production time in ps (e.g. 100); used with --dt-fs to set steps")
    ap.add_argument("--dt-fs", type=float, default=0.5,
                    help="Timestep in fs (default 0.5)")
    ap.add_argument("--seed", type=int, default=284759, help="i-PI PRNG seed")
    ap.add_argument(
        "--port", type=int, default=None,
        help=f"i-PI socket port (default {IPI_PORT_DEFAULT}, use for SLURM array isolation)",
    )
    ap.add_argument("--molecules", type=int, default=64,
                    help="Water molecules (default 64 = ice Ih 2x2x2)")
    ap.add_argument("--beads", type=int, default=4,
                    help="RPMD beads (default 4)")
    ap.add_argument("--device", default="cuda",
                    help="Device for UMA (cuda/cpu)")
    ap.add_argument("--model", default="uma-s-1p1")
    ap.add_argument(
        "--data", type=Path, default=None,
        help="LAMMPS data file (default: lammps/data.ice_uma_{molecules})",
    )
    ap.add_argument("--task-name", default="omol",
                    help="UMA task name (omol, omat, oc20, ...)")

    # Python client only
    ap.add_argument(
        "--uma-driver", type=Path, default=UMA_DRIVER_DEFAULT,
        help=f"Path to UMA i-PI PES driver (default: {UMA_DRIVER_DEFAULT})",
    )

    # i-PI configuration
    ap.add_argument("--ipi-input", type=Path, default=IPI_DIR / "input.xml")
    ap.add_argument(
        "--ipi-log", type=Path, default=IPI_DIR / "i-pi_run.log",
        help="Log file for i-PI stdout/stderr (default: ipi/i-pi_run.log).",
    )
    ap.add_argument(
        "--ipi-thermostat", type=str, default="pile_g",
        choices=["pile_g", "pile_l"],
        help="i-PI ring-polymer thermostat: pile_g (default) or pile_l.",
    )
    ap.add_argument(
        "--ipi-tau-fs", type=float, default=1000.0,
        help="i-PI thermostat tau in fs (default: 1000)",
    )
    ap.add_argument(
        "--lammps-run-steps",
        type=int,
        default=None,
        metavar="N",
        help=(
            "LAMMPS client ``run N`` step cap (default: a huge value; i-PI ``<total_steps>`` "
            "from --steps/--prod is authoritative).  Lower this only if LAMMPS rejects the default."
        ),
    )
    ap.add_argument(
        "--order-csv",
        type=Path,
        default=None,
        help=(
            "After the run, write ice order CSV via ipi_order_from_traj.py. "
            "Default estimator matches OpenMM UMA RPMD (path-integral / mean over beads of Q6 and q_tet)."
        ),
    )
    ap.add_argument(
        "--order-centroid",
        action="store_true",
        help=(
            "With --order-csv: use centroid-bead estimator instead of path-integral "
            "(not comparable to OpenMM ice_order_*_rpmd.csv)."
        ),
    )
    ap.add_argument(
        "--order-no-wrap",
        action="store_true",
        help="With --order-csv: pass --no-wrap to ipi_order_from_traj (no PBC remap).",
    )
    ap.add_argument(
        "--order-molecular-wrap",
        action="store_true",
        help=(
            "With --order-csv: pass --molecular-wrap to ipi_order_from_traj (legacy H2O group shift). "
            "Default post-process uses per-atom wrap (LAMMPS/Fairchem parity)."
        ),
    )
    ap.add_argument(
        "--order-thermo-prefer-ring-temperature",
        action="store_true",
        help=(
            "With --order-csv: pass --thermo-prefer-ring-temperature to ipi_order_from_traj "
            "(only when vtraj velocity files are missing: prefer ring ``temperature`` over "
            "``temperature(nm=0)`` for merged T_K)."
        ),
    )
    ap.add_argument(
        "--no-ipi-minimize",
        action="store_true",
        help=(
            "Skip LAMMPS+UMA energy minimization before i-PI; write init.xyz from raw LAMMPS data "
            "(legacy behaviour). Default is to minimize so the i-PI leg matches OpenMM's minimized "
            "starting geometry."
        ),
    )
    ap.add_argument(
        "--ipi-min-etol",
        type=float,
        default=_IPI_MIN_ETOL_DEFAULT,
        metavar="EV",
        help="LAMMPS minimize energy tolerance (metal units, eV). Default 0 (ignore).",
    )
    ap.add_argument(
        "--ipi-min-ftol",
        type=float,
        default=_IPI_MIN_FTOL_DEFAULT,
        metavar="EV_PER_ANG",
        help="LAMMPS minimize force tolerance (metal, eV/Å). Default 1e-6 (matches in.ice_uma.lmp).",
    )
    ap.add_argument(
        "--ipi-min-maxiter",
        type=int,
        default=_IPI_MIN_MAXITER_DEFAULT,
        metavar="N",
        help="LAMMPS minimize max iterations. Default 2000 (OpenMM LocalEnergyMinimizer cap).",
    )
    ap.add_argument(
        "--ipi-min-maxeval",
        type=int,
        default=_IPI_MIN_MAXEVAL_DEFAULT,
        metavar="N",
        help="LAMMPS minimize max energy/force evaluations. Default 50000 (in.ice_uma.lmp).",
    )
    args = ap.parse_args()
    if args.order_no_wrap and args.order_molecular_wrap:
        ap.error("use only one of --order-no-wrap and --order-molecular-wrap")

    _acquire_run_lock()

    # ── Resolve data path ────────────────────────────────────────────────
    data_path = args.data or (LAMMPS_DIR / f"data.ice_uma_{args.molecules}")
    if not data_path.is_file():
        print(
            f"Missing data file: {data_path}. Example (64 molecules, 2x2x2): "
            f"python lammps/build_lammps_ice_data.py --nx 2 --ny 2 --nz 2 -o {data_path}",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.client == "python":
        uma_driver_path = args.uma_driver.resolve()
        if not uma_driver_path.is_file():
            print(f"Missing UMA driver: {uma_driver_path}", file=sys.stderr)
            sys.exit(1)

    # ── LAMMPS data -> i-PI init.xyz (optional UMA minimize, OpenMM parity) ─
    init_xyz = IPI_DIR / "init.xyz"
    converter = IPI_DIR / "convert_lammps_to_ipi_xyz.py"
    if not converter.is_file():
        print(f"Missing converter: {converter}", file=sys.stderr)
        sys.exit(1)

    client_lammps_data = data_path.resolve()
    minimized_client_data = IPI_DIR / "lammps_data_ipi_client.data"

    if args.no_ipi_minimize:
        minimized_client_data.unlink(missing_ok=True)
        subprocess.check_call(
            [
                sys.executable,
                str(converter),
                "--data",
                str(data_path),
                "-o",
                str(init_xyz),
            ],
            cwd=str(_SCRIPT_DIR),
        )
        print(f"Wrote {init_xyz} ({3 * args.molecules} atoms) [--no-ipi-minimize, raw data]", flush=True)
    else:
        try:
            _lammps_uma_minimize_to_ipi_init(
                data_path=data_path,
                init_xyz=init_xyz,
                model=args.model,
                device=args.device,
                task_name=args.task_name,
                ipi_log=args.ipi_log,
                etol=float(args.ipi_min_etol),
                ftol=float(args.ipi_min_ftol),
                maxiter=int(args.ipi_min_maxiter),
                maxeval=int(args.ipi_min_maxeval),
            )
            client_lammps_data = minimized_client_data.resolve()
        except ImportError:
            if args.client == "lammps":
                print(
                    "ERROR: --client lammps with default i-PI minimize requires fairchem-lammps + lammps.",
                    file=sys.stderr,
                    flush=True,
                )
                sys.exit(1)
            print(
                "WARNING: UMA minimization skipped (import error); using raw LAMMPS data -> init.xyz.",
                file=sys.stderr,
                flush=True,
            )
            minimized_client_data.unlink(missing_ok=True)
            subprocess.check_call(
                [
                    sys.executable,
                    str(converter),
                    "--data",
                    str(data_path),
                    "-o",
                    str(init_xyz),
                ],
                cwd=str(_SCRIPT_DIR),
            )
            print(f"Wrote {init_xyz} ({3 * args.molecules} atoms) [raw fallback]", flush=True)

    # ── Clean stale i-PI outputs ─────────────────────────────────────────
    for stale in IPI_DIR.glob("ice_*i-pi*"):
        stale.unlink()
        print(f"  Removed stale: {stale.name}")
    for stale in IPI_DIR.glob("RESTART*"):
        stale.unlink()
        print(f"  Removed stale: {stale.name}")

    port = args.port if args.port is not None else IPI_PORT_DEFAULT

    # ── Resolve steps ────────────────────────────────────────────────────
    if args.steps is not None:
        steps = args.steps
    elif args.prod is not None:
        steps = int(args.prod * 1000.0 / args.dt_fs)
    else:
        steps = 1000

    sim_time_ps = steps * float(args.dt_fs) / 1000.0
    lammps_run_steps = (
        args.lammps_run_steps
        if args.lammps_run_steps is not None
        else LAMMPS_IPI_CLIENT_RUN_CAP_DEFAULT
    )
    print(
        f"Resolved: {steps} i-PI <total_steps> ≈ {sim_time_ps:.6f} ps (dt_fs={args.dt_fs} fs); "
        f"LAMMPS client run cap: {lammps_run_steps} (i-PI ends the job when done). "
        "If trajectories end at a much smaller Step:, the run was interrupted — not a --prod cap.",
        flush=True,
    )

    # ── Write i-PI input.xml ─────────────────────────────────────────────
    _write_ipi_input(
        args.ipi_input,
        molecules=args.molecules,
        beads=args.beads,
        total_steps=steps,
        dt_fs=args.dt_fs,
        seed=args.seed,
        port=port,
        thermostat_mode=args.ipi_thermostat,
        tau_fs=args.ipi_tau_fs,
    )
    print(f"Wrote {args.ipi_input} ({args.molecules} mol, {args.beads} beads)")

    # ── Port check ───────────────────────────────────────────────────────
    bind_host = "127.0.0.1"
    if _tcp_port_in_use(bind_host, port):
        print(
            f"ERROR: {bind_host}:{port} is already in use.\n"
            f"  Diagnose: ss -tlnp | grep {port}\n"
            "  Stop:     pkill -f 'i-pi.*input.xml'",
            file=sys.stderr,
        )
        sys.exit(1)

    # ── Launch i-PI server ───────────────────────────────────────────────
    ipi_env = _make_ipi_env()
    ipi_server_env = _ipi_subprocess_env(ipi_env)

    args.ipi_log.parent.mkdir(parents=True, exist_ok=True)
    _backup_ipi_log_if_present(args.ipi_log)
    ipi_log_f = open(args.ipi_log, "w", encoding="utf-8", buffering=1)
    ipi_proc: subprocess.Popen | None = None
    ipi_holder: list[subprocess.Popen | None] = [None]
    _ipi_signal_handlers_installed = False
    prev_sigint: signal.Handlers | None = None
    prev_sigterm: signal.Handlers | None = None

    def _orchestrator_signal_handler(signum: int, frame) -> None:  # noqa: ARG001
        """Terminate i-PI on parent signals so the client does not outlive the server alone."""
        p = ipi_holder[0]
        if p is not None and p.poll() is None:
            try:
                p.terminate()
            except OSError:
                pass
        if signum == signal.SIGINT:
            raise KeyboardInterrupt
        raise SystemExit(143 if signum == signal.SIGTERM else 128 + signum)

    try:
        prev_sigint = signal.signal(signal.SIGINT, _orchestrator_signal_handler)
        prev_sigterm = signal.signal(signal.SIGTERM, _orchestrator_signal_handler)
        _ipi_signal_handlers_installed = True
    except ValueError:
        # Only the main thread may register signals.
        prev_sigint = None
        prev_sigterm = None

    try:
        ipi_proc = subprocess.Popen(
            ["i-pi", str(args.ipi_input)],
            cwd=str(IPI_DIR),
            stdout=ipi_log_f,
            stderr=subprocess.STDOUT,
            text=True,
            env=ipi_server_env,
            start_new_session=True,
        )
        ipi_holder[0] = ipi_proc
        print(
            f"Started i-PI (pid {ipi_proc.pid}) from {IPI_DIR}; "
            f"log: {args.ipi_log.resolve()}",
            flush=True,
        )

        if not _wait_for_ipi_ready(args.ipi_log, timeout_s=60):
            print(
                "i-PI did not become ready in time (check log for errors).\n"
                "  Common causes: broken PyFFTW (orchestrator auto-disables when detected; "
                "else: pip uninstall pyfftw), port in use, bad input.xml\n"
                f"  Log: {args.ipi_log.resolve()}",
                file=sys.stderr,
            )
            if ipi_proc is not None:
                ipi_proc.terminate()
            sys.exit(1)
        print(f"i-PI ready: {bind_host}:{port}", flush=True)
        time.sleep(2.0)

        # ── Launch force client ──────────────────────────────────────────
        print(f"Starting {args.client} client ...", flush=True)

        if args.client == "lammps":
            try:
                _run_client_lammps(
                    data_path=client_lammps_data,
                    port=port,
                    lammps_run_steps=lammps_run_steps,
                    model=args.model,
                    device=args.device,
                    task_name=args.task_name,
                    ipi_env=ipi_env,
                    ipi_log=args.ipi_log,
                )
            except Exception as exc:
                if _is_lammps_ipi_normal_shutdown(exc):
                    print(
                        "LAMMPS: i-PI closed the socket after the run "
                        "(fix ipi ERROR 'Got EXIT message from i-PI' is expected). "
                        "Continuing with trajectory verify and --order-csv if requested.",
                        flush=True,
                    )
                else:
                    raise
            if ipi_proc is not None:
                rc = ipi_proc.poll()
                if rc is not None and rc != 0:
                    print(
                        f"ERROR: i-PI subprocess exited with code {rc} before or while "
                        "the LAMMPS client finished. See i-PI log.",
                        file=sys.stderr,
                        flush=True,
                    )
                    sys.exit(1)
        else:
            _run_client_python(
                init_xyz=init_xyz,
                uma_driver_path=uma_driver_path,
                port=port,
                model=args.model,
                device=args.device,
                task_name=args.task_name,
                ipi_env=ipi_env,
                ipi_log=args.ipi_log,
            )

        _verify_traj_matches_requested_steps(
            IPI_DIR,
            steps,
            IPI_TRAJ_OUTPUT_STRIDE,
            fatal_if_short=(args.client == "lammps"),
        )

        if args.order_csv is not None:
            if args.order_no_wrap:
                order_pbc = "none"
            elif args.order_molecular_wrap:
                order_pbc = "molecular"
            else:
                order_pbc = "atom"
            _run_ipi_order_csv_postprocess(
                script_dir=_SCRIPT_DIR,
                ipi_dir=IPI_DIR,
                beads=args.beads,
                dt_fs=float(args.dt_fs),
                order_csv=args.order_csv,
                centroid=bool(args.order_centroid),
                pbc_wrap=order_pbc,
                thermo_prefer_ring_temperature=bool(args.order_thermo_prefer_ring_temperature),
            )

    finally:
        ipi_holder[0] = None
        if _ipi_signal_handlers_installed and prev_sigint is not None and prev_sigterm is not None:
            try:
                signal.signal(signal.SIGINT, prev_sigint)
                signal.signal(signal.SIGTERM, prev_sigterm)
            except ValueError:
                pass
        if ipi_proc is not None:
            if ipi_proc.poll() is None:
                try:
                    ipi_proc.wait(timeout=5.0)
                except subprocess.TimeoutExpired:
                    pass
            if ipi_proc.poll() is None:
                ipi_proc.terminate()
                try:
                    ipi_proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    ipi_proc.kill()
                    ipi_proc.wait(timeout=5)
            print("i-PI stopped")
        ipi_log_f.close()

    # Report outputs
    for pattern in ("ice_*traj*.xyz", "ice_*out*", "ice__i-pi.md"):
        for p in sorted(IPI_DIR.glob(pattern)):
            print(f"Output: {p}")


if __name__ == "__main__":
    main()

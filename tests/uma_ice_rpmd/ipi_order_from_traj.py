#!/usr/bin/env python3
"""
Compute Q6 order parameter from i-PI trajectory (centroid or single-bead xyz).

i-PI xyz format per frame:
  Line 1: natoms
  Line 2: # [CELL(abcABC): a b c alpha beta gamma] or other comment
  Lines 3+: symbol x y z (O and H)

Trajectory comments often include ``positions{atomic_unit}`` and ``cell{atomic_unit}``:
coordinates are then in **Bohr**; ``CELL(abcABC)`` lengths are also **Bohr** and are
converted to Å for PBC (MIC / wrap). If only ``cell{angstrom}`` appears (e.g. init.xyz
from :file:`convert_lammps_to_ipi_xyz.py`), a,b,c are already Å. Without a ``cell{…}``
tag, a,b,c are treated as Å (legacy i-PI quirk where numbers were Å despite labels).

By default each frame uses **per-atom** Cartesian wrapping into the primary
orthorhombic cell (same rule as ``ice_order_parameters.wrap_cartesian_orthorhombic``),
matching **LAMMPS / Fairchem** ``wrap_positions`` and OpenMM parity workflows.
i-PI xyz is usually **unwrapped**; without wrapping, MIC-based neighbor shells for
order parameters are unreliable. Use ``--molecular-wrap`` for the legacy O,H,H
group translation (one ``n*L`` shift per water). Use ``--no-wrap`` only for
coordinates already in ``[0,L)``.

When --beads > 1, reads traj_0 .. traj_{N-1} and computes **path-integral** order:
for each oxygen, average Q6 / q_tet over beads, then report mean/std/percentiles
over oxygens (same as OpenMM RPMD order CSV). Derives bead paths from --traj
(e.g. ipi/ice__i-pi.traj_0.xyz -> ipi/ice__i-pi.traj_{0..N-1}.xyz).

**OpenMM RPMD comparison:** omit ``--centroid`` so this matches the path-integral
estimator (observable averaged over beads). ``--centroid`` averages Cartesian
coordinates per atom after MIC alignment of each H₂O to bead 0; use it only when
you explicitly want a centroid structure.

Optional ``--thermo`` (default: use ``ipi/ice__i-pi.md`` if present): merge
``PE_kj_mol`` from i-PI properties.  **Centroid kinetic ``T_K``** matches OpenMM
when bead velocity trajectories exist (``ice__i-pi.vtraj_*.xyz`` from the
``<trajectory …>velocities</trajectory>`` line in ``input.xml``): same formula as
``rpmd_thermo_utils.centroid_kinetic_energy_and_temperature`` (mean bead velocity
per atom, then ½m|v_centroid|² and equipartition). If ``vtraj`` files are missing,
``T_K`` falls back to i-PI ``<properties>`` (``temperature(nm=0)`` then
``temperature``; use ``--thermo-prefer-ring-temperature`` to swap that order when
nm=0 is pathological).

Output CSV matches lammps_order_from_dump schema for plot_rpmd_comparison.py.

Usage:
  cd tests/uma_ice_rpmd
  python ipi_order_from_traj.py --traj ipi/ice_traj.xyz -o pipeline_out/ice_order_ipi_rpmd.csv
  python ipi_order_from_traj.py --traj ipi/ice__i-pi.traj_0.xyz --beads 4 -o pipeline_out/ice_order_ipi_rpmd.csv
  # Path-integral (default): same bead-averaging philosophy as OpenMM RPMD order CSV
  python ipi_order_from_traj.py --traj ipi/ice__i-pi.traj_00.xyz --beads 32 -o pipeline_out/ice_order_ipi_rpmd.csv
  # Centroid structure (MIC-aligned waters before mean), then Q6/q_tet once per frame
  python ipi_order_from_traj.py --traj ipi/ice__i-pi.traj_00.xyz --beads 32 --centroid -o pipeline_out/ice_order_ipi_centroid.csv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from ice_order_parameters import (
    ICE_ORDER_PI_CSV_HEADER,
    format_ice_order_pi_csv_row,
    ice_order_metrics_centroid_beads,
    ice_order_metrics_path_integral,
    wrap_cartesian_orthorhombic,
    wrap_water_molecules_orthorhombic,
)
from ipi_thermo_utils import load_thermo_lookups, lookup_thermo
from rpmd_thermo_utils import centroid_kinetic_energy_and_temperature


def _abc_to_box_orth(a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Convert a,b,c,angles (deg) to (Lx,Ly,Lz) for orthorhombic or triclinic.
    For orthogonal (alpha=beta=gamma=90): Lx=a, Ly=b, Lz=c.
    """
    if abs(alpha - 90) < 0.1 and abs(beta - 90) < 0.1 and abs(gamma - 90) < 0.1:
        return np.array([a, b, c], dtype=np.float64)
    # General triclinic: approximate orthorhombic box lengths
    rad = np.pi / 180.0
    cos_a = np.cos(alpha * rad)
    cos_b = np.cos(beta * rad)
    cos_g = np.cos(gamma * rad)
    vol = a * b * c * np.sqrt(
        1 - cos_a**2 - cos_b**2 - cos_g**2 + 2 * cos_a * cos_b * cos_g
    )
    # Use lengths as rough Lx,Ly,Lz (for orthorhombic projection)
    return np.array([a, b, c], dtype=np.float64)


BOHR_TO_ANG = 0.529177249  # 1 Bohr in Angstrom
# CODATA: atomic unit of velocity ℏ/(m_e a_0) in m/s; OpenMM reports nm/ps → scale ×10⁻³
AU_VELOCITY_M_PER_S = 2.18769126364e6
AU_VELOCITY_TO_NM_PS = AU_VELOCITY_M_PER_S * 1e-3

_CELL_TAG_RE = re.compile(r"cell\{([A-Za-z_]+)\}")
_POSITIONS_ATOMIC_RE = re.compile(r"positions\{atomic_unit\}")
_VELOCITIES_ATOMIC_RE = re.compile(r"velocities\{atomic_unit\}")


def _cell_abc_to_angstrom(comment: str, a: float, b: float, c: float) -> tuple[float, float, float]:
    """Interpret CELL(abcABC) lengths using ``cell{unit}`` in the comment line."""
    m = _CELL_TAG_RE.search(comment)
    if m:
        tag = m.group(1).lower()
        if tag in ("atomic_unit", "atomicunit", "au"):
            return (a * BOHR_TO_ANG, b * BOHR_TO_ANG, c * BOHR_TO_ANG)
        return (a, b, c)
    # No cell{…} tag: numeric abc are treated as Å (common for hand-written init.xyz)
    return (a, b, c)


def _positions_are_bohr(comment: str) -> bool:
    """True when ``positions{atomic_unit}`` is declared."""
    if _POSITIONS_ATOMIC_RE.search(comment):
        return True
    # Legacy: some lines only say atomic_unit without explicit positions{…}
    if "positions{" not in comment and "atomic_unit" in comment.lower():
        return True
    return False


def _velocities_are_atomic(comment: str) -> bool:
    """True when i-PI declares Cartesian velocities in atomic units (convert to nm/ps)."""
    if _VELOCITIES_ATOMIC_RE.search(comment):
        return True
    if "velocities{" not in comment and "atomic_unit" in comment.lower():
        return True
    return False


def _pos_traj_to_vel_traj_path(traj_path: Path) -> Path:
    """``ice__i-pi.traj_00.xyz`` → ``ice__i-pi.vtraj_00.xyz`` (i-PI ``filename=vtraj``)."""
    return traj_path.with_name(traj_path.name.replace(".traj_", ".vtraj_"))


def _masses_da_from_symbols(symbols: list[str]) -> np.ndarray:
    """Dalton masses for O/H (same order as trajectory); unknown symbols → 1.0 Da."""
    masses: list[float] = []
    for sym in symbols:
        u = sym.upper()
        if u == "O":
            masses.append(15.99943)
        elif u == "H":
            masses.append(1.00794)
        else:
            masses.append(1.0)
    return np.array(masses, dtype=np.float64)


def _masses_da_from_z(z: np.ndarray) -> np.ndarray:
    """Match OpenMM convention: O→15.99943 Da, H→1.00794 Da (``z`` as in ice_order_parameters)."""
    m = np.empty(z.shape[0], dtype=np.float64)
    for i, zi in enumerate(z):
        m[i] = 15.99943 if int(zi) == 8 else 1.00794
    return m


def _traj_paths_from_prefix(traj_path: Path, beads: int) -> list[Path]:
    """Infer bead trajectory paths from one bead file (e.g. ice__i-pi.traj_0 or traj_00).

    i-PI 3 may zero-pad indices (traj_00 .. traj_31). We detect padding width from the
    given path and generate bead_0 .. bead_{beads-1} with the same width.
    """
    parent = traj_path.parent
    suffix = traj_path.suffix
    stem = traj_path.stem
    m = re.match(r"^(.*_)(\d+)$", stem)
    if m:
        prefix_part, idx_str = m.group(1), m.group(2)
        width = len(idx_str)
        return [parent / f"{prefix_part}{b:0{width}d}{suffix}" for b in range(beads)]
    base = re.sub(r"_\d+$", "", stem)
    return [parent / f"{base}_{b}{suffix}" for b in range(beads)]


def _iter_ipi_xyz_frames(traj_path: Path, dt_fs: float = 0.1):
    """Yield (step, time_ps, pos_ang, box_orth, z) per frame.
    step/time parsed from comment if present; else frame index.
    If the comment has ``Step: N`` (i-PI), ``time_ps`` is ``N * dt_fs / 1000`` unless
    a ``time{…}`` field overrides. *dt_fs* must match the i-PI MD timestep (often 0.1 fs).
    Positions and cell lengths in Bohr are converted to Å per comment tags.
    """
    cell_re = re.compile(
        r"#\s*CELL[\(\{\[]?abcABC[\)\}\]]?:\s*([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)"
        r"\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)"
    )
    step_re = re.compile(r"Step:\s*(\d+)")
    time_re = re.compile(r"time\{\w+\}:\s*([\d.eE+-]+)")

    text = traj_path.read_text()
    lines = text.splitlines()
    i = 0
    frame_idx = 0
    while i < len(lines):
        try:
            natoms = int(lines[i])
        except (ValueError, IndexError):
            i += 1
            continue
        i += 1
        if i >= len(lines):
            break
        comment = lines[i]
        i += 1

        step = frame_idx
        m_step = step_re.search(comment)
        if m_step:
            step = int(m_step.group(1))
        # Fallback: one output frame per dt (wrong if traj stride > 1 and no Step: line)
        time_ps = float(frame_idx) * float(dt_fs) / 1000.0
        m = time_re.search(comment)
        if m:
            time_ps = float(m.group(1))
        elif m_step:
            time_ps = float(step) * float(dt_fs) / 1000.0

        box_orth = np.array([15.6, 14.7, 18.1], dtype=np.float64)  # fallback from init.xyz
        m = cell_re.search(comment)
        if m:
            a, b, c = float(m.group(1)), float(m.group(2)), float(m.group(3))
            alpha = float(m.group(4)) if len(m.groups()) >= 6 else 90.0
            beta = float(m.group(5)) if len(m.groups()) >= 6 else 90.0
            gamma = float(m.group(6)) if len(m.groups()) >= 6 else 90.0
            a, b, c = _cell_abc_to_angstrom(comment, a, b, c)
            box_orth = _abc_to_box_orth(a, b, c, alpha, beta, gamma)

        pos_bohr = _positions_are_bohr(comment)

        pos_list = []
        z_list = []
        for _ in range(natoms):
            if i >= len(lines):
                break
            parts = lines[i].split()
            i += 1
            if len(parts) >= 4:
                sym = parts[0]
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                pos_list.append([x, y, z])
                z_list.append(8 if sym.upper() == "O" else 1)

        if len(pos_list) == natoms:
            pos_ang = np.array(pos_list, dtype=np.float64)
            if pos_bohr:
                pos_ang = pos_ang * BOHR_TO_ANG
            z = np.array(z_list, dtype=np.int64)
            yield step, time_ps, pos_ang, box_orth, z
        frame_idx += 1


def _iter_ipi_xyz_velocity_frames(traj_path: Path, dt_fs: float = 0.1):
    """Yield (step, time_ps, vel_nmps, z) per frame; velocities converted to nm/ps for OpenMM parity."""
    cell_re = re.compile(
        r"#\s*CELL[\(\{\[]?abcABC[\)\}\]]?:\s*([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)"
        r"\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)"
    )
    step_re = re.compile(r"Step:\s*(\d+)")
    time_re = re.compile(r"time\{\w+\}:\s*([\d.eE+-]+)")

    text = traj_path.read_text()
    lines = text.splitlines()
    i = 0
    frame_idx = 0
    while i < len(lines):
        try:
            natoms = int(lines[i])
        except (ValueError, IndexError):
            i += 1
            continue
        i += 1
        if i >= len(lines):
            break
        comment = lines[i]
        i += 1

        step = frame_idx
        m_step = step_re.search(comment)
        if m_step:
            step = int(m_step.group(1))
        time_ps = float(frame_idx) * float(dt_fs) / 1000.0
        m = time_re.search(comment)
        if m:
            time_ps = float(m.group(1))
        elif m_step:
            time_ps = float(step) * float(dt_fs) / 1000.0

        vel_atomic = _velocities_are_atomic(comment)

        vel_list = []
        z_list = []
        for _ in range(natoms):
            if i >= len(lines):
                break
            parts = lines[i].split()
            i += 1
            if len(parts) >= 4:
                sym = parts[0]
                vx, vy, vz = float(parts[1]), float(parts[2]), float(parts[3])
                vel_list.append([vx, vy, vz])
                z_list.append(8 if sym.upper() == "O" else 1)

        if len(vel_list) == natoms:
            vel_nmps = np.array(vel_list, dtype=np.float64)
            if vel_atomic:
                vel_nmps = vel_nmps * AU_VELOCITY_TO_NM_PS
            z = np.array(z_list, dtype=np.int64)
            yield step, time_ps, vel_nmps, z
        frame_idx += 1


def _iter_ipi_order_frames(
    pos_traj_path: Path,
    beads: int,
    dt_fs: float,
    *,
    use_centroid_vel: bool,
):
    """Yield (step, time_ps, bead_poses_ang, box_orth, z, vel_stack_or_none).

    *vel_stack_or_none* is shape ``(n_beads, n_atoms, 3)`` in nm/ps when centroid
    velocities are available; else ``None`` (caller uses thermo for ``T_K``).
    """
    pos_paths = (
        [pos_traj_path]
        if beads <= 1
        else _traj_paths_from_prefix(pos_traj_path, beads)
    )
    for p in pos_paths:
        if not p.is_file():
            raise FileNotFoundError(f"Position trajectory not found: {p}")
    pos_iters = [_iter_ipi_xyz_frames(p, dt_fs=dt_fs) for p in pos_paths]

    vel_iters: list | None = None
    if use_centroid_vel:
        vel_paths = [_pos_traj_to_vel_traj_path(p) for p in pos_paths]
        if all(vp.is_file() for vp in vel_paths):
            vel_iters = [_iter_ipi_xyz_velocity_frames(vp, dt_fs=dt_fs) for vp in vel_paths]
        else:
            vel_iters = None

    warned_mismatch = False
    warned_vel_shorter = False
    while True:
        try:
            pos_rows = [next(pi) for pi in pos_iters]
        except StopIteration:
            break
        step = pos_rows[0][0]
        time_ps = pos_rows[0][1]
        bead_poses_ang = [r[2] for r in pos_rows]
        box_orth = pos_rows[0][3]
        z = pos_rows[0][4]

        vel_stack = None
        if vel_iters is not None:
            try:
                vel_rows = [next(vi) for vi in vel_iters]
            except StopIteration:
                vel_iters = None
                if not warned_vel_shorter:
                    print(
                        "Warning: velocity trajectory ended before positions; "
                        "remaining frames use T_K from thermo only.",
                        file=sys.stderr,
                    )
                    warned_vel_shorter = True
            else:
                for vr in vel_rows:
                    if vr[0] != step and not warned_mismatch:
                        print(
                            "Warning: position/velocity trajectory step mismatch; "
                            "centroid T_K may be misaligned.",
                            file=sys.stderr,
                        )
                        warned_mismatch = True
                        break
                vel_stack = np.stack([r[2] for r in vel_rows], axis=0)

        yield step, time_ps, bead_poses_ang, box_orth, z, vel_stack


def _iter_ipi_xyz_frames_beads(traj_path: Path, beads: int, dt_fs: float = 0.1):
    """Yield (step, time_ps, list[pos_ang per bead], box_orth, z) per frame."""
    paths = _traj_paths_from_prefix(traj_path, beads)
    for p in paths:
        if not p.is_file():
            raise FileNotFoundError(f"Bead trajectory not found: {p}")
    bead_frames: list[list[tuple[int, float, np.ndarray, np.ndarray, np.ndarray]]] = []
    for p in paths:
        bead_frames.append(list(_iter_ipi_xyz_frames(p, dt_fs=dt_fs)))
    n_frames = min(len(f) for f in bead_frames)
    if n_frames == 0:
        print("Warning: no frames parsed from bead trajectories", file=sys.stderr)
    for i in range(n_frames):
        step, time_ps, _foo, box_orth, z = bead_frames[0][i]
        poses = [bead_frames[b][i][2] for b in range(beads)]
        yield step, time_ps, poses, box_orth, z


def main() -> None:
    ap = argparse.ArgumentParser(description="i-PI trajectory → ice order CSV")
    ap.add_argument("--traj", type=Path, default=_SCRIPT_DIR / "ipi" / "ice_traj.xyz")
    ap.add_argument("-o", "--output", type=Path, default=None)
    ap.add_argument("--every", type=int, default=1, help="Write every N frames")
    ap.add_argument(
        "--dt-fs",
        type=float,
        default=0.1,
        help="i-PI MD timestep in fs; used as time_ps = step * dt_fs / 1000 when Step: is in the traj comment",
    )
    ap.add_argument(
        "--beads",
        type=int,
        default=1,
        help="RPMD beads; when >1, read traj_0..traj_{N-1} and path-integral average of Q6/q_tet",
    )
    ap.add_argument(
        "--centroid",
        action="store_true",
        help=(
            "Average Cartesian positions over all beads per frame, then compute Q6/q_tet on "
            "that centroid structure (default: average observable over beads, standard RPMD)."
        ),
    )
    ap.add_argument(
        "--no-wrap",
        action="store_true",
        help="Disable PBC wrapping: use raw Cartesian coords.",
    )
    ap.add_argument(
        "--molecular-wrap",
        action="store_true",
        help="Wrap each H2O as a rigid translation (O in [0,L)); legacy mode. Default is per-atom wrap.",
    )
    ap.add_argument(
        "--thermo",
        type=Path,
        default=None,
        metavar="PATH",
        help=(
            "i-PI properties file (e.g. ipi/ice__i-pi.md) for T_K and PE_kj_mol. "
            "If omitted, loads ipi/ice__i-pi.md when present (use --no-thermo to disable)."
        ),
    )
    ap.add_argument(
        "--no-thermo",
        action="store_true",
        help="Do not merge thermo from ipi/ice__i-pi.md even if it exists.",
    )
    ap.add_argument(
        "--thermo-prefer-ring-temperature",
        action="store_true",
        help=(
            "When bead velocity trajectories (vtraj) are missing, for merged T_K prefer "
            "i-PI ``temperature`` (ring kinetic T) over ``temperature(nm=0)``. "
            "Ignored for T_K when vtraj files are present (centroid T from velocities)."
        ),
    )
    args = ap.parse_args()
    if args.no_wrap and args.molecular_wrap:
        print("Error: use only one of --no-wrap and --molecular-wrap", file=sys.stderr)
        sys.exit(2)

    if args.no_wrap:
        pbc_wrap_mode = "none"
    elif args.molecular_wrap:
        pbc_wrap_mode = "molecular"
    else:
        pbc_wrap_mode = "atom"

    if not args.traj.is_file():
        print(f"Missing trajectory: {args.traj}", file=sys.stderr)
        sys.exit(1)

    out_path = args.output or args.traj.with_name(args.traj.stem + "_order.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    thermo_by_time: dict[float, tuple[str, str]] = {}
    thermo_by_step: dict[int, tuple[str, str]] = {}
    default_md = _SCRIPT_DIR / "ipi" / "ice__i-pi.md"
    if not args.no_thermo:
        md_path = args.thermo if args.thermo is not None else default_md
        if md_path.is_file():
            thermo_by_time, thermo_by_step = load_thermo_lookups(
                md_path,
                prefer_ring_temperature=bool(args.thermo_prefer_ring_temperature),
            )
        elif args.thermo is not None:
            print(f"Missing --thermo file: {md_path}", file=sys.stderr)
            sys.exit(1)

    pos_paths = (
        [args.traj]
        if args.beads <= 1
        else _traj_paths_from_prefix(args.traj, args.beads)
    )
    vel_paths_ok = all(_pos_traj_to_vel_traj_path(p).is_file() for p in pos_paths)
    if not vel_paths_ok:
        print(
            "Note: bead velocity files (*.vtraj_*) not found; T_K from i-PI .md properties "
            "(use run_ipi_lammps_uma_rpmd input.xml with <trajectory …>velocities</trajectory> for "
            "centroid T_K parity with OpenMM).",
            file=sys.stderr,
        )

    masses_template: np.ndarray | None = None
    ndof_centroid_ke = 0

    frame_iter = _iter_ipi_order_frames(
        args.traj,
        args.beads,
        args.dt_fs,
        use_centroid_vel=True,
    )

    with out_path.open("w", encoding="utf-8", newline="\n") as f:
        f.write(ICE_ORDER_PI_CSV_HEADER)
        f.flush()
        for frame_idx, (step, time_ps, bead_poses_ang, box_orth, z, vel_stack) in enumerate(
            frame_iter
        ):
            if frame_idx % args.every != 0:
                continue
            if masses_template is None:
                masses_template = _masses_da_from_z(z)
                ndof_centroid_ke = int(3 * np.count_nonzero(masses_template > 0))
            if pbc_wrap_mode == "atom":
                bead_poses_ang = [
                    wrap_cartesian_orthorhombic(p, box_orth) for p in bead_poses_ang
                ]
            elif pbc_wrap_mode == "molecular":
                bead_poses_ang = [
                    wrap_water_molecules_orthorhombic(p, box_orth) for p in bead_poses_ang
                ]
            if args.centroid:
                res = ice_order_metrics_centroid_beads(bead_poses_ang, box_orth, z=z)
            else:
                res = ice_order_metrics_path_integral(bead_poses_ang, box_orth, z=z)
            tk_s, pe_s = lookup_thermo(
                time_ps=time_ps,
                step=int(step),
                by_time=thermo_by_time,
                by_step=thermo_by_step,
            )
            if vel_stack is not None and masses_template is not None and ndof_centroid_ke > 0:
                _, centroid_temp = centroid_kinetic_energy_and_temperature(
                    vel_stack,
                    masses_template,
                    ndof_centroid_ke=ndof_centroid_ke,
                )
                tk_s = f"{centroid_temp:.4f}"
            f.write(format_ice_order_pi_csv_row(step, time_ps, tk_s, pe_s, res))
            f.flush()
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()

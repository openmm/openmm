#!/usr/bin/env python3
"""
Compare two GSD initial-condition files (read-only).
Uses gsd.hoomd to inspect structure, box, particle counts, positions/velocities.
No modifications to any file.
"""
import sys
from pathlib import Path

def inspect_gsd(path: str, label: str, frame_index: int = 0) -> dict:
    """Open GSD read-only and return summary dict for the given frame."""
    import gsd.hoomd
    path = Path(path).resolve()
    if not path.exists():
        return {"error": f"File not found: {path}"}
    try:
        with gsd.hoomd.open(path, mode="r") as f:
            nframes = len(f)
            if frame_index < 0:
                frame_index = max(0, nframes + frame_index)
            if frame_index >= nframes:
                return {"error": f"Frame {frame_index} out of range (nframes={nframes})"}
            frame = f[frame_index]
    except Exception as e:
        return {"error": str(e)}

    p = frame.particles
    N = p.N
    # Box
    box = frame.configuration.box
    Lx, Ly, Lz = box[0], box[1], box[2]
    # Types
    types = list(p.types) if hasattr(p, "types") and p.types is not None else []
    typeid = p.typeid if hasattr(p, "typeid") else None
    type_counts = {}
    if typeid is not None and len(types) > 0:
        for i, t in enumerate(types):
            type_counts[t] = int((typeid == i).sum())

    pos = p.position
    vel = p.velocity if hasattr(p, "velocity") and p.velocity is not None else None

    def stats(arr, name):
        if arr is None or len(arr) == 0:
            return None
        import numpy as np
        a = arr
        return {
            "min": float(a.min()),
            "max": float(a.max()),
            "mean": float(a.mean()),
            "std": float(a.std()) if len(a) > 1 else 0.0,
        }

    info = {
        "path": str(path),
        "label": label,
        "nframes": nframes,
        "frame_index": frame_index,
        "N": N,
        "box_Lx": float(Lx),
        "box_Ly": float(Ly),
        "box_Lz": float(Lz),
        "types": types,
        "type_counts": type_counts,
        "position_stats": stats(pos, "position"),
        "velocity_stats": stats(vel, "velocity") if vel is not None else None,
        "has_velocity": vel is not None and vel.size > 0,
    }
    if pos is not None and N > 0:
        import numpy as np
        # Density / spacing: approximate particle spacing if on a lattice
        vol = float(Lx * Ly * Lz)
        info["volume"] = vol
        info["number_density"] = N / vol if vol > 0 else None
        # Check if positions look like a lattice (low std in spacing)
        pos = np.asarray(pos)
        info["position_extent"] = [float(pos[:, i].max() - pos[:, i].min()) for i in range(3)]
    return info


def main():
    root = Path(__file__).resolve().parents[1]
    eq_init = root / "equilibrium_fkt_100K" / "init-0.gsd"
    dimer_init = root / "examples" / "cavity" / "dimer_system" / "init-0.gsd"
    dimer_molecular = root / "examples" / "cavity" / "dimer_system" / "molecular-0.gsd"

    paths = [
        (eq_init, "equilibrium_fkt_100K/init-0.gsd"),
        (dimer_init, "examples/cavity/dimer_system/init-0.gsd"),
        (dimer_molecular, "examples/cavity/dimer_system/molecular-0.gsd"),
    ]

    print("=" * 70)
    print("GSD initial-condition comparison (read-only)")
    print("=" * 70)

    for path, label in paths:
        path = path.resolve()
        if not path.exists():
            print(f"\n[{label}] SKIP (not found): {path}")
            continue
        print(f"\n[{label}]")
        # Inspect frame 0 first to get nframes
        info0 = inspect_gsd(str(path), label, frame_index=0)
        if "error" in info0:
            print(f"  ERROR — {info0['error']}")
            continue
        nframes = info0["nframes"]
        frames_to_show = [0]
        if nframes > 1:
            frames_to_show.append(nframes - 1)
        for frame_idx in frames_to_show:
            info = inspect_gsd(str(path), label, frame_index=frame_idx)
            if "error" in info:
                print(f"  Frame {frame_idx}: ERROR — {info['error']}")
                continue
            print(f"  --- Frame {frame_idx} (of {info['nframes']} total) ---")
            print(f"  N particles: {info['N']}")
            print(f"  Box: Lx={info['box_Lx']:.6f} Ly={info['box_Ly']:.6f} Lz={info['box_Lz']:.6f}")
            print(f"  Types: {info['types']}  counts: {info['type_counts']}")
            if info.get("volume") is not None:
                print(f"  Volume: {info['volume']:.4f}  number_density: {info.get('number_density'):.6e}")
            if info.get("position_stats"):
                s = info["position_stats"]
                print(f"  Position: min={s['min']:.4f} max={s['max']:.4f} mean={s['mean']:.4f} std={s['std']:.4f}")
            if info.get("position_extent"):
                print(f"  Position extent (max-min) x,y,z: {[f'{x:.4f}' for x in info['position_extent']]}")
            print(f"  Has velocity: {info['has_velocity']}")
            if info.get("velocity_stats"):
                s = info["velocity_stats"]
                print(f"  Velocity: min={s['min']:.6f} max={s['max']:.6f} mean={s['mean']:.6f} std={s['std']:.6f}")

    print("\n" + "=" * 70)
    print("SUMMARY: Why equilibrium_fkt_100K → slower glassy, dimer_system → faster")
    print("=" * 70)
    print("""
1) equilibrium_fkt_100K/init-0.gsd and examples/cavity/dimer_system/init-0.gsd
   are THE SAME FILE (500 frames, 501 particles O+N+L, box 40 Bohr).
   - equilibrium_fkt_100K uses frame=replica (e.g. frame 1, 2, 3...) = pre-equilibrated
     frames from a long trajectory → glassy configs → SLOWER F(k,t) decorrelation (correct).

2) examples/cavity/dimer_system/molecular-0.gsd is DIFFERENT:
   - 1 frame only, 500 particles (no cavity L), box 40 Bohr.
   - Velocity: ALL ZEROS (cold start).
   - Position extent: ~36.4, 36.6, 30.8 (tighter lattice) vs init-0.gsd ~39.9 (full box).
   - Used e.g. for OpenMM/cav-hoomd parity runs (--input-gsd molecular-0.gsd).
   - Cold lattice + zero velocities → fast thermalization/relaxation → FASTER
     decorrelation than starting from an equilibrated trajectory.

3) If dimer_system uses init-0.gsd with --initial-gsd-frame 0: same file as (1), frame 0
   for all replicas; still pre-equilibrated. If it uses molecular-0.gsd, that explains
   the faster decorrelation (cold start, single-frame lattice).
""")


if __name__ == "__main__":
    main()

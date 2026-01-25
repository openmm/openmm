#!/usr/bin/env python3
"""
Overlay cavity position plots on rendered protein frames.

This script composes a top matplotlib plot (mode energy, interaction, self polarization)
with each rendered frame, producing a new set of frames and an mp4.

Inputs:
  - Rendered frame images (PNG/TGA) from PyMOL/VMD
  - Cavity data (.npz or structured .npy)

Usage:
  python overlay_cavity_plot.py \
      --trajectory-dir trajectory \
      --frames-dir movie_frames \
      --cavity-data protein_cavity_lambda0.0100.npz \
      --output protein_movie_with_cavity.mp4
"""

import argparse
import re
import subprocess
from pathlib import Path

import numpy as np
from PIL import Image

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_time_ps_from_name(path):
    match = re.search(r"_t([0-9.]+)ps", Path(path).name)
    if match:
        return float(match.group(1))
    return None


def load_cavity_data(path):
    path = Path(path)
    if path.suffix == ".npy":
        stream = np.load(path)
        time_ps = stream["time_ps"]
        dipole_nm = stream["dipole_nm"]
        cavity_nm = stream["cavity_nm"]
        return time_ps, dipole_nm, cavity_nm
    if path.suffix == ".npz":
        with np.load(path, allow_pickle=True) as data:
            time_ps = data["time_ps"]
            dipole_nm = data["dipole_nm"]
            cavity_nm = data["cavity_nm"]
        return time_ps, dipole_nm, cavity_nm
    raise ValueError(f"Unsupported cavity data format: {path}")


def build_plot_image(
    time_ps,
    cavity_nm,
    current_time_ps,
    width_px,
    plot_height_px,
):
    dpi = 100
    fig_w = width_px / dpi
    fig_h = plot_height_px / dpi
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)

    ax.plot(time_ps, cavity_nm[:, 0], label="q_x", color="#1f77b4", linewidth=1.2)
    ax.plot(time_ps, cavity_nm[:, 1], label="q_y", color="#ff7f0e", linewidth=1.2)
    ax.plot(time_ps, cavity_nm[:, 2], label="q_z", color="#2ca02c", linewidth=1.2)

    ax.axvline(current_time_ps, color="black", linestyle="--", linewidth=1.0)

    current_x = np.interp(current_time_ps, time_ps, cavity_nm[:, 0])
    current_y = np.interp(current_time_ps, time_ps, cavity_nm[:, 1])
    current_z = np.interp(current_time_ps, time_ps, cavity_nm[:, 2])
    ax.scatter([current_time_ps] * 3, [current_x, current_y, current_z], color="black", s=10)

    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Cavity position (nm)")
    ax.legend(loc="upper right", frameon=False, ncol=3, fontsize=8)
    ax.grid(alpha=0.2)

    fig.tight_layout(pad=0.6)

    fig.canvas.draw()
    # Matplotlib >= 3.8: use buffer_rgba()
    rgba = np.asarray(fig.canvas.buffer_rgba())
    plot_img = rgba[:, :, :3]
    plt.close(fig)

    return Image.fromarray(plot_img)


def main():
    parser = argparse.ArgumentParser(
        description="Overlay cavity energy plot on protein frames",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--trajectory-dir", type=str, default="trajectory")
    parser.add_argument("--frames-dir", type=str, default="movie_frames")
    parser.add_argument("--frames-glob", type=str, default="frame_*.png")
    parser.add_argument("--cavity-data", type=str, required=True)
    parser.add_argument("--metadata-json", type=str, default=None)
    parser.add_argument("--output", type=str, default="protein_movie_with_cavity.mp4")
    parser.add_argument("--output-frames-dir", type=str, default="movie_frames_with_plot")
    parser.add_argument("--plot-height", type=int, default=280)
    parser.add_argument("--fps", type=int, default=24)
    parser.add_argument("--start-frame", type=int, default=0)
    parser.add_argument("--end-frame", type=int, default=None)
    parser.add_argument("--skip-compile", action="store_true")
    args = parser.parse_args()

    traj_dir = Path(args.trajectory_dir)
    frames_dir = Path(args.frames_dir)
    output_frames_dir = Path(args.output_frames_dir)
    output_frames_dir.mkdir(exist_ok=True)

    frame_images = sorted(frames_dir.glob(args.frames_glob))
    if not frame_images:
        raise RuntimeError(f"No frames found in {frames_dir} with pattern {args.frames_glob}")

    traj_frames = sorted(traj_dir.glob("frame_*.pdb"))
    frame_times = [parse_time_ps_from_name(p) for p in traj_frames]
    have_times = all(t is not None for t in frame_times) and len(frame_times) == len(frame_images)

    time_ps, dipole_nm, cavity_nm = load_cavity_data(args.cavity_data)

    if args.end_frame is None or args.end_frame >= len(frame_images):
        end_frame = len(frame_images) - 1
    else:
        end_frame = args.end_frame
    if args.start_frame < 0 or end_frame < args.start_frame:
        raise ValueError("Invalid start/end frame range")

    # Determine time per frame
    if have_times:
        frame_time_ps = frame_times
    else:
        frame_time_ps = np.linspace(time_ps[0], time_ps[-1], len(frame_images))

    # Compose frames
    for idx in range(args.start_frame, end_frame + 1):
        frame_path = frame_images[idx]
        frame_img = Image.open(frame_path).convert("RGB")
        width_px, height_px = frame_img.size

        plot_img = build_plot_image(
            time_ps=time_ps,
            cavity_nm=cavity_nm,
            current_time_ps=frame_time_ps[idx],
            width_px=width_px,
            plot_height_px=args.plot_height,
        )

        combined = Image.new("RGB", (width_px, height_px + args.plot_height), (255, 255, 255))
        combined.paste(plot_img, (0, 0))
        combined.paste(frame_img, (0, args.plot_height))

        out_path = output_frames_dir / f"{frame_path.stem}.png"
        combined.save(out_path)

        if (idx + 1) % 10 == 0 or idx == end_frame:
            print(f"Rendered overlay {idx + 1 - args.start_frame}/{end_frame - args.start_frame + 1}")

    if args.skip_compile:
        return

    # Compile with ffmpeg
    ffmpeg_cmd = [
        "ffmpeg", "-y",
        "-framerate", str(args.fps),
        "-start_number", str(args.start_frame),
        "-i", str(output_frames_dir / "frame_%06d.png"),
        "-frames:v", str(end_frame - args.start_frame + 1),
        "-c:v", "libx264",
        "-preset", "slow",
        "-crf", "18",
        "-pix_fmt", "yuv420p",
        str(Path(args.output).resolve()),
    ]
    subprocess.run(ffmpeg_cmd, check=True)
    print(f"Movie saved: {Path(args.output).resolve()}")


if __name__ == "__main__":
    main()

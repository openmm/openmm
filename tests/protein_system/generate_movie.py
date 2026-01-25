#!/usr/bin/env python3
"""
Generate Protein-Water Movie from MD Trajectory
================================================

This script creates a movie visualization of a protein simulation with:
- Protein shown as opaque surface (colored by secondary structure or chain)
- First solvation shell waters shown as transparent spheres/sticks
- Smooth camera rotation (optional)

Requirements:
    - pymol (conda install -c conda-forge pymol-open-source)
    - ffmpeg (for video compilation)
    - MDAnalysis (optional, for trajectory preprocessing)

Usage:
    python generate_movie.py --trajectory-dir trajectory --output movie.mp4
    python generate_movie.py --help
"""

import os
import sys
import argparse
import subprocess
import glob
from pathlib import Path
import shutil

# Try to import PyMOL
try:
    import pymol
    from pymol import cmd
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False
    print("Warning: PyMOL not available. Will generate VMD script instead.")


def find_trajectory_frames(trajectory_dir):
    """Find all PDB frame files in the trajectory directory."""
    traj_path = Path(trajectory_dir)
    
    # Look for frame files
    frame_files = sorted(glob.glob(str(traj_path / "frame_*.pdb")))
    
    if not frame_files:
        raise FileNotFoundError(f"No frame_*.pdb files found in {trajectory_dir}")
    
    # Also look for topology
    topology = traj_path / "topology.pdb"
    if not topology.exists():
        # Use first frame as topology
        topology = Path(frame_files[0])
    
    print(f"Found {len(frame_files)} frames")
    print(f"Topology: {topology}")
    
    return frame_files, topology


def generate_pymol_movie(
    frame_files,
    topology,
    output_path,
    solvation_cutoff=5.0,
    protein_color="spectrum",
    water_transparency=0.7,
    resolution=(1920, 1080),
    fps=24,
    rotation_per_frame=1.0,
    ray_trace=True,
    temp_dir="movie_frames",
    start_frame=0,
    end_frame=None,
    skip_render=False,
    skip_compile=False,
    keep_frames=False,
):
    """
    Generate movie using PyMOL.
    
    Parameters
    ----------
    frame_files : list
        List of PDB frame file paths
    topology : Path
        Path to topology PDB
    output_path : str
        Output movie file path
    solvation_cutoff : float
        Distance cutoff for first solvation shell (Angstroms)
    protein_color : str
        Color scheme for protein ("spectrum", "ss", "chain", "white")
    water_transparency : float
        Transparency of water molecules (0-1, higher = more transparent)
    resolution : tuple
        Output resolution (width, height)
    fps : int
        Frames per second in output movie
    rotation_per_frame : float
        Degrees of rotation per frame (0 for no rotation)
    ray_trace : bool
        Use ray tracing for higher quality (slower)
    temp_dir : str
        Directory for temporary frame images
    """
    # Create temp directory for frames
    temp_path = Path(temp_dir).resolve()
    temp_path.mkdir(exist_ok=True)
    output_path = Path(output_path).resolve()
    
    total_frames = len(frame_files)
    if end_frame is None or end_frame >= total_frames:
        end_frame = total_frames - 1
    if start_frame < 0 or end_frame < start_frame:
        raise ValueError("Invalid start/end frame range")

    if skip_render:
        if not temp_path.exists():
            raise RuntimeError(f"Frame directory not found: {temp_path}")
        if not (temp_path / f"frame_{start_frame:06d}.png").exists():
            raise RuntimeError(
                f"No frame images found for start frame {start_frame}. "
                f"Render frames first or set --start-frame to an existing frame."
            )
    else:
        if not PYMOL_AVAILABLE:
            raise RuntimeError("PyMOL is required for movie generation")
        
        # Initialize PyMOL in headless mode
        pymol.finish_launching(['pymol', '-cq'])
    
        # Set up rendering
        cmd.set("ray_trace_mode", 1)
        cmd.set("ray_shadows", 1)
        cmd.set("ray_trace_gain", 0.1)
        cmd.set("ambient", 0.4)
        cmd.set("direct", 0.6)
        cmd.set("specular", 0.3)
        cmd.set("shininess", 30)
        cmd.set("antialias", 2)
        cmd.set("orthoscopic", 1)
        
        # Background
        cmd.bg_color("white")
        
        # Load first frame to set up view
        print("Setting up initial view...")
        cmd.load(str(frame_files[0]), "frame")
    
        def select_first_valid(selection_name, expr_candidates):
            last_error = None
            for candidate in expr_candidates:
                try:
                    cmd.select(selection_name, candidate)
                    return candidate
                except Exception as exc:
                    last_error = exc
                    continue
            raise RuntimeError(
                f"Failed to create selection '{selection_name}'. "
                f"Tried: {', '.join(expr_candidates)}. "
                f"Last error: {last_error}"
            )

        # Select protein and waters with fallbacks for selection syntax
        protein_expr = select_first_valid(
            "protein_sel",
            [
                "polymer.protein",
                "polymer",
                "protein",
                "not resn HOH",
            ],
        )
        water_expr = select_first_valid(
            "water_sel",
            [
                f"resn HOH and within {solvation_cutoff} of (protein_sel)",
                f"resn HOH within {solvation_cutoff} of (protein_sel)",
                f"resn HOH and (around {solvation_cutoff} and protein_sel)",
                (
                    f"(resn HOH or resn WAT or resn TIP3 or resn TIP4) "
                    f"and within {solvation_cutoff} of (protein_sel)"
                ),
            ],
        )

        # Hide everything first
        cmd.hide("everything")
        
        # Show protein as surface
        cmd.show("surface", "protein_sel")
        
        # Color protein
        if protein_color == "spectrum":
            cmd.spectrum("count", "rainbow", "protein_sel")
        elif protein_color == "ss":
            cmd.color("red", "protein_sel and ss h")  # Helices
            cmd.color("yellow", "protein_sel and ss s")  # Sheets
            cmd.color("green", "protein_sel and ss l+''")  # Loops
        elif protein_color == "chain":
            cmd.spectrum("chain", "rainbow", "protein_sel")
        else:
            cmd.color(protein_color, "protein_sel")
        
        # Show water as small spheres
        cmd.show("spheres", "water_sel and name O")
        cmd.set("sphere_scale", 0.3, "water_sel")
        cmd.color("lightblue", "water_sel")
        cmd.set("sphere_transparency", water_transparency, "water_sel")
        
        # Set up view - zoom to protein with some padding
        cmd.zoom("protein_sel", buffer=5)
        
        # Store the view
        cmd.get_view()
        
        # Set viewport
        cmd.viewport(resolution[0], resolution[1])

    if not skip_render:
        print(f"Rendering frames {start_frame}..{end_frame} (total {end_frame - start_frame + 1})")
        
        for i, frame_file in enumerate(frame_files):
            if i < start_frame or i > end_frame:
                continue
            # Load frame (replace coordinates)
            if i > 0:
                cmd.delete("frame")
                cmd.load(str(frame_file), "frame")
                
                # Reselect waters for this frame (they may have moved)
                cmd.select("water_sel", water_expr)
                
                # Reapply representations
                cmd.hide("everything")
                cmd.show("surface", "protein_sel")
                cmd.show("spheres", "water_sel and name O")
                
                # Reapply colors
                if protein_color == "spectrum":
                    cmd.spectrum("count", "rainbow", "protein_sel")
                elif protein_color == "ss":
                    cmd.color("red", "protein_sel and ss h")
                    cmd.color("yellow", "protein_sel and ss s")
                    cmd.color("green", "protein_sel and ss l+''")
                elif protein_color == "chain":
                    cmd.spectrum("chain", "rainbow", "protein_sel")
                else:
                    cmd.color(protein_color, "protein_sel")
                
                cmd.color("lightblue", "water_sel")
                cmd.set("sphere_transparency", water_transparency, "water_sel")
            
            # Optional rotation
            if rotation_per_frame > 0:
                cmd.turn("y", rotation_per_frame)
            
            # Render frame
            frame_path = temp_path / f"frame_{i:06d}.png"
            if i == start_frame or i % 10 == 0 or i == end_frame:
                print(f"  Rendering frame {i} -> {frame_path}")
            if frame_path.exists():
                if i == start_frame or i % 10 == 0 or i == end_frame:
                    print(f"  Skipping existing frame {i}")
                continue
            
            if ray_trace:
                cmd.ray(resolution[0], resolution[1])
            
            cmd.png(str(frame_path), width=resolution[0], height=resolution[1], dpi=150, ray=0)
            
            # Progress
            if (i + 1) % 10 == 0 or i == end_frame:
                print(f"  Rendered {i + 1 - start_frame}/{end_frame - start_frame + 1} frames")
        
        # Close PyMOL
        cmd.quit()
    
    # Compile movie with ffmpeg
    if skip_compile:
        return
    print("\nCompiling movie with ffmpeg...")
    ffmpeg_cmd = [
        "ffmpeg", "-y",
        "-framerate", str(fps),
        "-start_number", str(start_frame),
        "-i", str(temp_path / "frame_%06d.png"),
    ]
    if end_frame is not None:
        ffmpeg_cmd.extend(["-frames:v", str(end_frame - start_frame + 1)])
    ffmpeg_cmd.extend([
        "-c:v", "libx264",
        "-preset", "slow",
        "-crf", "18",
        "-pix_fmt", "yuv420p",
        str(output_path),
    ])
    
    try:
        subprocess.run(ffmpeg_cmd, check=True, capture_output=True)
        print(f"Movie saved: {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"ffmpeg error: {e.stderr.decode()}")
        raise
    except FileNotFoundError:
        print("ffmpeg not found. Frame images saved in:", temp_path)
        print("Compile manually with: ffmpeg -framerate 24 -i frame_%06d.png -c:v libx264 -pix_fmt yuv420p movie.mp4")
        return
    
    # Clean up temp files
    if not keep_frames and not skip_render:
        print("Cleaning up temporary files...")
        shutil.rmtree(temp_path)
    
    print("Done!")


def generate_vmd_script(
    frame_files,
    topology,
    output_script,
    solvation_cutoff=5.0,
    resolution=(1920, 1080),
    rotation_per_frame=1.0,
    temp_dir="vmd_frames",
):
    """
    Generate a VMD Tcl script for movie rendering (alternative to PyMOL).
    """
    script = f'''# VMD Movie Generation Script
# Generated by generate_movie.py
#
# Usage: vmd -e {output_script}

# Load topology
mol new {topology} type pdb waitfor all

# Load trajectory frames
'''
    
    for frame in frame_files[1:]:
        script += f'mol addfile {frame} type pdb waitfor all\n'
    
    script += f'''
# Set up representations
mol delrep 0 top

# Protein surface
mol representation NewCartoon
mol color Structure
mol selection {{protein}}
mol material Opaque
mol addrep top

# Protein surface overlay
mol representation Surf 1.4 0.0 0.0
mol color Structure  
mol selection {{protein}}
mol material Transparent
mol addrep top

# First solvation shell waters
mol representation VDW 0.3 12.0
mol color ColorID 23
mol selection {{water and within {solvation_cutoff} of protein}}
mol material Transparent
mol addrep top

# Set up display
display projection Orthographic
display depthcue off
axes location Off
color Display Background white
display resize {resolution[0]} {resolution[1]}

# Center on protein
mol top 0
display resetview

# Render movie frames
set nframes [molinfo top get numframes]
set framedir "{temp_dir}"
file mkdir $framedir

for {{set i 0}} {{$i < $nframes}} {{incr i}} {{
    animate goto $i
    display update
    
    # Optional rotation
    rotate y by {rotation_per_frame}
    
    # Render
    set filename [format "$framedir/frame_%06d.tga" $i]
    render TachyonInternal $filename
    
    puts "Rendered frame $i / $nframes"
}}

puts "Done! Compile with: ffmpeg -framerate 24 -i vmd_frames/frame_%06d.tga -c:v libx264 -pix_fmt yuv420p movie.mp4"
quit
'''
    
    with open(output_script, 'w') as f:
        f.write(script)
    
    print(f"VMD script saved: {output_script}")
    print(f"Run with: vmd -e {output_script}")


def run_vmd_render(
    frame_files,
    topology,
    output_path,
    solvation_cutoff=5.0,
    resolution=(1920, 1080),
    fps=24,
    rotation_per_frame=1.0,
    temp_dir="vmd_frames",
    start_frame=0,
    end_frame=None,
    skip_render=False,
    skip_compile=False,
    keep_frames=False,
):
    """
    Render frames with VMD in headless mode and compile with ffmpeg.
    """
    total_frames = len(frame_files)
    if end_frame is None or end_frame >= total_frames:
        end_frame = total_frames - 1
    if start_frame < 0 or end_frame < start_frame:
        raise ValueError("Invalid start/end frame range")

    temp_path = Path(temp_dir).resolve()
    temp_path.mkdir(exist_ok=True)
    output_path = Path(output_path).resolve()

    if skip_render:
        if not temp_path.exists():
            raise RuntimeError(f"Frame directory not found: {temp_path}")
        if not (temp_path / f"frame_{start_frame:06d}.tga").exists():
            raise RuntimeError(
                f"No VMD frame images found for start frame {start_frame}. "
                f"Render frames first or set --temp-dir to the VMD frames directory."
            )
    else:
        script_path = output_path.with_suffix(".tcl")
        generate_vmd_script(
            frame_files=frame_files,
            topology=topology,
            output_script=str(script_path),
            solvation_cutoff=solvation_cutoff,
            resolution=resolution,
            rotation_per_frame=rotation_per_frame,
            temp_dir=temp_path,
        )
        # Run VMD in text mode (headless)
        vmd_cmd = ["vmd", "-dispdev", "text", "-e", str(script_path)]
        subprocess.run(vmd_cmd, check=True)

    if skip_compile:
        return

    # Compile movie with ffmpeg
    print("\nCompiling movie with ffmpeg...")
    ffmpeg_cmd = [
        "ffmpeg", "-y",
        "-framerate", str(fps),
        "-start_number", str(start_frame),
        "-i", str(temp_path / "frame_%06d.tga"),
    ]
    if end_frame is not None:
        ffmpeg_cmd.extend(["-frames:v", str(end_frame - start_frame + 1)])
    ffmpeg_cmd.extend([
        "-c:v", "libx264",
        "-preset", "slow",
        "-crf", "18",
        "-pix_fmt", "yuv420p",
        str(output_path),
    ])
    try:
        subprocess.run(ffmpeg_cmd, check=True, capture_output=True)
        print(f"Movie saved: {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"ffmpeg error: {e.stderr.decode()}")
        raise

    if not keep_frames and not skip_render:
        print("Cleaning up temporary files...")
        shutil.rmtree(temp_path)


def generate_matplotlib_movie(
    frame_files,
    output_path,
    solvation_cutoff=5.0,
    fps=24,
):
    """
    Generate a simple movie using matplotlib (fallback if PyMOL not available).
    This is lower quality but doesn't require external dependencies.
    """
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.animation as animation
    except ImportError:
        raise RuntimeError("matplotlib is required for this visualization method")
    
    try:
        import MDAnalysis as mda
    except ImportError:
        raise RuntimeError("MDAnalysis is required for trajectory loading")
    
    print("Loading trajectory with MDAnalysis...")
    
    # Load trajectory
    u = mda.Universe(str(frame_files[0]))
    for frame in frame_files[1:]:
        u.load_new(str(frame))
    
    # Select atoms
    protein = u.select_atoms("protein")
    
    print(f"Protein atoms: {len(protein)}")
    print(f"Total frames: {len(u.trajectory)}")
    
    # Set up figure
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    def update(frame_idx):
        ax.clear()
        
        # Go to frame
        u.trajectory[frame_idx]
        
        # Get protein CA positions for simplified view
        ca = u.select_atoms("protein and name CA")
        pos = ca.positions
        
        # Get nearby waters
        waters = u.select_atoms(f"(resname HOH or resname WAT or resname TIP3 or resname TIP4) and within {solvation_cutoff} of protein")
        water_pos = waters.select_atoms("name O*").positions if len(waters) > 0 else []
        
        # Plot protein backbone
        ax.plot(pos[:, 0], pos[:, 1], pos[:, 2], 'b-', linewidth=2, alpha=0.8)
        ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c='blue', s=20, alpha=0.6)
        
        # Plot waters
        if len(water_pos) > 0:
            ax.scatter(water_pos[:, 0], water_pos[:, 1], water_pos[:, 2], 
                      c='lightblue', s=10, alpha=0.3, marker='o')
        
        # Set labels
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        ax.set_title(f'Frame {frame_idx + 1}/{len(u.trajectory)}')
        
        # Rotate view
        ax.view_init(elev=20, azim=frame_idx * 2)
        
        return ax,
    
    print("Generating animation...")
    anim = animation.FuncAnimation(
        fig, update, frames=len(u.trajectory),
        interval=1000/fps, blit=False
    )
    
    # Save
    print(f"Saving movie to {output_path}...")
    writer = animation.FFMpegWriter(fps=fps, bitrate=5000)
    anim.save(output_path, writer=writer)
    
    plt.close()
    print("Done!")


def main():
    parser = argparse.ArgumentParser(
        description="Generate protein-water movie from MD trajectory",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--trajectory-dir", "-t", type=str, default="trajectory",
        help="Directory containing frame_*.pdb files"
    )
    parser.add_argument(
        "--output", "-o", type=str, default="protein_movie.mp4",
        help="Output movie file"
    )
    parser.add_argument(
        "--solvation-cutoff", type=float, default=5.0,
        help="Distance cutoff for first solvation shell (Angstroms)"
    )
    parser.add_argument(
        "--protein-color", type=str, default="spectrum",
        choices=["spectrum", "ss", "chain", "white", "cyan", "green"],
        help="Protein coloring scheme"
    )
    parser.add_argument(
        "--water-transparency", type=float, default=0.7,
        help="Water transparency (0=opaque, 1=invisible)"
    )
    parser.add_argument(
        "--resolution", type=str, default="1920x1080",
        help="Output resolution (WIDTHxHEIGHT)"
    )
    parser.add_argument(
        "--fps", type=int, default=24,
        help="Frames per second"
    )
    parser.add_argument(
        "--rotation", type=float, default=1.0,
        help="Degrees of rotation per frame (0 for static view)"
    )
    parser.add_argument(
        "--no-raytrace", action="store_true",
        help="Disable ray tracing (faster but lower quality)"
    )
    parser.add_argument(
        "--method", type=str, default="auto",
        choices=["auto", "pymol", "vmd", "matplotlib"],
        help="Visualization method"
    )
    parser.add_argument(
        "--temp-dir", type=str, default="movie_frames",
        help="Temporary directory for frame images"
    )
    parser.add_argument(
        "--start-frame", type=int, default=0,
        help="Start frame index (0-based)"
    )
    parser.add_argument(
        "--end-frame", type=int, default=None,
        help="End frame index (inclusive). Default renders to last frame."
    )
    parser.add_argument(
        "--skip-render", action="store_true",
        help="Skip rendering frames (compile only)"
    )
    parser.add_argument(
        "--skip-compile", action="store_true",
        help="Skip ffmpeg compilation (render only)"
    )
    parser.add_argument(
        "--keep-frames", action="store_true",
        help="Keep rendered frames after compilation"
    )
    
    args = parser.parse_args()
    
    # Parse resolution
    try:
        width, height = map(int, args.resolution.split('x'))
        resolution = (width, height)
    except ValueError:
        print(f"Invalid resolution format: {args.resolution}")
        sys.exit(1)
    
    # Find trajectory files
    try:
        frame_files, topology = find_trajectory_frames(args.trajectory_dir)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    # Choose method
    method = args.method
    if method == "auto":
        if PYMOL_AVAILABLE:
            method = "pymol"
        else:
            method = "vmd"
            print("PyMOL not available, generating VMD script instead")
    
    print(f"\nUsing method: {method}")
    print(f"Solvation cutoff: {args.solvation_cutoff} Å")
    print(f"Output: {args.output}")
    
    if method == "pymol":
        generate_pymol_movie(
            frame_files=frame_files,
            topology=topology,
            output_path=args.output,
            solvation_cutoff=args.solvation_cutoff,
            protein_color=args.protein_color,
            water_transparency=args.water_transparency,
            resolution=resolution,
            fps=args.fps,
            rotation_per_frame=args.rotation,
            ray_trace=not args.no_raytrace,
            temp_dir=args.temp_dir,
            start_frame=args.start_frame,
            end_frame=args.end_frame,
            skip_render=args.skip_render,
            skip_compile=args.skip_compile,
            keep_frames=args.keep_frames,
        )
    elif method == "vmd":
        run_vmd_render(
            frame_files=frame_files,
            topology=topology,
            output_path=args.output,
            solvation_cutoff=args.solvation_cutoff,
            resolution=resolution,
            fps=args.fps,
            rotation_per_frame=args.rotation,
            temp_dir=args.temp_dir,
            start_frame=args.start_frame,
            end_frame=args.end_frame,
            skip_render=args.skip_render,
            skip_compile=args.skip_compile,
            keep_frames=args.keep_frames,
        )
    elif method == "matplotlib":
        generate_matplotlib_movie(
            frame_files=frame_files,
            output_path=args.output,
            solvation_cutoff=args.solvation_cutoff,
            fps=args.fps,
        )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Cavity Molecular Dynamics Simulation of a Protein in Explicit Solvent
====================================================================

System: Human pepsin 3b (PDB: 3UTL) + explicit TIP4P-Ew water
Cavity: Resonant with OH stretch (3656-3663 cm^-1)
Protocol: Equilibrate -> Activate coupling -> Observe dynamics

This script mirrors the water-case workflow and output format, but builds
the system by downloading a protein from RCSB and preparing it with PDBFixer.
"""

import sys
import tempfile
import urllib.error
import urllib.request
import json
import time
from pathlib import Path

import numpy as np

try:
    from openmm import openmm
    from openmm import unit
    from openmm import app
    from pdbfixer import PDBFixer
    print("✓ OpenMM and PDBFixer loaded successfully")
except ImportError as exc:
    print(f"Error importing required packages: {exc}")
    print("Install openmm, pdbfixer, and openmmforcefields.")
    sys.exit(1)

# Unit conversions (matching water-system conventions)
BOHR_TO_NM = 0.0529177           # 1 Bohr = 0.0529177 nm
HARTREE_TO_KJMOL = 2625.5        # 1 Hartree = 2625.5 kJ/mol
HARTREE_TO_CM = 219474.63        # 1 Hartree = 219474.63 cm^-1
AMU_TO_AU = 1822.888             # 1 amu = 1822.888 electron masses


def wavenumber_to_hartree(omega_cm):
    """Convert wavenumber (cm^-1) to Hartree (atomic units)."""
    return omega_cm / HARTREE_TO_CM


def download_pdb_file(pdb_id, timeout_s=30, retries=3, verbose=False):
    """
    Download a PDB file from RCSB and return the local path.

    Parameters
    ----------
    pdb_id : str
        PDB ID to download.
    timeout_s : float
        Timeout for HTTP request in seconds.
    retries : int
        Number of retries on failure.
    """
    pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    tmp_dir = Path(tempfile.mkdtemp(prefix=f"pdb_{pdb_id}_"))
    out_path = tmp_dir / f"{pdb_id}.pdb"

    for attempt in range(1, retries + 1):
        try:
            if verbose:
                print(f"  Downloading {url} (attempt {attempt}/{retries})...", flush=True)
            with urllib.request.urlopen(url, timeout=timeout_s) as response:
                out_path.write_bytes(response.read())
            if verbose:
                print(f"  Saved PDB to {out_path}", flush=True)
            return out_path
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt == retries:
                raise RuntimeError(
                    f"Failed to download PDB {pdb_id} after {retries} attempts. "
                    f"Last error: {exc}"
                ) from exc
            if verbose:
                print(f"  Download failed: {exc}. Retrying...", flush=True)
            time.sleep(2.0)

    raise RuntimeError(f"Unexpected failure downloading PDB {pdb_id}.")


def load_pdb_with_fixer(
    pdb_id="3UTL",
    pdb_path=None,
    fix_missing=True,
    add_hydrogens=True,
    ph=7.0,
    remove_heterogens=True,
    remove_water=True,
    verbose=False,
    download_timeout_s=30,
    download_retries=3,
):
    """
    Download or load a PDB structure and prepare it with PDBFixer.

    Parameters
    ----------
    pdb_id : str
        RCSB PDB ID to download (e.g., "3UTL").
    pdb_path : pathlib.Path or str
        Optional local PDB file path. If provided, download is skipped.
    fix_missing : bool
        Whether to add missing residues/atoms.
    add_hydrogens : bool
        Whether to add hydrogens at a specified pH.
    ph : float
        pH for protonation.
    remove_heterogens : bool
        Whether to remove heterogens (ligands, ions). Water removal is separate.
    remove_water : bool
        Whether to remove crystallographic water.
    """
    if pdb_path is None and pdb_id is None:
        raise ValueError("Provide either pdb_id or pdb_path.")

    if pdb_path is not None:
        fixer = PDBFixer(filename=str(pdb_path))
        pdb_tag = Path(pdb_path).stem
    else:
        pdb_tag = pdb_id
        if verbose:
            print("  Downloading PDB from RCSB...", flush=True)
        pdb_path = download_pdb_file(
            pdb_id=pdb_id,
            timeout_s=download_timeout_s,
            retries=download_retries,
            verbose=verbose,
        )
        fixer = PDBFixer(filename=str(pdb_path))

    print(f"\n--- PDBFixer: Loading {pdb_tag} ---")

    if remove_heterogens:
        fixer.removeHeterogens(keepWater=not remove_water)
        print("  Removed heterogens" + ("" if remove_water else " (kept waters)"))
    elif remove_water:
        fixer.removeHeterogens(keepWater=False)
        print("  Removed crystallographic water")

    if fix_missing:
        if verbose:
            print("  Finding missing residues/atoms...", flush=True)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        if verbose:
            print("  Adding missing atoms...", flush=True)
        fixer.addMissingAtoms()
        print("  Added missing residues/atoms")

    if add_hydrogens:
        if verbose:
            print("  Adding hydrogens...", flush=True)
        fixer.addMissingHydrogens(ph)
        print(f"  Added hydrogens at pH {ph:.1f}")

    return fixer


def prepare_protein_system(
    pdb_id="3UTL",
    pdb_path=None,
    ph=7.0,
    padding_nm=1.0,
    ionic_strength_molar=0.15,
    add_solvent=True,
    fix_missing=True,
    add_hydrogens=True,
    remove_heterogens=True,
    remove_water=True,
    forcefield_files=None,
    water_model="tip4pew",
    verbose=False,
    download_timeout_s=30,
    download_retries=3,
):
    """
    Build a protein + explicit solvent system using PDBFixer and OpenMM.

    Returns
    -------
    system : openmm.System
    topology : openmm.app.Topology
    positions : list of Vec3
    charges : list of float
        Charges for non-virtual particles (e).
    real_indices : np.ndarray
        Indices of non-virtual particles.
    """
    if forcefield_files is None:
        forcefield_files = ["amber14/protein.ff14SB.xml", "amber14/tip4pew.xml"]

    fixer = load_pdb_with_fixer(
        pdb_id=pdb_id,
        pdb_path=pdb_path,
        fix_missing=fix_missing,
        add_hydrogens=add_hydrogens,
        ph=ph,
        remove_heterogens=remove_heterogens,
        remove_water=remove_water,
        verbose=verbose,
        download_timeout_s=download_timeout_s,
        download_retries=download_retries,
    )

    modeller = app.Modeller(fixer.topology, fixer.positions)

    print(f"\n--- Building Force Field ---")
    print(f"  Force field: {', '.join(forcefield_files)}")
    if verbose:
        print("  Loading XML force fields...", flush=True)
    forcefield = app.ForceField(*forcefield_files)

    if add_solvent:
        print("\n--- Adding Explicit Solvent ---")
        if verbose:
            start_solvent = time.time()
            print("  Solvating system (this may take a while)...", flush=True)
        modeller.addSolvent(
            forcefield,
            model=water_model,
            padding=padding_nm * unit.nanometer,
            ionicStrength=ionic_strength_molar * unit.molar,
            neutralize=True,
        )
        if verbose:
            elapsed = time.time() - start_solvent
            print(f"  Solvation complete in {elapsed:.1f} s", flush=True)
        print(f"  Added solvent: {water_model}")
        print(f"  Padding: {padding_nm:.2f} nm")
        print(f"  Ionic strength: {ionic_strength_molar:.3f} M")

    print("\n--- Creating System ---")
    if verbose:
        start_system = time.time()
        print("  Building OpenMM System (PME, flexible water)...", flush=True)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=0.9 * unit.nanometer,
        constraints=None,
        rigidWater=False,
        ewaldErrorTolerance=0.0005,
    )
    if verbose:
        elapsed = time.time() - start_system
        print(f"  System creation complete in {elapsed:.1f} s", flush=True)

    # Identify non-virtual (real) particles for dipole output
    real_indices = np.array(
        [i for i in range(system.getNumParticles()) if not system.isVirtualSite(i)],
        dtype=int,
    )

    charges = []
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            for idx in real_indices:
                charge, sigma, epsilon = force.getParticleParameters(idx)
                charges.append(charge.value_in_unit(unit.elementary_charge))
            break

    positions = modeller.getPositions()

    print(f"  Total particles: {system.getNumParticles()}")
    print(f"  Real particles:  {len(real_indices)} (for dipole)")
    print(f"  Topology atoms:  {modeller.topology.getNumAtoms()}")

    return system, modeller.topology, positions, charges, real_indices


def add_cavity_particle(system, positions, omegac_au, photon_mass_amu):
    """Add cavity photon mode as a dummy particle."""
    print("\n--- Adding Cavity Particle ---")
    cavity_index = system.addParticle(photon_mass_amu * unit.amu)
    positions.append(openmm.Vec3(0, 0, 0) * unit.nanometer)

    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.addParticle(0.0, 0.1 * unit.nanometer, 0.0)
            break

    omegac_cm = omegac_au * HARTREE_TO_CM
    print(f"  Cavity particle index: {cavity_index}")
    print(f"  Omega_c: {omegac_au:.6f} a.u. = {omegac_cm:.1f} cm^-1")
    print(f"  Photon mass: {photon_mass_amu:.6f} amu")
    return cavity_index


def setup_cavity_coupling(system, cavity_index, omegac_au, lambda_coupling, photon_mass_amu):
    """Set up cavity coupling (force + finite-q displacer)."""
    print("\n--- Setting Up Cavity Coupling ---")
    print(f"  Lambda: {lambda_coupling}")
    print("  Coupling will be ON from t=0")

    cavity_force = openmm.CavityForce(cavity_index, omegac_au, lambda_coupling, photon_mass_amu)
    system.addForce(cavity_force)
    print("  CavityForce added")

    displacer = openmm.CavityParticleDisplacer(cavity_index, omegac_au, photon_mass_amu)
    displacer.setSwitchOnStep(0)
    displacer.setSwitchOnLambda(lambda_coupling)
    system.addForce(displacer)
    print("  CavityParticleDisplacer added")

    return cavity_force, displacer


def setup_thermostats(system, temperature_K, real_indices, tau_bussi_ps=1.0):
    """Apply Bussi thermostat to real particles (exclude cavity)."""
    print("\n--- Setting Up Thermostats ---")
    try:
        bussi = openmm.BussiThermostat(temperature_K, tau_bussi_ps)
        bussi.setApplyToAllParticles(False)
        for idx in real_indices:
            bussi.addParticle(int(idx))
        system.addForce(bussi)
        print(f"  BussiThermostat added for {bussi.getNumParticles()} particles")
        print(f"  Temperature: {temperature_K} K, Tau: {tau_bussi_ps} ps")
        return True
    except AttributeError:
        print("  BussiThermostat not available, using Langevin only")
        return False


def compute_dipole(state, charges, real_indices):
    """Compute total molecular dipole moment (excluding virtual sites)."""
    positions = state.getPositions(asNumpy=True)
    dipole = np.zeros(3)
    for charge, idx in zip(charges, real_indices):
        pos = positions[idx].value_in_unit(unit.nanometer)
        dipole += charge * np.array(pos)
    return dipole


def run_simulation(
    lambda_coupling=0.01,
    cavity_freq_cm=3663.0,
    pdb_id="3UTL",
    pdb_path=None,
    temperature_K=300.0,
    dt_ps=0.0005,
    equilibration_ps=100.0,
    production_ps=900.0,
    output_prefix="protein_cavity",
    padding_nm=1.0,
    ionic_strength_molar=0.15,
    ph=7.0,
    disable_dipole_output=False,
    verbose=False,
    download_timeout_s=30,
    download_retries=3,
    enable_cavity=True,
    stream_to_disk=False,
    write_pdb_snapshot=True,
    dipole_interval=1,
):
    """Run the protein cavity MD simulation with streaming output."""
    print("=" * 80)
    print("Protein Cavity Molecular Dynamics Simulation")
    print("=" * 80)

    omegac_au = wavenumber_to_hartree(cavity_freq_cm)
    photon_mass_amu = 1.0 / AMU_TO_AU

    # Compute effective dipole sampling rate for diagnostics
    sample_dt_ps = dt_ps * dipole_interval
    sample_freq_ps = 1.0 / sample_dt_ps  # ps^-1
    nyquist_cm = sample_freq_ps * 1e12 / (2.0 * 3e10)  # cm^-1

    print("\nSimulation Parameters:")
    print(f"  PDB: {pdb_id if pdb_path is None else Path(pdb_path).name}")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Cavity frequency: {cavity_freq_cm:.0f} cm^-1 (OH stretch)")
    print(f"  Lambda coupling: {lambda_coupling}")
    print(f"  Timestep: {dt_ps} ps")
    print(f"  Equilibration: {equilibration_ps} ps")
    print(f"  Production: {production_ps} ps")
    print(f"  Solvent padding: {padding_nm:.2f} nm")
    print(f"  Ionic strength: {ionic_strength_molar:.3f} M")
    print(f"  Dipole interval: every {dipole_interval} steps ({sample_dt_ps:.4f} ps)")
    print(f"  Nyquist frequency: {nyquist_cm:.0f} cm^-1")

    equilibration_steps = int(equilibration_ps / dt_ps)
    production_steps = int(production_ps / dt_ps)
    total_steps = equilibration_steps + production_steps

    system, topology, positions, charges, real_indices = prepare_protein_system(
        pdb_id=pdb_id,
        pdb_path=pdb_path,
        ph=ph,
        padding_nm=padding_nm,
        ionic_strength_molar=ionic_strength_molar,
        add_solvent=True,
        fix_missing=True,
        add_hydrogens=True,
        remove_heterogens=True,
        remove_water=True,
        forcefield_files=None,
        water_model="tip4pew",
        verbose=verbose,
        download_timeout_s=download_timeout_s,
        download_retries=download_retries,
    )

    cavity_index = None
    if enable_cavity:
        cavity_index = add_cavity_particle(system, positions, omegac_au, photon_mass_amu)
        cavity_force, displacer = setup_cavity_coupling(
            system, cavity_index, omegac_au, lambda_coupling, photon_mass_amu
        )
    else:
        print("\n--- Cavity Disabled ---")
        print("  No cavity particle or coupling forces will be added")
    setup_thermostats(system, temperature_K, real_indices=real_indices, tau_bussi_ps=1.0)

    print("\n--- Creating Integrator ---")
    # Use VerletIntegrator when BussiThermostat is active to avoid
    # conflicting random number seeds between two stochastic thermostats.
    # Bussi handles thermostatting; Verlet handles time integration.
    has_bussi = any(
        isinstance(f, openmm.BussiThermostat) for f in
        [system.getForce(i) for i in range(system.getNumForces())]
    )
    if has_bussi:
        integrator = openmm.VerletIntegrator(dt_ps * unit.picosecond)
        print(f"  VerletIntegrator: dt = {dt_ps} ps (thermostat: Bussi)")
    else:
        friction = 0.5  # ps^-1
        integrator = openmm.LangevinMiddleIntegrator(
            temperature_K * unit.kelvin,
            friction / unit.picosecond,
            dt_ps * unit.picosecond,
        )
        print(f"  LangevinMiddleIntegrator: friction = {friction} ps^-1, dt = {dt_ps} ps")

    print("\n--- Creating Simulation Context ---")
    platform = openmm.Platform.getPlatformByName("CUDA")
    properties = {"Precision": "mixed"}
    print("  Using CUDA platform (GPU acceleration)")

    context = openmm.Context(system, integrator, platform, properties)
    context.setPositions(positions)
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)

    print("\n--- Minimizing Energy ---")
    openmm.LocalEnergyMinimizer.minimize(context, maxIterations=1000)
    state = context.getState(getEnergy=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    print(f"  Initial PE: {pe:.1f} kJ/mol")

    dipole_times = []
    dipole_values = []
    cavity_values = []
    pdb_written = False

    output_file = f"{output_prefix}_lambda{lambda_coupling:.4f}.npz"
    pdb_file = f"{output_prefix}_lambda{lambda_coupling:.4f}.pdb"
    metadata_file = f"{output_prefix}_lambda{lambda_coupling:.4f}_metadata.json"

    # Number of dipole data points in production
    n_dipole_points = production_steps // dipole_interval
    report_interval = 1000  # steps between progress reports

    if stream_to_disk and not disable_dipole_output:
        stream_path = f"{output_prefix}_lambda{lambda_coupling:.4f}_stream.npy"
        stream_dtype = np.dtype(
            [
                ("time_ps", np.float64),
                ("dipole_nm", np.float64, (3,)),
                ("cavity_nm", np.float64, (3,)),
            ]
        )
        stream_mm = np.lib.format.open_memmap(
            stream_path,
            mode="w+",
            dtype=stream_dtype,
            shape=(n_dipole_points,),
        )
        prod_index = 0
        print("  Streaming dipole to disk (single memmap):")
        print(f"    stream: {stream_path}")
        print(f"    max data points: {n_dipole_points}")

    print("\n--- Running Simulation ---")
    print(f"  Total steps: {total_steps}")
    print(f"  Equilibration: {equilibration_steps} steps ({equilibration_ps} ps)")
    print(f"  Production: {production_steps} steps ({production_ps} ps)")
    print(f"  Output file: {output_file}")
    print(f"  Progress reports every {report_interval} steps")

    start_time = time.time()
    last_report_time = start_time
    last_report_step = 0

    # --- Equilibration phase: run in large batches ---
    print(f"\n  Phase 1: Equilibration ({equilibration_steps} steps)...")
    equil_done = 0
    while equil_done < equilibration_steps:
        batch = min(report_interval, equilibration_steps - equil_done)
        integrator.step(batch)
        equil_done += batch

        step = equil_done - 1  # 0-indexed current step
        pct = 100 * equil_done / total_steps
        current_time = time.time()
        dt_report = current_time - last_report_time
        steps_since_last = equil_done - last_report_step
        inst_rate = steps_since_last / max(dt_report, 1e-8)
        ns_per_day = (inst_rate * dt_ps * 86400) / 1000
        eta = (total_steps - equil_done) / max(inst_rate, 1e-8)
        last_report_time = current_time
        last_report_step = equil_done

        state = context.getState(getEnergy=True, getPositions=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        dipole_now = compute_dipole(state, charges, real_indices)
        sim_time_ps = equil_done * dt_ps
        print(
            f"[{pct:5.1f}%] EQUIL | t={sim_time_ps:7.1f} ps | PE: {pe:10.1f} kJ/mol | "
            f"d_xy: ({dipole_now[0]:6.2f}, {dipole_now[1]:6.2f}) e·nm | "
            f"Speed: {ns_per_day:5.2f} ns/day | ETA: {eta/60:5.1f}m",
            flush=True,
        )

    # --- Production phase: step in dipole_interval batches, record dipole ---
    print(f"\n  Phase 2: Production ({production_steps} steps, recording every {dipole_interval} steps)...")
    prod_done = 0
    last_report_prod = 0
    while prod_done < production_steps:
        # Take dipole_interval steps at once
        integrator.step(dipole_interval)
        prod_done += dipole_interval

        # Record dipole
        if not disable_dipole_output:
            state = context.getState(getPositions=True)
            pos = state.getPositions(asNumpy=True)
            dipole = compute_dipole(state, charges, real_indices)
            if cavity_index is not None:
                q_cav = pos[cavity_index].value_in_unit(unit.nanometer)
            else:
                q_cav = np.zeros(3)

            if stream_to_disk:
                stream_mm[prod_index]["time_ps"] = prod_done * dt_ps
                stream_mm[prod_index]["dipole_nm"] = dipole
                stream_mm[prod_index]["cavity_nm"] = q_cav
                prod_index += 1
            else:
                dipole_times.append(prod_done * dt_ps)
                dipole_values.append(dipole.copy())
                cavity_values.append(q_cav.copy())

        # Progress report
        if prod_done - last_report_prod >= report_interval or prod_done >= production_steps:
            last_report_prod = prod_done
            total_done = equilibration_steps + prod_done
            pct = 100 * total_done / total_steps
            current_time = time.time()
            dt_report = current_time - last_report_time
            steps_since_last = total_done - last_report_step
            inst_rate = steps_since_last / max(dt_report, 1e-8)
            ns_per_day = (inst_rate * dt_ps * 86400) / 1000
            eta = (total_steps - total_done) / max(inst_rate, 1e-8)
            last_report_time = current_time
            last_report_step = total_done

            state = context.getState(getEnergy=True, getPositions=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            pos = state.getPositions(asNumpy=True)
            if cavity_index is not None:
                q_cav_report = pos[cavity_index].value_in_unit(unit.nanometer)
            else:
                q_cav_report = np.zeros(3)
            dipole_now = compute_dipole(state, charges, real_indices)

            if not disable_dipole_output:
                metadata = {
                    "lambda_coupling": lambda_coupling,
                    "cavity_freq_cm": cavity_freq_cm,
                    "omega_c_au": omegac_au,
                    "dt_ps": dt_ps,
                    "dipole_dt_ps": sample_dt_ps,
                    "temperature_K": temperature_K,
                    "photon_mass_amu": photon_mass_amu,
                    "pdb_id": pdb_id,
                    "padding_nm": padding_nm,
                    "ionic_strength_molar": ionic_strength_molar,
                    "forcefield": "amber14/protein.ff14SB.xml + amber14/tip4pew.xml",
                    "water_model": "tip4pew",
                    "cavity_enabled": enable_cavity,
                    "dipole_interval": dipole_interval,
                    "status": "running",
                    "step": total_done,
                    "total_steps": total_steps,
                }
                if stream_to_disk:
                    with open(metadata_file, "w") as handle:
                        json.dump(metadata, handle, indent=2)
                elif len(dipole_times) > 0:
                    np.savez(
                        output_file,
                        time_ps=np.array(dipole_times),
                        dipole_nm=np.array(dipole_values),
                        cavity_nm=np.array(cavity_values),
                        metadata=metadata,
                    )
                if write_pdb_snapshot and not pdb_written:
                    # Trim positions to topology atoms (exclude cavity particle)
                    n_topo = topology.getNumAtoms()
                    pdb_pos = pos[:n_topo]
                    with open(pdb_file, "w") as handle:
                        app.PDBFile.writeFile(topology, pdb_pos, handle)
                    pdb_written = True
                    print(f"         Wrote PDB snapshot: {pdb_file}")

            sim_time_ps = total_done * dt_ps
            print(
                f"[{pct:5.1f}%] PROD | t={sim_time_ps:7.1f} ps | PE: {pe:10.1f} kJ/mol | "
                f"d_xy: ({dipole_now[0]:6.2f}, {dipole_now[1]:6.2f}) e·nm | "
                f"q_cav: ({q_cav_report[0]:6.3f}, {q_cav_report[1]:6.3f}) nm | "
                f"Speed: {ns_per_day:5.2f} ns/day | ETA: {eta/60:5.1f}m",
                flush=True,
            )

    elapsed = time.time() - start_time
    print("\n--- Simulation Complete ---")
    print(f"  Total time: {elapsed/60:.1f} min ({elapsed/3600:.2f} hr)")
    print(f"  Performance: {total_steps / elapsed:.1f} steps/s")

    if not disable_dipole_output:
        metadata = {
            "lambda_coupling": lambda_coupling,
            "cavity_freq_cm": cavity_freq_cm,
            "omega_c_au": omegac_au,
            "dt_ps": dt_ps,
            "dipole_dt_ps": sample_dt_ps,
            "temperature_K": temperature_K,
            "photon_mass_amu": photon_mass_amu,
            "pdb_id": pdb_id,
            "padding_nm": padding_nm,
            "ionic_strength_molar": ionic_strength_molar,
            "forcefield": "amber14/protein.ff14SB.xml + amber14/tip4pew.xml",
            "water_model": "tip4pew",
            "cavity_enabled": enable_cavity,
            "dipole_interval": dipole_interval,
            "status": "complete",
            "step": total_steps,
            "total_steps": total_steps,
            "elapsed_time_s": elapsed,
        }
        if stream_to_disk:
            time_ps = np.array(stream_mm[:prod_index]["time_ps"])
            dipole_nm = np.array(stream_mm[:prod_index]["dipole_nm"])
            cavity_nm = np.array(stream_mm[:prod_index]["cavity_nm"])
            np.savez(
                output_file,
                time_ps=time_ps,
                dipole_nm=dipole_nm,
                cavity_nm=cavity_nm,
                metadata=metadata,
            )
            with open(metadata_file, "w") as handle:
                json.dump(metadata, handle, indent=2)
            print(f"  Saved: {output_file}")
            print(f"  Data points: {prod_index}")
        else:
            np.savez(
                output_file,
                time_ps=np.array(dipole_times),
                dipole_nm=np.array(dipole_values),
                cavity_nm=np.array(cavity_values),
                metadata=metadata,
            )
            print(f"  Saved: {output_file}")
            print(f"  Data points: {len(dipole_times)}")
        if write_pdb_snapshot and not pdb_written:
            state = context.getState(getPositions=True)
            pos = state.getPositions(asNumpy=True)
            n_topo = topology.getNumAtoms()
            pdb_pos = pos[:n_topo]
            with open(pdb_file, "w") as handle:
                app.PDBFile.writeFile(topology, pdb_pos, handle)
            print(f"  Wrote PDB snapshot: {pdb_file}")
        if stream_to_disk and prod_index > 0:
            print(f"  Production time: {stream_mm[prod_index - 1]['time_ps']:.1f} ps")
        elif dipole_times:
            print(f"  Production time: {dipole_times[-1]:.1f} ps")
    else:
        print("  Dipole output: DISABLED (benchmark mode)")


def main():
    """Main entry point with command-line argument parsing."""
    import argparse

    parser = argparse.ArgumentParser(description="Protein Cavity MD Simulation")
    parser.add_argument("--lambda", type=float, default=0.01, dest="lambda_coupling",
                        help="Coupling strength (default: 0.01)")
    parser.add_argument("--cavity-freq", type=float, default=3663.0, dest="cavity_freq",
                        help="Cavity frequency in cm^-1 (default: 3663 for OH stretch)")
    parser.add_argument("--pdb-id", type=str, default="3UTL",
                        help="RCSB PDB ID to download (default: 3UTL)")
    parser.add_argument("--pdb-path", type=str, default=None,
                        help="Optional local PDB path (overrides --pdb-id)")
    parser.add_argument("--temp", type=float, default=300.0,
                        help="Temperature in K (default: 300)")
    parser.add_argument("--dt", type=float, default=0.0005,
                        help="Timestep in ps (default: 0.0005 = 0.5 fs)")
    parser.add_argument("--equil", type=float, default=100.0,
                        help="Equilibration time in ps (default: 100)")
    parser.add_argument("--prod", type=float, default=900.0,
                        help="Production time in ps (default: 900)")
    parser.add_argument("--output", type=str, default="protein_cavity",
                        help="Output prefix (default: protein_cavity)")
    parser.add_argument("--padding", type=float, default=1.0,
                        help="Solvent padding in nm (default: 1.0)")
    parser.add_argument("--ionic-strength", type=float, default=0.15,
                        help="Ionic strength in M (default: 0.15)")
    parser.add_argument("--ph", type=float, default=7.0,
                        help="pH for protonation (default: 7.0)")
    parser.add_argument("--test", action="store_true",
                        help="Quick test run (reduced runtime)")
    parser.add_argument("--no-dipole", action="store_true", dest="no_dipole",
                        help="Disable dipole output for speed benchmark")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose setup timing output")
    parser.add_argument("--download-timeout", type=float, default=30.0,
                        help="PDB download timeout in seconds (default: 30)")
    parser.add_argument("--download-retries", type=int, default=3,
                        help="PDB download retries (default: 3)")
    parser.add_argument("--no-cavity", action="store_true",
                        help="Disable cavity particle and coupling forces")
    parser.add_argument("--stream-to-disk", action="store_true",
                        help="Stream dipole/cavity to disk instead of growing arrays")
    parser.add_argument("--no-pdb", action="store_true",
                        help="Disable PDB snapshot output")
    parser.add_argument("--dipole-interval", type=int, default=1, dest="dipole_interval",
                        help="Record dipole every N steps (default: 1). "
                             "Higher values speed up production. "
                             "E.g., 4 with dt=0.5fs gives 2fs sampling (Nyquist ~8333 cm^-1).")

    args = parser.parse_args()

    if args.test:
        print("\n*** QUICK TEST MODE ***\n")
        args.equil = 1.0
        args.prod = 1.0

    run_simulation(
        lambda_coupling=args.lambda_coupling,
        cavity_freq_cm=args.cavity_freq,
        pdb_id=args.pdb_id,
        pdb_path=args.pdb_path,
        temperature_K=args.temp,
        dt_ps=args.dt,
        equilibration_ps=args.equil,
        production_ps=args.prod,
        output_prefix=args.output,
        padding_nm=args.padding,
        ionic_strength_molar=args.ionic_strength,
        ph=args.ph,
        disable_dipole_output=args.no_dipole,
        verbose=args.verbose,
        download_timeout_s=args.download_timeout,
        download_retries=args.download_retries,
        enable_cavity=not args.no_cavity,
        stream_to_disk=args.stream_to_disk,
        write_pdb_snapshot=not args.no_pdb,
        dipole_interval=args.dipole_interval,
    )


if __name__ == "__main__":
    main()

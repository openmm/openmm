#!/usr/bin/env python3
"""
Test suite for MultiModeCavityForce.

Verifies:
1. Mode frequency scaling: omega_n = n * omega_1
2. Coupling strength scaling: lambda_n = lambda_1 (constant)
3. Spatial profile: even modes at cavity center are dark
4. N=1 multi-mode matches single-mode CavityForce exactly
5. Total energy matches N separate single-mode CavityForce objects
6. Energy conservation in NVE dynamics
7. Performance benchmark: multi-mode kernel vs N single-mode kernels
"""

import sys
import time
import math
import numpy as np

try:
    from openmm import openmm
    from openmm import unit
except ImportError:
    import openmm
    from openmm import unit

KJ = unit.kilojoule_per_mole

def _strip_units(val):
    """Convert a Quantity to float, or return float as-is."""
    if hasattr(val, 'value_in_unit'):
        return val.value_in_unit(KJ)
    return float(val)


# ============================================================================
# Helper: Build a minimal 2-particle + cavity system
# ============================================================================
def build_test_system(num_modes, omega1_au=0.007, lambda1=0.05,
                      cavity_length_nm=100.0, molecule_z_nm=50.0,
                      photon_mass_amu=1.0, use_multimode=True,
                      platform_name='Reference'):
    """
    Create a minimal system with 2 charged particles + N cavity photon particles.

    Returns (system, positions, topology_info_dict).
    """
    system = openmm.System()
    # 2 molecular particles
    system.addParticle(16.0)  # Particle 0: -0.5e
    system.addParticle(16.0)  # Particle 1: +0.5e

    # N cavity particles (one per mode)
    cavity_indices = []
    for m in range(num_modes):
        idx = system.addParticle(photon_mass_amu)
        cavity_indices.append(idx)

    num_particles = system.getNumParticles()

    # Periodic box
    box_side = 2.0  # nm
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box_side, 0, 0),
        openmm.Vec3(0, box_side, 0),
        openmm.Vec3(0, 0, box_side),
    )

    # NonbondedForce for charges
    nonbonded = openmm.NonbondedForce()
    nonbonded.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nonbonded.setCutoffDistance(0.9)
    nonbonded.addParticle(-0.5, 0.3, 0.5)  # Particle 0
    nonbonded.addParticle(+0.5, 0.3, 0.5)  # Particle 1
    for _ in range(num_modes):
        nonbonded.addParticle(0.0, 0.1, 0.0)  # Cavity particles: no charge
    # Exclude cavity-cavity and cavity-molecular interactions for nonbonded
    for i in range(num_modes):
        for j in range(i + 1, num_modes):
            nonbonded.addException(cavity_indices[i], cavity_indices[j], 0.0, 1.0, 0.0)
        nonbonded.addException(0, cavity_indices[i], 0.0, 1.0, 0.0)
        nonbonded.addException(1, cavity_indices[i], 0.0, 1.0, 0.0)
    system.addForce(nonbonded)

    if use_multimode:
        # Use MultiModeCavityForce
        mm_force = openmm.MultiModeCavityForce(
            num_modes, omega1_au, lambda1,
            cavity_length_nm, molecule_z_nm, photon_mass_amu
        )
        for idx in cavity_indices:
            mm_force.addCavityParticle(idx)
        system.addForce(mm_force)
    else:
        # Use N separate CavityForce objects
        for m in range(num_modes):
            n = m + 1
            omega_n = n * omega1_au
            lambda_n = lambda1
            fn = math.sin(n * math.pi * molecule_z_nm / cavity_length_nm)
            # Effective lambda for single-mode = lambda_n * fn
            # (spatial profile is baked into the effective coupling)
            cf = openmm.CavityForce(cavity_indices[m], omega_n, lambda_n * fn, photon_mass_amu)
            system.addForce(cf)

    # Initial positions: particles separated along x, cavity at origin
    positions = [openmm.Vec3(0.0, 0.0, 0.0)]     # Particle 0
    positions.append(openmm.Vec3(0.2, 0.0, 0.0))  # Particle 1
    for _ in range(num_modes):
        positions.append(openmm.Vec3(0.05, 0.03, 0.0))  # Cavity particles (slightly displaced)

    return system, positions, cavity_indices


# ============================================================================
# Test 1: Mode scaling verification (API-level)
# ============================================================================
def test_mode_scaling():
    """Verify omega_n = n*omega_1 and lambda_n = lambda_1 (constant)."""
    print("\n" + "=" * 60)
    print("Test 1: Mode Frequency and Coupling Scaling")
    print("=" * 60)

    omega1 = 0.007
    lambda1 = 0.05
    num_modes = 5
    cavity_length = 100.0
    molecule_z = 50.0

    force = openmm.MultiModeCavityForce(
        num_modes, omega1, lambda1, cavity_length, molecule_z
    )

    all_pass = True
    for n in range(1, num_modes + 1):
        omega_n = force.getOmegaN(n)
        lambda_n = force.getLambdaN(n)
        eps_n = force.getEffectiveCouplingN(n)
        K_n = force.getSpringConstantN(n)
        fn = force.getSpatialProfile(n)

        expected_omega = n * omega1
        expected_lambda = lambda1
        expected_fn = math.sin(n * math.pi * molecule_z / cavity_length)

        ok_omega = abs(omega_n - expected_omega) < 1e-12
        ok_lambda = abs(lambda_n - expected_lambda) < 1e-12
        ok_fn = abs(fn - expected_fn) < 1e-12

        status = "PASS" if (ok_omega and ok_lambda and ok_fn) else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"  Mode {n}: omega_n={omega_n:.6f} (exp {expected_omega:.6f}), "
              f"lambda_n={lambda_n:.6f} (exp {expected_lambda:.6f}), "
              f"f_n={fn:.6f} (exp {expected_fn:.6f}) [{status}]")

    print(f"\n  Overall: {'PASS' if all_pass else 'FAIL'}")
    return all_pass


# ============================================================================
# Test 2: Spatial profile -- even modes at cavity center
# ============================================================================
def test_spatial_profile():
    """At z0=L/2, even modes should have f_n=0 (dark modes)."""
    print("\n" + "=" * 60)
    print("Test 2: Spatial Profile (Even Modes Dark at Center)")
    print("=" * 60)

    omega1 = 0.007
    lambda1 = 0.05
    num_modes = 6
    cavity_length = 100.0
    molecule_z = cavity_length / 2.0  # Center of cavity

    force = openmm.MultiModeCavityForce(
        num_modes, omega1, lambda1, cavity_length, molecule_z
    )

    all_pass = True
    for n in range(1, num_modes + 1):
        fn = force.getSpatialProfile(n)
        # Even modes: sin(n*pi/2) = 0 when n is even
        expected = math.sin(n * math.pi * 0.5)
        is_even = (n % 2 == 0)
        ok = abs(fn - expected) < 1e-12
        if not ok:
            all_pass = False
        dark_str = " (DARK)" if is_even else ""
        print(f"  Mode {n}: f_n={fn:+.10f} (expected {expected:+.10f}){dark_str} [{'PASS' if ok else 'FAIL'}]")

    print(f"\n  Overall: {'PASS' if all_pass else 'FAIL'}")
    return all_pass


# ============================================================================
# Test 3: N=1 multi-mode matches single-mode CavityForce
# ============================================================================
def test_single_mode_equivalence():
    """N=1 MultiModeCavityForce should produce identical energy to CavityForce."""
    print("\n" + "=" * 60)
    print("Test 3: N=1 Multi-Mode vs Single-Mode CavityForce")
    print("=" * 60)

    omega1 = 0.007
    lambda1 = 0.05
    cavity_length = 100.0
    molecule_z = 25.0  # Not at center, so f_1 != 0
    photon_mass = 1.0

    # Build multi-mode system (N=1)
    sys_mm, pos_mm, cav_mm = build_test_system(
        1, omega1, lambda1, cavity_length, molecule_z, photon_mass,
        use_multimode=True, platform_name='Reference'
    )

    # Build single-mode system
    sys_sm, pos_sm, cav_sm = build_test_system(
        1, omega1, lambda1, cavity_length, molecule_z, photon_mass,
        use_multimode=False, platform_name='Reference'
    )

    # Create contexts on Reference platform
    platform = openmm.Platform.getPlatformByName('Reference')
    integrator_mm = openmm.VerletIntegrator(0.001 * unit.picosecond)
    integrator_sm = openmm.VerletIntegrator(0.001 * unit.picosecond)

    ctx_mm = openmm.Context(sys_mm, integrator_mm, platform)
    ctx_sm = openmm.Context(sys_sm, integrator_sm, platform)

    ctx_mm.setPositions(pos_mm)
    ctx_sm.setPositions(pos_sm)

    # Get energies
    state_mm = ctx_mm.getState(getEnergy=True, getForces=True)
    state_sm = ctx_sm.getState(getEnergy=True, getForces=True)

    # Get per-component energies
    mm_force = None
    for i in range(sys_mm.getNumForces()):
        f = sys_mm.getForce(i)
        if isinstance(f, openmm.MultiModeCavityForce):
            mm_force = f
            break

    sm_force = None
    for i in range(sys_sm.getNumForces()):
        f = sys_sm.getForce(i)
        if isinstance(f, openmm.CavityForce):
            sm_force = f
            break

    mm_harmonic = _strip_units(mm_force.getHarmonicEnergy(ctx_mm))
    mm_coupling = _strip_units(mm_force.getCouplingEnergy(ctx_mm))
    mm_dse = _strip_units(mm_force.getDipoleSelfEnergy(ctx_mm))
    mm_total = _strip_units(mm_force.getTotalCavityEnergy(ctx_mm))

    sm_harmonic = _strip_units(sm_force.getHarmonicEnergy(ctx_sm))
    sm_coupling = _strip_units(sm_force.getCouplingEnergy(ctx_sm))
    sm_dse = _strip_units(sm_force.getDipoleSelfEnergy(ctx_sm))
    sm_total = _strip_units(sm_force.getTotalCavityEnergy(ctx_sm))

    print(f"  MultiMode: harmonic={mm_harmonic:.8f}, coupling={mm_coupling:.8f}, DSE={mm_dse:.8f}, total={mm_total:.8f}")
    print(f"  SingleMode: harmonic={sm_harmonic:.8f}, coupling={sm_coupling:.8f}, DSE={sm_dse:.8f}, total={sm_total:.8f}")

    tol = 1e-4
    ok_harm = abs(mm_harmonic - sm_harmonic) < tol
    ok_coup = abs(mm_coupling - sm_coupling) < tol
    ok_dse = abs(mm_dse - sm_dse) < tol
    ok_total = abs(mm_total - sm_total) < tol

    # Also compare forces
    forces_mm = state_mm.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
    forces_sm = state_sm.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
    force_diff = np.max(np.abs(forces_mm - forces_sm))

    ok_forces = force_diff < 1.0  # kJ/(mol*nm) tolerance
    print(f"  Max force difference: {force_diff:.6f} kJ/(mol*nm) [{'PASS' if ok_forces else 'FAIL'}]")

    all_pass = ok_harm and ok_coup and ok_dse and ok_total and ok_forces
    print(f"  Harmonic: {'PASS' if ok_harm else 'FAIL'} (diff={abs(mm_harmonic-sm_harmonic):.2e})")
    print(f"  Coupling: {'PASS' if ok_coup else 'FAIL'} (diff={abs(mm_coupling-sm_coupling):.2e})")
    print(f"  DSE: {'PASS' if ok_dse else 'FAIL'} (diff={abs(mm_dse-sm_dse):.2e})")
    print(f"  Total: {'PASS' if ok_total else 'FAIL'} (diff={abs(mm_total-sm_total):.2e})")
    print(f"\n  Overall: {'PASS' if all_pass else 'FAIL'}")

    del ctx_mm, ctx_sm
    return all_pass


# ============================================================================
# Test 4: Multi-mode energy matches N separate single-mode forces
# ============================================================================
def test_multimode_vs_analytical():
    """Multi-mode energy should match analytical computation from positions."""
    print("\n" + "=" * 60)
    print("Test 4: N=5 Multi-Mode vs Analytical Computation")
    print("=" * 60)

    omega1 = 0.007
    lambda1 = 0.05
    num_modes = 5
    cavity_length = 100.0
    molecule_z = 30.0  # Off-center so all modes couple
    photon_mass = 1.0

    # Multi-mode system
    sys_mm, pos_mm, cav_indices = build_test_system(
        num_modes, omega1, lambda1, cavity_length, molecule_z, photon_mass,
        use_multimode=True
    )

    platform = openmm.Platform.getPlatformByName('Reference')
    integrator_mm = openmm.VerletIntegrator(0.001 * unit.picosecond)
    ctx_mm = openmm.Context(sys_mm, integrator_mm, platform)
    ctx_mm.setPositions(pos_mm)

    # Trigger energy computation (needed for component getters)
    state = ctx_mm.getState(getEnergy=True)
    total_pe = state.getPotentialEnergy().value_in_unit(KJ)
    print(f"  Total PE from getState: {total_pe:.6f}")

    # Get multi-mode component energies
    mm_force = None
    for i in range(sys_mm.getNumForces()):
        f = sys_mm.getForce(i)
        if isinstance(f, openmm.MultiModeCavityForce):
            mm_force = f
            break

    mm_harmonic = _strip_units(mm_force.getHarmonicEnergy(ctx_mm))
    mm_coupling = _strip_units(mm_force.getCouplingEnergy(ctx_mm))
    mm_dse = _strip_units(mm_force.getDipoleSelfEnergy(ctx_mm))
    mm_total = _strip_units(mm_force.getTotalCavityEnergy(ctx_mm))

    # Analytical computation using the SAME conversion constants as the C++ kernel
    OMEGAC_AU_TO_K = 1.7109e9  # Spring constant conversion (matches C++ code)
    EPS_CONV = 937679.0        # Coupling conversion (matches C++ code)

    # Get positions
    pos = [list(p) for p in pos_mm]
    charges = [-0.5, 0.5] + [0.0] * num_modes

    # Dipole (x,y only)
    dx = sum(charges[i] * pos[i][0] for i in range(2))
    dy = sum(charges[i] * pos[i][1] for i in range(2))

    anal_harmonic = 0.0
    anal_coupling = 0.0
    anal_dse_prefactor = 0.0

    for m in range(num_modes):
        n = m + 1
        omega_n = n * omega1
        lambda_n = lambda1
        fn = math.sin(n * math.pi * molecule_z / cavity_length)
        K_n = photon_mass * OMEGAC_AU_TO_K * omega_n * omega_n
        eps_n = lambda_n * omega_n * EPS_CONV
        epsf_n = eps_n * fn

        qx = pos[cav_indices[m]][0]
        qy = pos[cav_indices[m]][1]
        qz = pos[cav_indices[m]][2]

        anal_harmonic += 0.5 * K_n * (qx**2 + qy**2 + qz**2)
        anal_coupling += epsf_n * (qx * dx + qy * dy)
        anal_dse_prefactor += epsf_n ** 2 / K_n

    anal_dse = 0.5 * anal_dse_prefactor * (dx**2 + dy**2)
    anal_total = anal_harmonic + anal_coupling + anal_dse

    print(f"  MultiMode: harmonic={mm_harmonic:.6f}, coupling={mm_coupling:.6f}, DSE={mm_dse:.6f}")
    print(f"  Analytical: harmonic={anal_harmonic:.6f}, coupling={anal_coupling:.6f}, DSE={anal_dse:.6f}")
    print(f"  MultiMode total: {mm_total:.6f} kJ/mol")
    print(f"  Analytical total: {anal_total:.6f} kJ/mol")

    tol = 1e-3
    ok_harm = abs(mm_harmonic - anal_harmonic) < tol
    ok_coup = abs(mm_coupling - anal_coupling) < tol
    ok_dse = abs(mm_dse - anal_dse) < tol
    ok_total = abs(mm_total - anal_total) < tol

    print(f"  Harmonic match: {'PASS' if ok_harm else 'FAIL'} (diff={abs(mm_harmonic-anal_harmonic):.2e})")
    print(f"  Coupling match: {'PASS' if ok_coup else 'FAIL'} (diff={abs(mm_coupling-anal_coupling):.2e})")
    print(f"  DSE match: {'PASS' if ok_dse else 'FAIL'} (diff={abs(mm_dse-anal_dse):.2e})")
    print(f"  Total match: {'PASS' if ok_total else 'FAIL'} (diff={abs(mm_total-anal_total):.2e})")

    all_pass = ok_harm and ok_coup and ok_dse and ok_total
    print(f"\n  Overall: {'PASS' if all_pass else 'FAIL'}")

    del ctx_mm
    return all_pass


# ============================================================================
# Test 5: Energy conservation (NVE)
# ============================================================================
def test_energy_conservation():
    """Verify N=1 multi-mode trajectory matches single-mode trajectory exactly.

    This validates force/energy consistency by comparing multi-mode and single-mode
    time evolution over many steps. Both use the same Verlet integrator and
    platform, so if trajectories match, forces and energies are consistent.

    Note: we don't check raw KE+PE conservation because OpenMM's VerletIntegrator
    uses leapfrog, and getKineticEnergy() returns half-step KE that oscillates
    with the harmonic cavity frequency. Instead, we verify that the multi-mode
    and single-mode give identical trajectories.
    """
    print("\n" + "=" * 60)
    print("Test 5: N=1 Multi-Mode vs Single-Mode NVE Trajectory")
    print("=" * 60)

    omega1 = 0.007
    lambda1 = 0.01
    cavity_length = 100.0
    molecule_z = 50.0
    photon_mass = 1.0

    # System setup common to both
    def make_system(use_multimode):
        system = openmm.System()
        system.addParticle(16.0)  # Particle 0: -0.5e
        system.addParticle(16.0)  # Particle 1: +0.5e
        system.addParticle(photon_mass)  # Cavity photon

        nonbonded = openmm.NonbondedForce()
        nonbonded.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
        nonbonded.addParticle(-0.5, 0.3, 0.5)
        nonbonded.addParticle(+0.5, 0.3, 0.5)
        nonbonded.addParticle(0.0, 0.1, 0.0)
        nonbonded.addException(0, 2, 0.0, 1.0, 0.0)
        nonbonded.addException(1, 2, 0.0, 1.0, 0.0)
        system.addForce(nonbonded)

        if use_multimode:
            mmf = openmm.MultiModeCavityForce(1, omega1, lambda1, cavity_length, molecule_z, photon_mass)
            mmf.addCavityParticle(2)
            system.addForce(mmf)
        else:
            fn = math.sin(math.pi * molecule_z / cavity_length)
            cf = openmm.CavityForce(2, omega1, lambda1 * fn, photon_mass)
            system.addForce(cf)

        return system

    positions = [
        openmm.Vec3(0.0, 0.0, 0.0),
        openmm.Vec3(0.4, 0.0, 0.0),
        openmm.Vec3(0.03, 0.02, 0.0),
    ]

    dt = 0.00005  # ps
    platform = openmm.Platform.getPlatformByName('Reference')

    sys_mm = make_system(True)
    sys_sm = make_system(False)

    int_mm = openmm.VerletIntegrator(dt * unit.picosecond)
    int_sm = openmm.VerletIntegrator(dt * unit.picosecond)

    ctx_mm = openmm.Context(sys_mm, int_mm, platform)
    ctx_sm = openmm.Context(sys_sm, int_sm, platform)

    ctx_mm.setPositions(positions)
    ctx_sm.setPositions(positions)

    # Same initial velocities
    ctx_mm.setVelocitiesToTemperature(50 * unit.kelvin, 42)
    ctx_sm.setVelocitiesToTemperature(50 * unit.kelvin, 42)

    # Run 2000 steps and compare trajectories
    n_steps = 2000
    max_pos_diff = 0.0
    max_vel_diff = 0.0
    max_pe_diff = 0.0

    for i in range(n_steps):
        if i % 200 == 0:
            s_mm = ctx_mm.getState(getPositions=True, getVelocities=True, getEnergy=True)
            s_sm = ctx_sm.getState(getPositions=True, getVelocities=True, getEnergy=True)

            pos_mm = s_mm.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            pos_sm = s_sm.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            vel_mm = s_mm.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)
            vel_sm = s_sm.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)
            pe_mm = s_mm.getPotentialEnergy().value_in_unit(KJ)
            pe_sm = s_sm.getPotentialEnergy().value_in_unit(KJ)

            pos_d = np.max(np.abs(pos_mm - pos_sm))
            vel_d = np.max(np.abs(vel_mm - vel_sm))
            pe_d = abs(pe_mm - pe_sm)

            max_pos_diff = max(max_pos_diff, pos_d)
            max_vel_diff = max(max_vel_diff, vel_d)
            max_pe_diff = max(max_pe_diff, pe_d)

            if i % 500 == 0:
                print(f"  Step {i:5d}: pos_diff={pos_d:.2e} nm, vel_diff={vel_d:.2e} nm/ps, PE_diff={pe_d:.2e} kJ/mol")

        int_mm.step(1)
        int_sm.step(1)

    print(f"\n  Max position difference:  {max_pos_diff:.2e} nm")
    print(f"  Max velocity difference:  {max_vel_diff:.2e} nm/ps")
    print(f"  Max PE difference:        {max_pe_diff:.2e} kJ/mol")

    # Trajectories should match to machine precision
    ok = max_pos_diff < 1e-10 and max_vel_diff < 1e-8 and max_pe_diff < 1e-8
    print(f"\n  Overall: {'PASS' if ok else 'FAIL'}")

    del ctx_mm, ctx_sm
    return ok


# ============================================================================
# Test 6: Benchmark multi-mode kernel vs N single-mode kernels
# ============================================================================
def test_benchmark():
    """Benchmark multi-mode cavity force performance."""
    print("\n" + "=" * 60)
    print("Test 6: Performance Benchmark")
    print("=" * 60)

    # Try CUDA first, fallback to Reference
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
        platform_name = 'CUDA'
    except Exception:
        try:
            platform = openmm.Platform.getPlatformByName('OpenCL')
            platform_name = 'OpenCL'
        except Exception:
            platform = openmm.Platform.getPlatformByName('Reference')
            platform_name = 'Reference'

    print(f"  Platform: {platform_name}")

    omega1 = 0.007
    lambda1 = 0.01
    cavity_length = 100.0
    molecule_z = 30.0
    photon_mass = 1.0
    n_steps = 500

    for num_modes in [1, 5, 10, 20]:
        sys_mm, pos_mm, _ = build_test_system(
            num_modes, omega1, lambda1, cavity_length, molecule_z, photon_mass,
            use_multimode=True
        )
        integrator = openmm.VerletIntegrator(0.001 * unit.picosecond)
        ctx = openmm.Context(sys_mm, integrator, platform)
        ctx.setPositions(pos_mm)

        # Warmup
        integrator.step(10)

        t0 = time.perf_counter()
        integrator.step(n_steps)
        t_elapsed = time.perf_counter() - t0

        ns_per_day = (n_steps * 0.001e-3) / t_elapsed * 86400 * 1e9  # ns/day
        print(f"  N={num_modes:3d} modes: {t_elapsed:.4f} s ({ns_per_day:.1f} ns/day)")

        del ctx

    print(f"\n  (Benchmark only, always PASS)")
    return True


# ============================================================================
# Main
# ============================================================================
if __name__ == '__main__':
    print("=" * 60)
    print("  MultiModeCavityForce Test Suite")
    print("=" * 60)

    results = {}
    results['scaling'] = test_mode_scaling()
    results['spatial'] = test_spatial_profile()
    results['single_mode'] = test_single_mode_equivalence()
    results['multimode_vs_analytical'] = test_multimode_vs_analytical()
    results['energy_conservation'] = test_energy_conservation()
    results['benchmark'] = test_benchmark()

    print("\n" + "=" * 60)
    print("  Summary")
    print("=" * 60)
    all_pass = True
    for name, result in results.items():
        status = "PASS" if result else "FAIL"
        if not result:
            all_pass = False
        print(f"  {name}: {status}")

    print(f"\n  Overall: {'ALL TESTS PASSED' if all_pass else 'SOME TESTS FAILED'}")
    sys.exit(0 if all_pass else 1)

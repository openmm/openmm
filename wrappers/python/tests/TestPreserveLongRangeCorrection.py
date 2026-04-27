"""Tests for the preserveLongRangeCorrection argument to updateParametersInContext."""

import unittest
import openmm as mm


def _build_grid_positions(grid_size, box_size):
    positions = []
    for i in range(grid_size):
        for j in range(grid_size):
            for k in range(grid_size):
                positions.append(
                    mm.Vec3(
                        i * box_size / grid_size,
                        j * box_size / grid_size,
                        k * box_size / grid_size,
                    )
                )
    return positions


class TestPreserveLongRangeCorrectionNonbonded(unittest.TestCase):
    """
    Test NonbondedForce.updateParametersInContext(preserveLongRangeCorrection=True).

    When only charges are updated (sigma/epsilon unchanged), the dispersion
    correction is the same regardless of whether we preserve or recompute it.
    """

    def setUp(self):
        grid_size = 4
        self.num_particles = grid_size**3
        box_size = grid_size * 0.7
        cutoff = box_size / 3

        def make_context():
            system = mm.System()
            force = mm.NonbondedForce()
            for idx in range(self.num_particles):
                system.addParticle(1.0)
                # Alternate two particle types so the LRC is non-trivial.
                if idx % 2 == 0:
                    force.addParticle(1.0, 1.1, 0.5)
                else:
                    force.addParticle(-1.0, 1.0, 1.0)
            force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
            force.setCutoffDistance(cutoff)
            force.setUseDispersionCorrection(True)
            system.setDefaultPeriodicBoxVectors(
                mm.Vec3(box_size, 0, 0),
                mm.Vec3(0, box_size, 0),
                mm.Vec3(0, 0, box_size),
            )
            system.addForce(force)
            integrator = mm.VerletIntegrator(0.01)
            platform = mm.Platform.getPlatformByName("Reference")
            context = mm.Context(system, integrator, platform)
            context.setPositions(_build_grid_positions(grid_size, box_size))
            return context, force

        self.context1, self.force1 = make_context()
        self.context2, self.force2 = make_context()

    def test_preserve_gives_same_energy_as_recompute(self):
        """
        Updating only charges with preserveLongRangeCorrection=True should
        give the same energy as a normal update.
        """
        # Scale charges by 0.5, leaving sigma and epsilon unchanged.
        for i in range(self.num_particles):
            charge, sigma, epsilon = self.force1.getParticleParameters(i)
            self.force1.setParticleParameters(i, charge * 0.5, sigma, epsilon)
            self.force2.setParticleParameters(i, charge * 0.5, sigma, epsilon)

        self.force1.updateParametersInContext(self.context1)
        self.force2.updateParametersInContext(
            self.context2, preserveLongRangeCorrection=True
        )

        e1 = self.context1.getState(getEnergy=True).getPotentialEnergy()._value
        e2 = self.context2.getState(getEnergy=True).getPotentialEnergy()._value
        self.assertAlmostEqual(e1, e2, places=5)

    def test_energy_actually_changed(self):
        """
        Verify the charge update produced a meaningful energy difference.
        """
        e_initial = self.context1.getState(getEnergy=True).getPotentialEnergy()._value
        for i in range(self.num_particles):
            charge, sigma, epsilon = self.force1.getParticleParameters(i)
            self.force1.setParticleParameters(i, 0.0, sigma, epsilon)
        self.force1.updateParametersInContext(self.context1)
        e_final = self.context1.getState(getEnergy=True).getPotentialEnergy()._value
        self.assertNotAlmostEqual(e_initial, e_final, places=2)


class TestPreserveLongRangeCorrectionCustomNonbonded(unittest.TestCase):
    """
    Test CustomNonbondedForce.updateParametersInContext(preserveLongRangeCorrection=True).

    We use two equal-sized groups of particles and swap their parameter
    assignments.  The LRC coefficient depends only on the count of each
    particle-class pair, which is unchanged by the swap, so both code paths
    should yield the same energy.
    """

    def setUp(self):
        grid_size = 4
        self.num_particles = grid_size**3
        box_size = grid_size * 0.7
        cutoff = box_size / 3

        params_a = [1.1, 0.5]  # sigma, eps for type A
        params_b = [1.0, 1.0]  # sigma, eps for type B

        def make_context():
            system = mm.System()
            force = mm.CustomNonbondedForce(
                "4*eps*((sigma/r)^12-(sigma/r)^6);"
                "sigma=0.5*(sigma1+sigma2);"
                "eps=sqrt(eps1*eps2)"
            )
            force.addPerParticleParameter("sigma")
            force.addPerParticleParameter("eps")
            for idx in range(self.num_particles):
                system.addParticle(1.0)
                force.addParticle(params_a if idx % 2 == 0 else params_b)
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            force.setCutoffDistance(cutoff)
            force.setUseLongRangeCorrection(True)
            system.setDefaultPeriodicBoxVectors(
                mm.Vec3(box_size, 0, 0),
                mm.Vec3(0, box_size, 0),
                mm.Vec3(0, 0, box_size),
            )
            system.addForce(force)
            integrator = mm.VerletIntegrator(0.01)
            platform = mm.Platform.getPlatformByName("Reference")
            context = mm.Context(system, integrator, platform)
            context.setPositions(_build_grid_positions(grid_size, box_size))
            return context, force

        self.context1, self.force1 = make_context()
        self.context2, self.force2 = make_context()
        self.params_a = params_a
        self.params_b = params_b

    def test_preserve_gives_same_energy_as_recompute(self):
        """
        Swapping A/B parameter assignments preserves the LRC coefficient
        (same number of each particle class) but changes the pair energy.
        Both update paths should give the same energy.
        """
        # Swap A and B parameter assignments: even->B, odd->A.
        for i in range(self.num_particles):
            new_params = self.params_b if i % 2 == 0 else self.params_a
            self.force1.setParticleParameters(i, new_params)
            self.force2.setParticleParameters(i, new_params)

        self.force1.updateParametersInContext(self.context1)
        self.force2.updateParametersInContext(
            self.context2, preserveLongRangeCorrection=True
        )

        e1 = self.context1.getState(getEnergy=True).getPotentialEnergy()._value
        e2 = self.context2.getState(getEnergy=True).getPotentialEnergy()._value
        self.assertAlmostEqual(e1, e2, places=5)

    def test_global_param_lrc_caching(self):
        """
        With preserveLongRangeCorrection=True and a global parameter keying the
        LRC cache:
        - Energies must agree with preserveLongRangeCorrection=False at each
          discrete state (catches stale longRangeCorrectionData).
        - Revisiting a state produces the same energy as the first visit
          (confirms cache entries are stored and retrieved correctly).
        """
        grid_size = 4
        num_particles = grid_size**3
        box_size = grid_size * 0.7
        cutoff = box_size / 3

        states = [
            [1.1, 0.5],   # sigma, eps — state 0
            [1.0, 1.0],   # state 1
            [0.9, 0.8],   # state 2
        ]
        positions = _build_grid_positions(grid_size, box_size)

        def make_context(with_global_param):
            system = mm.System()
            force = mm.CustomNonbondedForce(
                "4*eps*((sigma/r)^12-(sigma/r)^6);"
                "sigma=0.5*(sigma1+sigma2);"
                "eps=sqrt(eps1*eps2)"
            )
            force.addPerParticleParameter("sigma")
            force.addPerParticleParameter("eps")
            if with_global_param:
                # stateIndex does not appear in the energy expression but acts as
                # a pure cache key: LRC is recomputed only when stateIndex changes.
                force.addGlobalParameter("stateIndex", 0)
            for _ in range(num_particles):
                system.addParticle(1.0)
                force.addParticle(states[0])
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            force.setCutoffDistance(cutoff)
            force.setUseLongRangeCorrection(True)
            system.setDefaultPeriodicBoxVectors(
                mm.Vec3(box_size, 0, 0),
                mm.Vec3(0, box_size, 0),
                mm.Vec3(0, 0, box_size),
            )
            system.addForce(force)
            integrator = mm.VerletIntegrator(0.01)
            platform = mm.Platform.getPlatformByName("Reference")
            context = mm.Context(system, integrator, platform)
            context.setPositions(positions)
            return context, force

        # Reference context: always recomputes LRC from scratch.
        ctx_ref, force_ref = make_context(with_global_param=False)
        # Cached context: preserves LRC data, keyed by stateIndex.
        ctx_cache, force_cache = make_context(with_global_param=True)

        def energy_ref(state_idx):
            for i in range(num_particles):
                force_ref.setParticleParameters(i, states[state_idx])
            force_ref.updateParametersInContext(ctx_ref)
            return ctx_ref.getState(getEnergy=True).getPotentialEnergy()._value

        def energy_cache(state_idx):
            for i in range(num_particles):
                force_cache.setParticleParameters(i, states[state_idx])
            ctx_cache.setParameter("stateIndex", float(state_idx))
            force_cache.updateParametersInContext(
                ctx_cache, preserveLongRangeCorrection=True
            )
            return ctx_cache.getState(getEnergy=True).getPotentialEnergy()._value

        # First pass: cold cache.  Energies must agree with the reference, which
        # would fail if longRangeCorrectionData were stale.
        ref_energies = [energy_ref(i) for i in range(len(states))]
        cold_energies = [energy_cache(i) for i in range(len(states))]

        for i, (e_ref, e_cold) in enumerate(zip(ref_energies, cold_energies)):
            self.assertAlmostEqual(
                e_ref, e_cold, places=5,
                msg=f"Cold-cache energy mismatch at state {i}"
            )

        # States must have distinct energies so the test is non-trivial.
        self.assertNotAlmostEqual(ref_energies[0], ref_energies[1], places=2)

        # Second pass: warm cache, reversed.  Energies must still match reference.
        for i in reversed(range(len(states))):
            e_warm = energy_cache(i)
            self.assertAlmostEqual(
                ref_energies[i], e_warm, places=5,
                msg=f"Warm-cache energy mismatch at state {i}"
            )

    def test_uniform_type_changes_energy(self):
        """
        Verify the test is non-trivial: making all particles the same type
        (removing class diversity) changes the energy from the mixed initial state.
        """
        e_initial = self.context1.getState(getEnergy=True).getPotentialEnergy()._value
        for i in range(self.num_particles):
            self.force1.setParticleParameters(i, self.params_b)
        self.force1.updateParametersInContext(self.context1)
        e_uniform = self.context1.getState(getEnergy=True).getPotentialEnergy()._value
        self.assertNotAlmostEqual(e_initial, e_uniform, places=2)


if __name__ == "__main__":
    unittest.main()

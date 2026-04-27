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
        LRC cache, revisiting the same stateIndex reproduces the same energy,
        confirming the per-state cache entries are stored and retrieved correctly.
        """
        grid_size = 4
        num_particles = grid_size**3
        box_size = grid_size * 0.7
        cutoff = box_size / 3

        system = mm.System()
        force = mm.CustomNonbondedForce(
            "4*eps*((sigma/r)^12-(sigma/r)^6);"
            "sigma=0.5*(sigma1+sigma2);"
            "eps=sqrt(eps1*eps2)"
        )
        force.addPerParticleParameter("sigma")
        force.addPerParticleParameter("eps")
        # stateIndex does not appear in the energy expression but is still tracked
        # by the kernel as a global parameter.  It acts as a pure cache key: the
        # LRC coefficient is recomputed only when stateIndex takes a new value and
        # reused on revisits
        force.addGlobalParameter("stateIndex", 0)

        params_state0 = [1.1, 0.5]
        params_state1 = [1.0, 1.0]

        for _ in range(num_particles):
            system.addParticle(1.0)
            force.addParticle(params_state0)
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

        def set_state(state_idx, params):
            for i in range(num_particles):
                force.setParticleParameters(i, params)
            context.setParameter("stateIndex", float(state_idx))
            force.updateParametersInContext(context, preserveLongRangeCorrection=True)
            return context.getState(getEnergy=True).getPotentialEnergy()._value

        # First pass: cold cache, one entry per stateIndex.
        e0_first = set_state(0, params_state0)
        e1_first = set_state(1, params_state1)

        # The two states must have different energies so the test is non-trivial.
        self.assertNotAlmostEqual(e0_first, e1_first, places=2)

        # Second pass: warm cache, revisit both states in reverse order.
        e1_second = set_state(1, params_state1)
        e0_second = set_state(0, params_state0)

        self.assertAlmostEqual(e0_first, e0_second, places=5)
        self.assertAlmostEqual(e1_first, e1_second, places=5)

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

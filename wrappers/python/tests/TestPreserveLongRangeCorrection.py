"""Tests for LRC handling in updateParametersInContext."""

import unittest
import openmm as mm


class _LRCTestBase(unittest.TestCase):
    """Common grid geometry and context factory used by all LRC tests."""

    GRID_SIZE = 4
    NUM_PARTICLES = GRID_SIZE**3
    BOX_SIZE = GRID_SIZE * 0.7
    CUTOFF = BOX_SIZE / 3

    @classmethod
    def _positions(cls):
        positions = []
        for i in range(cls.GRID_SIZE):
            for j in range(cls.GRID_SIZE):
                for k in range(cls.GRID_SIZE):
                    positions.append(
                        mm.Vec3(
                            i * cls.BOX_SIZE / cls.GRID_SIZE,
                            j * cls.BOX_SIZE / cls.GRID_SIZE,
                            k * cls.BOX_SIZE / cls.GRID_SIZE,
                        )
                    )
        return positions

    @classmethod
    def _make_context(cls, force):
        """
        Build a periodic System containing *force* and return a ready Context.

        The caller is responsible for adding particles to *force* before calling
        this method; the system is sized to match force.getNumParticles().
        """
        system = mm.System()
        for _ in range(force.getNumParticles()):
            system.addParticle(1.0)
        system.setDefaultPeriodicBoxVectors(
            mm.Vec3(cls.BOX_SIZE, 0, 0),
            mm.Vec3(0, cls.BOX_SIZE, 0),
            mm.Vec3(0, 0, cls.BOX_SIZE),
        )
        system.addForce(force)
        integrator = mm.VerletIntegrator(0.01)
        platform = mm.Platform.getPlatformByName("Reference")
        context = mm.Context(system, integrator, platform)
        context.setPositions(cls._positions())
        return context


class TestAutoDetectLRCNonbonded(_LRCTestBase):
    """
    Test NonbondedForce.updateParametersInContext auto-detection of LJ changes.

    The dispersion correction is recomputed only when sigma or epsilon changes.
    Charge-only updates leave the LRC intact.
    """

    def _make_force(self, sigma=1.0, epsilon=1.0):
        force = mm.NonbondedForce()
        for idx in range(self.NUM_PARTICLES):
            charge = 1.0 if idx % 2 == 0 else -1.0
            force.addParticle(charge, sigma, epsilon)
        force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
        force.setCutoffDistance(self.CUTOFF)
        force.setUseDispersionCorrection(True)
        return force

    def setUp(self):
        self.force1 = self._make_force()
        self.force2 = self._make_force()
        self.context1 = self._make_context(self.force1)
        self.context2 = self._make_context(self.force2)

    def test_charge_only_update_preserves_lrc(self):
        """
        Updating only charges auto-detects no LJ change; energies in two
        otherwise identical contexts must agree after the update.
        """
        for i in range(self.NUM_PARTICLES):
            charge, sigma, epsilon = self.force1.getParticleParameters(i)
            self.force1.setParticleParameters(i, charge * 0.5, sigma, epsilon)
            self.force2.setParticleParameters(i, charge * 0.5, sigma, epsilon)

        self.force1.updateParametersInContext(self.context1)
        self.force2.updateParametersInContext(self.context2)

        e1 = self.context1.getState(getEnergy=True).getPotentialEnergy()._value
        e2 = self.context2.getState(getEnergy=True).getPotentialEnergy()._value
        self.assertAlmostEqual(e1, e2, places=5)

    def test_lj_update_recomputes_lrc(self):
        """
        Updating sigma/epsilon must trigger an LRC recompute; the resulting
        energy must match a freshly initialized context with those parameters.
        """
        sigma_new, epsilon_new = 1.5, 0.8

        # Reference: context built from scratch with the new LJ parameters.
        force_ref = self._make_force(sigma=sigma_new, epsilon=epsilon_new)
        ctx_ref = self._make_context(force_ref)
        e_ref = ctx_ref.getState(getEnergy=True).getPotentialEnergy()._value

        # Update context: start from the default parameters and update in-place.
        for i in range(self.force1.getNumParticles()):
            charge, _, _ = self.force1.getParticleParameters(i)
            self.force1.setParticleParameters(i, charge, sigma_new, epsilon_new)
        self.force1.updateParametersInContext(self.context1)
        e_updated = self.context1.getState(getEnergy=True).getPotentialEnergy()._value

        self.assertAlmostEqual(e_ref, e_updated, places=5)

    def test_energy_actually_changed(self):
        """Verify the charge update produced a meaningful energy difference."""
        e_initial = self.context1.getState(getEnergy=True).getPotentialEnergy()._value
        for i in range(self.NUM_PARTICLES):
            charge, sigma, epsilon = self.force1.getParticleParameters(i)
            self.force1.setParticleParameters(i, 0.0, sigma, epsilon)
        self.force1.updateParametersInContext(self.context1)
        e_final = self.context1.getState(getEnergy=True).getPotentialEnergy()._value
        self.assertNotAlmostEqual(e_initial, e_final, places=2)


class TestPreserveLongRangeCorrectionCustomNonbonded(_LRCTestBase):
    """
    Test CustomNonbondedForce.updateParametersInContext(preserveLongRangeCorrection=True).

    We use two equal-sized groups of particles and swap their parameter
    assignments.  The LRC coefficient depends only on the count of each
    particle-class pair, which is unchanged by the swap, so both code paths
    should yield the same energy.
    """

    PARAMS_A = [1.1, 0.5]  # sigma, eps for type A
    PARAMS_B = [1.0, 1.0]  # sigma, eps for type B

    LJ_EXPRESSION = (
        "4*eps*((sigma/r)^12-(sigma/r)^6);"
        "sigma=0.5*(sigma1+sigma2);"
        "eps=sqrt(eps1*eps2)"
    )

    def _make_force(self, initial_params=None, global_param=None):
        """
        Return a CutoffPeriodic LJ CustomNonbondedForce with LRC enabled.

        Particles alternate between PARAMS_A and PARAMS_B unless *initial_params*
        is given, in which case all particles get that single parameter set.
        An optional *global_param* name is added as a pure cache key.
        """
        force = mm.CustomNonbondedForce(self.LJ_EXPRESSION)
        force.addPerParticleParameter("sigma")
        force.addPerParticleParameter("eps")
        if global_param is not None:
            force.addGlobalParameter(global_param, 0.0)
        for idx in range(self.NUM_PARTICLES):
            if initial_params is not None:
                force.addParticle(initial_params)
            else:
                force.addParticle(self.PARAMS_A if idx % 2 == 0 else self.PARAMS_B)
        force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        force.setCutoffDistance(self.CUTOFF)
        force.setUseLongRangeCorrection(True)
        return force

    def setUp(self):
        self.force1 = self._make_force()
        self.force2 = self._make_force()
        self.context1 = self._make_context(self.force1)
        self.context2 = self._make_context(self.force2)

    def test_preserve_gives_same_energy_as_recompute(self):
        """
        Swapping A/B parameter assignments preserves the LRC coefficient
        (same number of each particle class) but changes the pair energy.
        Both update paths should give the same energy.
        """
        for i in range(self.NUM_PARTICLES):
            new_params = self.PARAMS_B if i % 2 == 0 else self.PARAMS_A
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
        states = [
            [1.1, 0.5],  # sigma, eps — state 0
            [1.0, 1.0],  # state 1
            [0.9, 0.8],  # state 2
        ]

        # Reference context: always recomputes LRC from scratch.
        force_ref = self._make_force(initial_params=states[0])
        ctx_ref = self._make_context(force_ref)

        # Cached context: preserves LRC data, keyed by stateIndex.
        # stateIndex does not appear in the energy expression; it is a pure cache key.
        force_cache = self._make_force(
            initial_params=states[0], global_param="stateIndex"
        )
        ctx_cache = self._make_context(force_cache)

        def energy_ref(state_idx):
            for i in range(self.NUM_PARTICLES):
                force_ref.setParticleParameters(i, states[state_idx])
            force_ref.updateParametersInContext(ctx_ref)
            return ctx_ref.getState(getEnergy=True).getPotentialEnergy()._value

        def energy_cache(state_idx):
            for i in range(self.NUM_PARTICLES):
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
                e_ref, e_cold, places=5, msg=f"Cold-cache energy mismatch at state {i}"
            )

        # States must have distinct energies so the test is non-trivial.
        self.assertNotAlmostEqual(ref_energies[0], ref_energies[1], places=2)

        # Second pass: warm cache, reversed.  Energies must still match reference.
        for i in reversed(range(len(states))):
            e_warm = energy_cache(i)
            self.assertAlmostEqual(
                ref_energies[i],
                e_warm,
                places=5,
                msg=f"Warm-cache energy mismatch at state {i}",
            )

    def test_uniform_type_changes_energy(self):
        """
        Verify the test is non-trivial: making all particles the same type
        (removing class diversity) changes the energy from the mixed initial state.
        """
        e_initial = self.context1.getState(getEnergy=True).getPotentialEnergy()._value
        for i in range(self.NUM_PARTICLES):
            self.force1.setParticleParameters(i, self.PARAMS_B)
        self.force1.updateParametersInContext(self.context1)
        e_uniform = self.context1.getState(getEnergy=True).getPotentialEnergy()._value
        self.assertNotAlmostEqual(e_initial, e_uniform, places=2)


if __name__ == "__main__":
    unittest.main()

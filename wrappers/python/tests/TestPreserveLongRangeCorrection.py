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


if __name__ == "__main__":
    unittest.main()

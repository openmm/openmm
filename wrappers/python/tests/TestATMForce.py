import unittest
from openmm import *
from openmm.app import *
from openmm.unit import *

class TestATMForce(unittest.TestCase):
    """Tests the ATMForce"""

    def test2ParticlesNonbonded(self):
        """Test for a Nonbonded force previously added to the System"""
        system = System()
        system.addParticle(1.0)
        system.addParticle(1.0)

        nbforce = NonbondedForce();
        nbforce.addParticle( 1.0, 1.0, 1.0)
        nbforce.addParticle(-1.0, 1.0, 1.0)

        system.addForce(nbforce)

        atmforce = ATMForce(0.5, 0.5, 0, 0, 0,   0, 0, 0,  1.0)
        atmforce.addParticle(Vec3(0., 0., 0.))
        atmforce.addParticle(Vec3(1., 0., 0.))

        atmforce.addForce(copy.copy(nbforce))
        system.removeForce(0)
        system.addForce(atmforce)

        integrator = VerletIntegrator(1.0)
        platform = Platform.getPlatform('Reference')
        context = Context(system, integrator, platform)

        positions = []
        positions.append(Vec3(0., 0., 0.))
        positions.append(Vec3(1., 0., 0.))
        context.setPositions(positions)

        state = context.getState(getEnergy = True, getForces = True)
        epot = state.getPotentialEnergy()
        
        (u1, u0, energy) = atmforce.getPerturbationEnergy(context)
        epert = u1 - u0

        #print("Potential energy = ", epot)
        #print("ATM perturbation energy = ", epert)
        
        epot_expected = -104.2320*kilojoules_per_mole
        epert_expected = 69.4062*kilojoules_per_mole
        assert( abs(epot-epot_expected) < 1.e-3*kilojoules_per_mole )
        assert( abs(epert-epert_expected) < 1.e-3*kilojoules_per_mole )

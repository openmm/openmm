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
        atmforce.addParticle()
        p = atmforce.addParticle()
        atmforce.setParticleTransformation(p, FixedDisplacement(Vec3(1., 0., 0.)))

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
        assert isinstance(atmforce.getParticleTransformation(0), FixedDisplacement)

    def test3ParticlesNonbondedSwap(self):
        """Test coordinate swap"""
        system = System()
        system.addParticle(1.0)
        system.addParticle(1.0)
        system.addParticle(1.0)

        nbforce = NonbondedForce();
        nbforce.addParticle( 0.0, 1.0, 1.0)
        nbforce.addParticle( 1.0, 1.0, 1.0)
        nbforce.addParticle(-1.0, 1.0, 1.0)

        atmforce = ATMForce(0.5, 0.5, 0, 0, 0,   0, 0, 0,  1.0)
        atmforce.addParticle() #particle 0 is not displaced

        #particle 1's coordinate is swapped with 2
        p = atmforce.addParticle()
        atmforce.setParticleTransformation(p, ParticleOffsetDisplacement(2,  1))

        #particle 2's coordinate is swapped with 1
        p = atmforce.addParticle()
        atmforce.setParticleTransformation(p, ParticleOffsetDisplacement(1,  2))

        atmforce.addForce(nbforce)
        system.addForce(atmforce)

        integrator = VerletIntegrator(1.0)
        platform = Platform.getPlatform('Reference')
        context = Context(system, integrator, platform)

        positions = []
        positions.append(Vec3( 0., 0., 0.))
        positions.append(Vec3( 1., 0., 0.))
        positions.append(Vec3(-1., 0., 0.))
        context.setPositions(positions)

        state = context.getState(getEnergy = True, getForces = True)
        epot = state.getPotentialEnergy()

        (u1, u0, energy) = atmforce.getPerturbationEnergy(context)
        epert = u1 - u0

        #print("Potential energy = ", epot)
        #print("ATM perturbation energy = ", epert)

        epot_expected = -69.52925*kilojoules_per_mole
        epert_expected = 0*kilojoules_per_mole
        assert( abs(epot-epot_expected) < 1.e-3*kilojoules_per_mole )
        assert( abs(epert-epert_expected) < 1.e-3*kilojoules_per_mole )

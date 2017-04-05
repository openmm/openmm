import unittest
import tempfile
from datetime import datetime, timedelta
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
import math, random

class TestIntegrators(unittest.TestCase):
    """Test Python Integrator classes"""

    def testMTSIntegratorExplicit(self):
        """Test the MTS integrator on an explicit solvent system"""
        # Create a periodic solvated system with PME
        pdb = PDBFile('systems/alanine-dipeptide-explicit.pdb')
        ff = ForceField('amber99sbildn.xml', 'tip3p.xml')
        system = ff.createSystem(pdb.topology, cutoffMethod=PME)

        # Split forces into groups
        for force in system.getForces():
            if force.__class__.__name__ == 'NonbondedForce':
                force.setForceGroup(1)
                force.setReciprocalSpaceForceGroup(2)
            else:
                force.setForceGroup(0)

        # Create an integrator
        integrator = MTSIntegrator(4*femtoseconds, [(2,1), (1,2), (0,8)])

        # Run a few steps of dynamics
        context = Context(system, integrator)
        context.setPositions(pdb.positions)
        integrator.step(10)

        # Ensure energy is well-behaved.
        state = context.getState(getEnergy=True)
        if not (state.getPotentialEnergy() / kilojoules_per_mole < 0.0):
            raise Exception('Potential energy of alanine dipeptide system with MTS integrator is blowing up: %s' % str(state.getPotentialEnergy()))

    def testMTSIntegratorConstraints(self):
        """Test the MTS integrator energy conservation on a system of constrained particles with no inner force (just constraints)"""

        # Create a constrained test system
        numParticles = 8
        numConstraints = 5
        system = System()
        force = NonbondedForce()
        for i in range(numParticles):
            system.addParticle(5.0 if i%2==0 else 10.0)
            force.addParticle((0.2 if i%2==0 else -0.2), 0.5, 5.0);
        system.addConstraint(0, 1, 1.0);
        system.addConstraint(1, 2, 1.0);
        system.addConstraint(2, 3, 1.0);
        system.addConstraint(4, 5, 1.0);
        system.addConstraint(6, 7, 1.0);
        system.addForce(force)

        # Create integrator where inner timestep just evaluates constraints
        integrator = MTSIntegrator(1*femtoseconds, [(1,1), (0,4)])
        integrator.setConstraintTolerance(1e-5);

        positions = [ (i/2., (i+1)/2., 0.) for i in range(numParticles) ]
        velocities = [ (random.random()-0.5, random.random()-0.5, random.random()-0.5) for i in range(numParticles) ]

        # Create Context
        platform = Platform.getPlatformByName('Reference')
        context = Context(system, integrator, platform)
        context.setPositions(positions)
        context.setVelocities(velocities)
        context.applyConstraints(1e-5)

        # Simulate it and see whether the constraints remain satisfied.
        CONSTRAINT_RELATIVE_TOLERANCE = 1.e-4 # relative constraint violation tolerance
        ENERGY_RELATIVE_TOLERANCE = 1.e-2 # relative energy violation tolerance
        for i in range(1000):
            state = context.getState(getPositions=True, getEnergy=True)
            positions = state.getPositions()
            for j in range(numConstraints):
                [particle1, particle2, constraint_distance] = system.getConstraintParameters(j)
                current_distance = 0.0 * nanometers**2
                for k in range(3):
                    current_distance += (positions[particle1][k] - positions[particle2][k])**2
                current_distance = sqrt(current_distance)
                # Fail test if outside of relative tolerance
                relative_violation = (current_distance - constraint_distance) / constraint_distance
                if (relative_violation > CONSTRAINT_RELATIVE_TOLERANCE):
                    raise Exception('Constrained distance is violated by relative tolerance of %f (constraint %s actual %s)' % (relative_violation, str(constraint_distance), str(current_distance)))
            # Check total energy
            total_energy = state.getPotentialEnergy() + state.getKineticEnergy()
            if (i == 1):
                initial_energy = total_energy
            elif (i > 1):
                relative_violation = abs((total_energy - initial_energy) / initial_energy)
                if (relative_violation > ENERGY_RELATIVE_TOLERANCE):
                    raise Exception('Total energy is violated by relative tolerance of %f on step %d (initial %s final %s)' % (relative_violation, i, str(initial_energy), str(total_energy)))
            # Take a step
            integrator.step(1)

    def testBadGroups(self):
        """Test the MTS integrator with bad force group substeps."""
        # Create a periodic solvated system with PME
        pdb = PDBFile('systems/alanine-dipeptide-explicit.pdb')
        ff = ForceField('amber99sbildn.xml', 'tip3p.xml')
        system = ff.createSystem(pdb.topology, cutoffMethod=PME)

        # Split forces into groups
        for force in system.getForces():
            if force.__class__.__name__ == 'NonbondedForce':
                force.setForceGroup(1)
                force.setReciprocalSpaceForceGroup(2)
            else:
                force.setForceGroup(0)

        with self.assertRaises(ValueError):
            # Create an integrator
            integrator = MTSIntegrator(4*femtoseconds, [(2,1), (1,3), (0,8)])

            # Run a few steps of dynamics
            context = Context(system, integrator)
            context.setPositions(pdb.positions)
            integrator.step(10)

if __name__ == '__main__':
    unittest.main()

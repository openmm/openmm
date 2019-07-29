import unittest
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

class TestMetadynamics(unittest.TestCase):
    """Test the Metadynamics class"""

    def testHarmonicOscillator(self):
        """Test running metadynamics on a harmonic oscillator."""
        system = System()
        system.addParticle(1.0)
        system.addParticle(1.0)
        force = HarmonicBondForce()
        force.addBond(0, 1, 1.0, 100000.0)
        system.addForce(force)
        cv = CustomBondForce('r')
        cv.addBond(0, 1)
        bias = BiasVariable(cv, 0.94, 1.06, 0.02)
        meta = Metadynamics(system, [bias], 300*kelvin, 2.0, 5.0, 10)
        integrator = LangevinIntegrator(300*kelvin, 10/picosecond, 0.001*picosecond)
        topology = Topology()
        chain = topology.addChain()
        residue = topology.addResidue('H2', chain)
        topology.addAtom('H1', element.hydrogen, residue)
        topology.addAtom('H2', element.hydrogen, residue)
        simulation = Simulation(topology, system, integrator, Platform.getPlatformByName('Reference'))
        simulation.context.setPositions([Vec3(0, 0, 0), Vec3(1, 0, 0)])
        meta.step(simulation, 200000)
        fe = meta.getFreeEnergy()
        center = bias.gridWidth//2
        fe -= fe[center]

        # Energies should be reasonably well converged over the central part of the range.

        for i in range(center-3, center+4):
            r = bias.minValue + i*(bias.maxValue-bias.minValue)/(bias.gridWidth-1)
            e = 0.5*100000.0*(r-1.0)**2*kilojoules_per_mole
            assert abs(fe[i]-e) < 1.0*kilojoules_per_mole
import unittest
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

class TestSimulatedTempering(unittest.TestCase):
    """Test the SimulatedTempering class"""

    def testHarmonicOscillator(self):
        """Test running simulated tempering on a harmonic oscillator."""
        system = System()
        system.addParticle(1.0)
        system.addParticle(1.0)
        force = HarmonicBondForce()
        force.addBond(0, 1, 1.0, 1000.0)
        system.addForce(force)
        integrator = LangevinIntegrator(300*kelvin, 10/picosecond, 0.001*picosecond)
        topology = Topology()
        chain = topology.addChain()
        residue = topology.addResidue('H2', chain)
        topology.addAtom('H1', element.hydrogen, residue)
        topology.addAtom('H2', element.hydrogen, residue)
        simulation = Simulation(topology, system, integrator, Platform.getPlatformByName('Reference'))
        st = SimulatedTempering(simulation, numTemperatures=10, minTemperature=200*kelvin, maxTemperature=400*kelvin, tempChangeInterval=5, reportInterval=10000)
        self.assertEqual(10, len(st.temperatures))
        self.assertEqual(200*kelvin, st.temperatures[0])
        self.assertEqual(400*kelvin, st.temperatures[-1])
        simulation.context.setPositions([Vec3(0, 0, 0), Vec3(1, 0, 0)])
        
        # Run for a little while to let the weights stabilize.
        
        st.step(10000)
        
        # Run for a while and record the temperatures and distances.
        
        distances = [[] for i in range(10)]
        count = 0
        for i in range(7000):
            st.step(5)
            pos = simulation.context.getState(getPositions=True).getPositions().value_in_unit(nanometers)
            r = norm(pos[0]-pos[1])
            distances[st.currentTemperature].append(r)
            count += 1

        # Check that it spent roughly equal time at each temperature, and that the distributions
        # are correct.

        for d, t in zip(distances, st.temperatures):
            n = len(d)
            assert count/20 < n < count/5
            meanDist = sum(d)/n
            assert 0.97 < meanDist < 1.03
            meanEnergy = sum([0.5*1000*(r-1)**2 for r in d])/n
            expectedEnergy = (0.5*MOLAR_GAS_CONSTANT_R*t).value_in_unit(kilojoules_per_mole)
            self.assertAlmostEqual(expectedEnergy, meanEnergy, delta=expectedEnergy*0.3)
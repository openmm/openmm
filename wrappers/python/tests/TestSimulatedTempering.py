import unittest
from openmm import *
from openmm.app import *
from openmm.unit import *

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
        simulation = Simulation(topology, system, integrator, Platform.getPlatform('Reference'))
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
        for i in range(20000):
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


    def testHarmonicOscillatorNPT(self):
        """Test running simulated tempering on a harmonic oscillator with a Monte Carlo barostat."""
        system = System()
        system.addParticle(1.0)
        system.addParticle(1.0)
        force = HarmonicBondForce()
        force.addBond(0, 1, 1.0, 1000.0)
        system.addForce(force)
        mcbarostat = MonteCarloBarostat(1*bar, 100*kelvin, 2)
        system.addForce(mcbarostat)
        force.setUsesPeriodicBoundaryConditions(True)
        integrator = LangevinIntegrator(100*kelvin, 10/picosecond, 0.001*picosecond)
        topology = Topology()
        chain = topology.addChain()
        residue = topology.addResidue('H2', chain)
        topology.addAtom('H1', element.hydrogen, residue)
        topology.addAtom('H2', element.hydrogen, residue)
        simulation = Simulation(topology, system, integrator, Platform.getPlatform('Reference'))
        simulation.context.setPositions([Vec3(0, 0, 0), Vec3(1, 0, 0)])

        # Check the temperature is correct before creating the simulated tempering simulation
        
        self.assertEqual(100*kelvin, integrator.getTemperature())
        self.assertEqual(100*kelvin, simulation.context.getParameter('MonteCarloTemperature')*kelvin)

        st = SimulatedTempering(simulation, numTemperatures=10, minTemperature=200*kelvin, maxTemperature=400*kelvin, tempChangeInterval=4, reportInterval=10000)

        # Check the temperatures of the integrator and barostat are set to the value of minTemperature

        self.assertEqual(10, len(st.temperatures))
        self.assertEqual(200*kelvin, st.temperatures[0])
        self.assertEqual(400*kelvin, st.temperatures[-1])
        self.assertEqual(200*kelvin, integrator.getTemperature())
        self.assertEqual(200*kelvin, simulation.context.getParameter('MonteCarloTemperature')*kelvin)

        # Let the simulation run and assert at every step T(mcbarostat) == T(integrator) == T(tempering)

        for i in range(100):
            st.step(2)
            self.assertEqual(st.temperatures[st.currentTemperature], integrator.getTemperature())
            self.assertEqual(st.temperatures[st.currentTemperature], simulation.context.getParameter('MonteCarloTemperature')*kelvin)



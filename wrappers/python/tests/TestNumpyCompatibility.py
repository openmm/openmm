import unittest
import numpy as np

from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit

class TestNumpyCompatibility(unittest.TestCase):
    
    def setUp(self):
        prmtop = app.AmberPrmtopFile('systems/water-box-216.prmtop')
    
        system = prmtop.createSystem(nonbondedMethod=app.PME, 
                                         nonbondedCutoff=0.9*unit.nanometers,
                                         constraints=app.HBonds, rigidWater=True, 
                                         ewaldErrorTolerance=0.0005)
        integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 
                                           2.0*unit.femtoseconds)
        self.simulation = app.Simulation(prmtop.topology, system, integrator,
                                    mm.Platform.getPlatformByName('Reference'))
    
    
    def test_setPositions(self):
        n_particles = self.simulation.context.getSystem().getNumParticles()
        input = np.random.randn(n_particles, 3)
        self.simulation.context.setPositions(input)
        output = self.simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
    
        np.testing.assert_array_almost_equal(input, output)
    
    def test_setPositions_units(self):
        n_particles = self.simulation.context.getSystem().getNumParticles()
        input = unit.Quantity(np.random.randn(n_particles, 3), unit.angstroms)
        self.simulation.context.setPositions(input)
        output = self.simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
    
        np.testing.assert_array_almost_equal(input.value_in_unit(unit.nanometers), output.value_in_unit(unit.nanometers))
    
    
    def test_setVelocities(self):
        n_particles = self.simulation.context.getSystem().getNumParticles()
        input = np.random.randn(n_particles, 3)
        self.simulation.context.setVelocities(input)
        output = self.simulation.context.getState(getVelocities=True).getVelocities(asNumpy=True)
    
        np.testing.assert_array_almost_equal(input, output)
    
    def test_setVelocities_units(self):
        n_particles = self.simulation.context.getSystem().getNumParticles()
        input = unit.Quantity(np.random.randn(n_particles, 3), unit.angstroms / unit.femtoseconds)
        self.simulation.context.setVelocities(input)
        output = self.simulation.context.getState(getVelocities=True).getVelocities(asNumpy=True)
    
        np.testing.assert_array_almost_equal(input.value_in_unit(unit.angstroms / unit.femtoseconds),
                                             output.value_in_unit(unit.angstroms / unit.femtoseconds))
    
    
    def test_tabulatedFunction(self):
        f = mm.CustomNonbondedForce('g(r)')
    
        r = np.linspace(0,10)
        g_of_r = np.sin(r)
    
        indx = f.addFunction('g', g_of_r, np.min(r), np.max(r)) 
    
        name, g_of_r_out, min_r_out, max_r_out = f.getFunctionParameters(indx)
    
        np.testing.assert_array_almost_equal(g_of_r, np.asarray(g_of_r_out))
        assert min_r_out == np.min(r)
        assert max_r_out == np.max(r)
    
    def test_CMAP(self):
        f = mm.CMAPTorsionForce()
        energy = np.random.randn(10*10)
        f.addMap(10, energy)
        size, energy_out = f.getMapParameters(0)
    
        assert size == 10
        np.testing.assert_array_almost_equal(energy, np.asarray(energy_out))

if __name__ == '__main__':
    unittest.main()

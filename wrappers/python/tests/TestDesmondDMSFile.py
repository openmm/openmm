import os
import unittest
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem

class TestDesmondDMSFile(unittest.TestCase):
    def setUp(self):
        """Set up the tests by loading the input files."""

        # alanine dipeptide with explicit water
        path = os.path.join(os.path.dirname(__file__), 'systems/alanine-dipeptide-explicit-amber99SBILDN-tip3p.dms')
        self.dms = DesmondDMSFile(path)

        #nabumetone OPLS force field
        path = os.path.join(os.path.dirname(__file__), 'systems/bcd-nabumetone_lig.dms')
        self.dms_opls1 = DesmondDMSFile(path)
        
        #beta-cyclodextrin/nabumetone complex OPLS force field
        path1 = os.path.join(os.path.dirname(__file__), 'systems/bcd-nabumetone_lig.dms')
        path2 = os.path.join(os.path.dirname(__file__), 'systems/bcd-nabumetone_rcpt.dms')
        self.dms_opls2 = DesmondDMSFile([path1,path2])
    
    def test_NonbondedMethod(self):
        """Test all six options for the nonbondedMethod parameter."""

        methodMap = {NoCutoff:NonbondedForce.NoCutoff, 
                     CutoffNonPeriodic:NonbondedForce.CutoffNonPeriodic, 
                     CutoffPeriodic:NonbondedForce.CutoffPeriodic, 
                     Ewald:NonbondedForce.Ewald,
                     PME:NonbondedForce.PME,
                     LJPME:NonbondedForce.LJPME}
        for method in methodMap:
            system = self.dms.createSystem(nonbondedMethod=method)
            forces = system.getForces()
            self.assertTrue(any(isinstance(f, NonbondedForce) and 
                                f.getNonbondedMethod()==methodMap[method] 
                                for f in forces))
                                
    def test_Cutoff(self):
        """Test to make sure the nonbondedCutoff parameter is passed correctly."""

        for method in [CutoffNonPeriodic, CutoffPeriodic, Ewald, PME, LJPME]:
            system = self.dms.createSystem(nonbondedMethod=method, 
                                            nonbondedCutoff=2*nanometer)
            cutoff_distance = 0.0*nanometer
            cutoff_check = 2.0*nanometer
            for force in system.getForces():
                if isinstance(force, NonbondedForce):
                    cutoff_distance = force.getCutoffDistance()
            self.assertEqual(cutoff_distance, cutoff_check)

    def test_EwaldErrorTolerance(self):
        """Test to make sure the ewaldErrorTolerance parameter is passed correctly."""

        for method in [Ewald, PME, LJPME]:
            system = self.dms.createSystem(nonbondedMethod=method,
                                            ewaldErrorTolerance=1e-6)
            tolerance = 0
            tolerance_check = 1e-6
            for force in system.getForces():
                if isinstance(force, NonbondedForce):
                    tolerance = force.getEwaldErrorTolerance()
            self.assertEqual(tolerance, tolerance_check)

    def test_RemoveCMMotion(self):
        """Test both options (True and False) for the removeCMMotion parameter."""

        for b in [True, False]:
            system = self.dms.createSystem(removeCMMotion=b)
            self.assertEqual(any(isinstance(f, CMMotionRemover) for f in system.getForces()), b)

    def test_HydrogenMass(self):
        """Test that altering the mass of hydrogens works correctly."""
        
        topology = self.dms.topology
        hydrogenMass = 4*amu
        system1 = self.dms.createSystem()
        system2 = self.dms.createSystem(hydrogenMass=hydrogenMass)
        for atom in topology.atoms():
            if atom.element == elem.hydrogen:
                self.assertNotEqual(hydrogenMass, system1.getParticleMass(atom.index))
                self.assertEqual(hydrogenMass, system2.getParticleMass(atom.index))
        totalMass1 = sum([system1.getParticleMass(i) for i in range(system1.getNumParticles())]).value_in_unit(amu)
        totalMass2 = sum([system2.getParticleMass(i) for i in range(system2.getNumParticles())]).value_in_unit(amu)
        self.assertAlmostEqual(totalMass1, totalMass2)

    def test_OPLS1(self):
        """Test OPLS potential energy of nabumetone """
        system = self.dms_opls1.createSystem(nonbondedMethod=NoCutoff, OPLS = True)
        integrator = LangevinIntegrator(300*kelvin, 1.0/picosecond, 0.0005*picoseconds)
        context = Context(system, integrator,Platform.getPlatformByName('Reference'))
        context.setPositions(self.dms_opls1.positions)
        ene = context.getState(getEnergy=True).getPotentialEnergy()
        self.assertAlmostEqual(ene.value_in_unit(kilojoules_per_mole), 94.21, places=2)
        
    def test_OPLS2(self):
        """Test OPLS potential energy of beta-cyclodextrin/nabumetone complex """
        system = self.dms_opls2.createSystem(nonbondedMethod=NoCutoff, OPLS = True)
        integrator = LangevinIntegrator(300*kelvin, 1.0/picosecond, 0.0005*picoseconds)
        context = Context(system, integrator,Platform.getPlatformByName('Reference'))
        context.setPositions(self.dms_opls2.positions)
        ene = context.getState(getEnergy=True).getPotentialEnergy()
        self.assertAlmostEqual(ene.value_in_unit(kilojoules_per_mole), 903.04, places=2)

if __name__ == '__main__':
    unittest.main()

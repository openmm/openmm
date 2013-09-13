import unittest
from validateConstraints import *
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem

class TestGromacsTopFile(unittest.TestCase):

    """Test the GromacsTopFile.createSystem() method."""
 
    def setUp(self):
        """Set up the tests by loading the input files."""

        # alanine dipeptide with explicit water
        self.top1 = GromacsTopFile('systems/explicit.top', unitCellDimensions=Vec3(6.223, 6.223, 6.223)*nanometers)
        # alanine dipeptide with implicit water
        self.top2 = GromacsTopFile('systems/implicit.top')

    def test_NonbondedMethod(self):
        """Test all five options for the nonbondedMethod parameter."""

        methodMap = {NoCutoff:NonbondedForce.NoCutoff, 
                     CutoffNonPeriodic:NonbondedForce.CutoffNonPeriodic, 
                     CutoffPeriodic:NonbondedForce.CutoffPeriodic, 
                     Ewald:NonbondedForce.Ewald, PME: NonbondedForce.PME}
        for method in methodMap:
            system = self.top1.createSystem(nonbondedMethod=method)
            forces = system.getForces()
            self.assertTrue(any(isinstance(f, NonbondedForce) and 
                                f.getNonbondedMethod()==methodMap[method] 
                                for f in forces))

    def test_Cutoff(self):
        """Test to make sure the nonbondedCutoff parameter is passed correctly."""

        for method in [CutoffNonPeriodic, CutoffPeriodic, Ewald, PME]:
            system = self.top1.createSystem(nonbondedMethod=method, 
                                               nonbondedCutoff=2*nanometer, 
                                               constraints=HBonds)
            cutoff_distance = 0.0*nanometer
            cutoff_check = 2.0*nanometer
            for force in system.getForces():
                if isinstance(force, NonbondedForce):
                    cutoff_distance = force.getCutoffDistance()
            self.assertEqual(cutoff_distance, cutoff_check)

    def test_EwaldErrorTolerance(self):
        """Test to make sure the ewaldErrorTolerance parameter is passed correctly."""

        for method in [Ewald, PME]:
            system = self.top1.createSystem(nonbondedMethod=method, 
                                               ewaldErrorTolerance=1e-6, 
                                               constraints=HBonds)
            tolerance = 0
            tolerance_check = 1e-6
            for force in system.getForces():
                if isinstance(force, NonbondedForce):
                    tolerance = force.getEwaldErrorTolerance()
            self.assertEqual(tolerance, tolerance_check)

    def test_RemoveCMMotion(self):
        """Test both options (True and False) for the removeCMMotion parameter."""

        for b in [True, False]:
            system = self.top1.createSystem(removeCMMotion=b)
            self.assertEqual(any(isinstance(f, CMMotionRemover) for f in system.getForces()), b)

    def test_RigidWaterAndConstraints(self):
        """Test all eight options for the constraints and rigidWater parameters."""

        topology = self.top1.topology
        for constraints_value in [None, HBonds, AllBonds, HAngles]:
            for rigidWater_value in [True, False]:
                system = self.top1.createSystem(constraints=constraints_value, 
                                                   rigidWater=rigidWater_value)
                validateConstraints(self, topology, system, 
                                    constraints_value, rigidWater_value)

    def test_ImplicitSolvent(self):
        """Test implicit solvent using the implicitSolvent parameter.

        """
        system = self.top2.createSystem(implicitSolvent=OBC2)
        self.assertTrue(any(isinstance(f, GBSAOBCForce) for f in system.getForces()))

    def test_ImplicitSolventParameters(self):
        """Test that solventDielectric and soluteDielectric are passed correctly.

        """
        system = self.top2.createSystem(implicitSolvent=OBC2, 
                                           solventDielectric=50.0, 
                                           soluteDielectric = 0.9)
        found_matching_solvent_dielectric=False
        found_matching_solute_dielectric=False
        for force in system.getForces():
            if isinstance(force, GBSAOBCForce):
                if force.getSolventDielectric() == 50.0:
                    found_matching_solvent_dielectric = True
                if force.getSoluteDielectric() == 0.9:
                    found_matching_solute_dielectric = True
            if isinstance(force, NonbondedForce):
                self.assertEqual(force.getReactionFieldDielectric(), 1.0)
        self.assertTrue(found_matching_solvent_dielectric and 
                        found_matching_solute_dielectric)

    def test_HydrogenMass(self):
        """Test that altering the mass of hydrogens works correctly."""
        
        topology = self.top1.topology
        hydrogenMass = 4*amu
        system1 = self.top1.createSystem()
        system2 = self.top1.createSystem(hydrogenMass=hydrogenMass)
        for atom in topology.atoms():
            if atom.element == elem.hydrogen:
                self.assertNotEqual(hydrogenMass, system1.getParticleMass(atom.index))
                self.assertEqual(hydrogenMass, system2.getParticleMass(atom.index))
        totalMass1 = sum([system1.getParticleMass(i) for i in range(system1.getNumParticles())]).value_in_unit(amu)
        totalMass2 = sum([system2.getParticleMass(i) for i in range(system2.getNumParticles())]).value_in_unit(amu)
        self.assertAlmostEqual(totalMass1, totalMass2)

if __name__ == '__main__':
    unittest.main()


import unittest
from validateConstraints import *
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem

class TestAmberPrmtopFile(unittest.TestCase):

    """Test the AmberPrmtopFile.createSystem() method."""
 
    def setUp(self):
        """Set up the tests by loading the input files."""

        # alanine dipeptide with explicit water
        self.prmtop1 = AmberPrmtopFile('systems/alanine-dipeptide-explicit.prmtop')
        # alanine dipeptide with implicit water
        self.prmtop2 = AmberPrmtopFile('systems/alanine-dipeptide-implicit.prmtop')

    def test_NonbondedMethod(self):
        """Test all five options for the nonbondedMethod parameter."""

        methodMap = {NoCutoff:NonbondedForce.NoCutoff, 
                     CutoffNonPeriodic:NonbondedForce.CutoffNonPeriodic, 
                     CutoffPeriodic:NonbondedForce.CutoffPeriodic, 
                     Ewald:NonbondedForce.Ewald, PME: NonbondedForce.PME}
        for method in methodMap:
            system = self.prmtop1.createSystem(nonbondedMethod=method)
            forces = system.getForces()
            self.assertTrue(any(isinstance(f, NonbondedForce) and 
                                f.getNonbondedMethod()==methodMap[method] 
                                for f in forces))

    def test_Cutoff(self):
        """Test to make sure the nonbondedCutoff parameter is passed correctly."""

        for method in [CutoffNonPeriodic, CutoffPeriodic, Ewald, PME]:
            system = self.prmtop1.createSystem(nonbondedMethod=method, 
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
            system = self.prmtop1.createSystem(nonbondedMethod=method, 
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
            system = self.prmtop1.createSystem(removeCMMotion=b)
            forces = system.getForces()
            self.assertEqual(any(isinstance(f, CMMotionRemover) for f in forces), b)

    def test_RigidWaterAndConstraints(self):
        """Test all eight options for the constraints and rigidWater parameters."""

        topology = self.prmtop1.topology
        for constraints_value in [None, HBonds, AllBonds, HAngles]:
            for rigidWater_value in [True, False]:
                system = self.prmtop1.createSystem(constraints=constraints_value, 
                                                   rigidWater=rigidWater_value)
                validateConstraints(self, topology, system, 
                                    constraints_value, rigidWater_value)

    def test_ImplicitSolvent(self):
        """Test the four types of implicit solvents using the implicitSolvent 
        parameter.

        """
        for implicitSolvent_value in [HCT, OBC1, OBC2, GBn]:
            system = self.prmtop2.createSystem(implicitSolvent=implicitSolvent_value)
            forces = system.getForces()
            if implicitSolvent_value in set([HCT, OBC1, GBn]):
                force_type = CustomGBForce
            else:
                force_type = GBSAOBCForce
            
            self.assertTrue(any(isinstance(f, force_type) for f in forces))

    def test_ImplicitSolventParameters(self):
        """Test that parameters are set correctly for the different types of implicit solvent."""
        methodMap = {NoCutoff:NonbondedForce.NoCutoff,
                     CutoffNonPeriodic:NonbondedForce.CutoffNonPeriodic}
        for implicitSolvent_value in [HCT, OBC1, OBC2, GBn]:
            for method in methodMap:
                system = self.prmtop2.createSystem(implicitSolvent=implicitSolvent_value, 
                                    solventDielectric=50.0, soluteDielectric=0.9, nonbondedMethod=method)
                found_matching_solvent_dielectric=False
                found_matching_solute_dielectric=False
                if implicitSolvent_value in set([HCT, OBC1, GBn]):
                    for force in system.getForces():
                        if isinstance(force, CustomGBForce):
                            self.assertEqual(force.getNonbondedMethod(), methodMap[method])
                            for j in range(force.getNumGlobalParameters()):
                                if (force.getGlobalParameterName(j) == 'solventDielectric' and
                                   force.getGlobalParameterDefaultValue(j) == 50.0):
                                    found_matching_solvent_dielectric = True
                                if (force.getGlobalParameterName(j) == 'soluteDielectric' and
                                   force.getGlobalParameterDefaultValue(j) == 0.9):
                                    found_matching_solute_dielectric = True
                        if isinstance(force, NonbondedForce):
                            self.assertEqual(force.getReactionFieldDielectric(), 1.0)
                            self.assertEqual(force.getNonbondedMethod(), methodMap[method])
                    self.assertTrue(found_matching_solvent_dielectric and 
                                    found_matching_solute_dielectric)
                else:
                    for force in system.getForces():
                        if isinstance(force, GBSAOBCForce):
                            self.assertEqual(force.getNonbondedMethod(), methodMap[method])
                            if force.getSolventDielectric() == 50.0:
                                found_matching_solvent_dielectric = True
                            if force.getSoluteDielectric() == 0.9:
                                found_matching_solute_dielectric = True
                        if isinstance(force, NonbondedForce):
                            self.assertEqual(force.getReactionFieldDielectric(), 1.0)
                            self.assertEqual(force.getNonbondedMethod(), methodMap[method])
                    self.assertTrue(found_matching_solvent_dielectric and 
                                    found_matching_solute_dielectric)

    def test_HydrogenMass(self):
        """Test that altering the mass of hydrogens works correctly."""
        
        topology = self.prmtop1.topology
        hydrogenMass = 4*amu
        system1 = self.prmtop1.createSystem()
        system2 = self.prmtop1.createSystem(hydrogenMass=hydrogenMass)
        for atom in topology.atoms():
            if atom.element == elem.hydrogen:
                self.assertNotEqual(hydrogenMass, system1.getParticleMass(atom.index))
                self.assertEqual(hydrogenMass, system2.getParticleMass(atom.index))
        totalMass1 = sum([system1.getParticleMass(i) for i in range(system1.getNumParticles())]).value_in_unit(amu)
        totalMass2 = sum([system2.getParticleMass(i) for i in range(system2.getNumParticles())]).value_in_unit(amu)
        self.assertAlmostEqual(totalMass1, totalMass2)

if __name__ == '__main__':
    unittest.main()

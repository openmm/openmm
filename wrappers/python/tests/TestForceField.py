import unittest
from validateConstraints import *
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem

class TestForceField(unittest.TestCase):
    """Test the ForceField.createSystem() method."""
 
    def setUp(self):
        """Set up the tests by loading the input pdb files and force field 
        xml files.

        """
        # alanine dipeptide with explicit water
        self.pdb1 = PDBFile('systems/alanine-dipeptide-explicit.pdb')
        self.forcefield1 = ForceField('amber99sb.xml', 'tip3p.xml')
        self.topology1 = self.pdb1.topology
        self.topology1.setUnitCellDimensions(Vec3(2, 2, 2))

        # alalnine dipeptide with implicit water
        self.pdb2 = PDBFile('systems/alanine-dipeptide-implicit.pdb')
        self.forcefield2 = ForceField('amber99sb.xml', 'amber99_obc.xml')


    def test_NonbondedMethod(self):
        """Test all five options for the nonbondedMethod parameter."""

        methodMap = {NoCutoff:NonbondedForce.NoCutoff, 
                     CutoffNonPeriodic:NonbondedForce.CutoffNonPeriodic, 
                     CutoffPeriodic:NonbondedForce.CutoffPeriodic, 
                     Ewald:NonbondedForce.Ewald, PME: NonbondedForce.PME}
        for method in methodMap:
            system = self.forcefield1.createSystem(self.pdb1.topology,
                                                  nonbondedMethod=method)
            forces = system.getForces()
            self.assertTrue(any(isinstance(f, NonbondedForce) and 
                                f.getNonbondedMethod()==methodMap[method] 
                                for f in forces))

    def test_DispersionCorrection(self):
        """Test to make sure the nonbondedCutoff parameter is passed correctly."""

        for useDispersionCorrection in [True, False]:
            system = self.forcefield1.createSystem(self.pdb1.topology,
                                                   nonbondedCutoff=2*nanometer, 
                                                   useDispersionCorrection=useDispersionCorrection)

            for force in system.getForces():
                if isinstance(force, NonbondedForce):
                    self.assertEqual(useDispersionCorrection, force.getUseDispersionCorrection())

    def test_Cutoff(self):
        """Test to make sure the nonbondedCutoff parameter is passed correctly."""

        for method in [CutoffNonPeriodic, CutoffPeriodic, Ewald, PME]:
            system = self.forcefield1.createSystem(self.pdb1.topology,
                                                   nonbondedMethod=method, 
                                                   nonbondedCutoff=2*nanometer, 
                                                   constraints=HBonds)
            cutoff_distance = 0.0*nanometer
            cutoff_check = 2.0*nanometer
            for force in system.getForces():
                if isinstance(force, NonbondedForce):
                    cutoff_distance = force.getCutoffDistance()
            self.assertEqual(cutoff_distance, cutoff_check)

    def test_RemoveCMMotion(self):
        """Test both options (True and False) for the removeCMMotion parameter."""
        for b in [True, False]:
            system = self.forcefield1.createSystem(self.pdb1.topology,removeCMMotion=b)
            forces = system.getForces()
            self.assertEqual(any(isinstance(f, CMMotionRemover) for f in forces), b)

    def test_RigidWaterAndConstraints(self):
        """Test all eight options for the constraints and rigidWater parameters."""
        
        topology = self.pdb1.topology
        for constraints_value in [None, HBonds, AllBonds, HAngles]:
            for rigidWater_value in [True, False]:
                system = self.forcefield1.createSystem(topology,
                                                       constraints=constraints_value, 
                                                       rigidWater=rigidWater_value)
                validateConstraints(self, topology, system, 
                                    constraints_value, rigidWater_value)

    def test_ImplicitSolvent(self):
        """Test the four types of implicit solvents using the implicitSolvent 
        parameter.

        """
        
        topology = self.pdb2.topology 
        system = self.forcefield2.createSystem(topology)
        forces = system.getForces()
        self.assertTrue(any(isinstance(f, GBSAOBCForce) for f in forces))

    def test_ImplicitSolventParameters(self):
        """Test that solventDielectric and soluteDielectric are passed correctly 
        for the different types of implicit solvent.

        """

        topology = self.pdb2.topology
        system = self.forcefield2.createSystem(topology, solventDielectric=50.0,
                                               soluteDielectric=0.9)
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
        
        topology = self.pdb1.topology
        hydrogenMass = 4*amu
        system1 = self.forcefield1.createSystem(topology)
        system2 = self.forcefield1.createSystem(topology, hydrogenMass=hydrogenMass)
        for atom in topology.atoms():
            if atom.element == elem.hydrogen:
                self.assertNotEqual(hydrogenMass, system1.getParticleMass(atom.index))
                self.assertEqual(hydrogenMass, system2.getParticleMass(atom.index))
        totalMass1 = sum([system1.getParticleMass(i) for i in range(system1.getNumParticles())]).value_in_unit(amu)
        totalMass2 = sum([system2.getParticleMass(i) for i in range(system2.getNumParticles())]).value_in_unit(amu)
        self.assertAlmostEqual(totalMass1, totalMass2)
    
    def test_Forces(self):
        """Compute forces and compare them to ones generated with a previous version of OpenMM to ensure they haven't changed."""
        
        pdb = PDBFile('systems/lysozyme-implicit.pdb')
        system = self.forcefield2.createSystem(pdb.topology)
        integrator = VerletIntegrator(0.001)
        context = Context(system, integrator)
        context.setPositions(pdb.positions)
        state1 = context.getState(getForces=True)
        state2 = XmlSerializer.deserialize(open('systems/lysozyme-implicit-forces.xml').read())
        numDifferences = 0
        for f1, f2, in zip(state1.getForces().value_in_unit(kilojoules_per_mole/nanometer), state2.getForces().value_in_unit(kilojoules_per_mole/nanometer)):
            diff = norm(f1-f2)
            if diff > 0.1 and diff/norm(f1) > 1e-3:
                numDifferences += 1
        self.assertTrue(numDifferences < system.getNumParticles()/20) # Tolerate occasional differences from numerical error

class AmoebaTestForceField(unittest.TestCase):
    """Test the ForceField.createSystem() method with the AMOEBA forcefield."""
 
    def setUp(self):
        """Set up the tests by loading the input pdb files and force field 
        xml files.

        """

        self.pdb1 = PDBFile('systems/amoeba-ion-in-water.pdb')
        self.forcefield1 = ForceField('amoeba2009.xml')
        self.topology1 = self.pdb1.topology


    def test_NonbondedMethod(self):
        """Test all five options for the nonbondedMethod parameter."""

        methodMap = {NoCutoff:AmoebaMultipoleForce.NoCutoff,
                     PME:AmoebaMultipoleForce.PME}

        for method in methodMap:
            system = self.forcefield1.createSystem(self.pdb1.topology,
                                                  nonbondedMethod=method)
            forces = system.getForces()
            self.assertTrue(any(isinstance(f, AmoebaMultipoleForce) and
                                f.getNonbondedMethod()==methodMap[method]
                                for f in forces))
    def test_Cutoff(self):
        """Test to make sure the nonbondedCutoff parameter is passed correctly."""

        cutoff_distance = 0.7*nanometer
        for method in [NoCutoff, PME]:
            system = self.forcefield1.createSystem(self.pdb1.topology,
                                                   nonbondedMethod=method,
                                                   nonbondedCutoff=cutoff_distance,
                                                   constraints=None)

            for force in system.getForces():
                if isinstance(force, AmoebaVdwForce):
                    self.assertEqual(force.getCutoff(), cutoff_distance)
                if isinstance(force, AmoebaMultipoleForce):
                    self.assertEqual(force.getCutoffDistance(), cutoff_distance)

    def test_DispersionCorrection(self):
        """Test to make sure the nonbondedCutoff parameter is passed correctly."""

        for useDispersionCorrection in [True, False]:
            system = self.forcefield1.createSystem(self.pdb1.topology,
                                                   nonbondedMethod=PME,
                                                   useDispersionCorrection=useDispersionCorrection)

            for force in system.getForces():
                if isinstance(force, AmoebaVdwForce):
                    self.assertEqual(useDispersionCorrection, force.getUseDispersionCorrection())


if __name__ == '__main__':
    unittest.main()


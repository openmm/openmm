import unittest
from validateConstraints import *
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk.openmm.app.gromacstopfile import _defaultGromacsIncludeDir
import simtk.openmm.app.element as elem

GROMACS_INCLUDE = _defaultGromacsIncludeDir()

@unittest.skipIf(not os.path.exists(GROMACS_INCLUDE), 'GROMACS is not installed')
class TestGromacsTopFile(unittest.TestCase):

    """Test the GromacsTopFile.createSystem() method."""

    def setUp(self):
        """Set up the tests by loading the input files."""

        # alanine dipeptide with explicit water
        self.top1 = GromacsTopFile('systems/explicit.top', unitCellDimensions=Vec3(6.223, 6.223, 6.223)*nanometers)
        # alanine dipeptide with implicit water
        self.top2 = GromacsTopFile('systems/implicit.top')

    def test_NonbondedMethod(self):
        """Test all six options for the nonbondedMethod parameter."""

        methodMap = {NoCutoff:NonbondedForce.NoCutoff,
                     CutoffNonPeriodic:NonbondedForce.CutoffNonPeriodic,
                     CutoffPeriodic:NonbondedForce.CutoffPeriodic,
                     Ewald:NonbondedForce.Ewald,
                     PME:NonbondedForce.PME,
                     LJPME:NonbondedForce.LJPME}
        for method in methodMap:
            system = self.top1.createSystem(nonbondedMethod=method)
            forces = system.getForces()
            self.assertTrue(any(isinstance(f, NonbondedForce) and
                                f.getNonbondedMethod()==methodMap[method]
                                for f in forces))

    def test_ff99SBILDN(self):
        """ Test Gromacs topology #define replacement as used in ff99SB-ILDN """
        top = GromacsTopFile('systems/aidilnaaaaa.top')
        gro = GromacsGroFile('systems/aidilnaaaaa.gro')
        system = top.createSystem()
        for force in system.getForces():
            if isinstance(force, PeriodicTorsionForce):
                force.setForceGroup(1)
        context = Context(system, VerletIntegrator(1*femtosecond),
                          Platform.getPlatformByName('Reference'))
        context.setPositions(gro.positions)
        ene = context.getState(getEnergy=True, groups=2**1).getPotentialEnergy()
        self.assertAlmostEqual(ene.value_in_unit(kilojoules_per_mole), 341.6905133582857)

    def test_SMOG(self):
        """ Test to ensure that SMOG models can be run without problems """
        top = GromacsTopFile('systems/2ci2.pdb.top')
        gro = GromacsGroFile('systems/2ci2.pdb.gro')
        system = top.createSystem()

        context = Context(system, VerletIntegrator(1*femtosecond),
                          Platform.getPlatformByName('Reference'))
        context.setPositions(gro.positions)
        ene = context.getState(getEnergy=True).getPotentialEnergy()
        self.assertAlmostEqual(ene.value_in_unit(kilojoules_per_mole), -346.940915296)

    def test_ionic(self):
        """Test simulating an ionic liquid"""
        gro = GromacsGroFile('systems/ionic.gro')
        top = GromacsTopFile('systems/ionic.top', periodicBoxVectors=gro.getPeriodicBoxVectors())
        system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.2)
        for f in system.getForces():
            if isinstance(f, CustomNonbondedForce):
                f.setUseLongRangeCorrection(True)

        context = Context(system, VerletIntegrator(1*femtosecond),
                          Platform.getPlatformByName('Reference'))
        context.setPositions(gro.positions)
        energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        self.assertAlmostEqual(energy, 3135.33, delta=energy*0.005)
        self.assertEqual(1400, system.getNumConstraints())

    def test_Cutoff(self):
        """Test to make sure the nonbondedCutoff parameter is passed correctly."""

        for method in [CutoffNonPeriodic, CutoffPeriodic, Ewald, PME, LJPME]:
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

        for method in [Ewald, PME, LJPME]:
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

    def test_VirtualParticle(self):
        """Test virtual particle works correctly."""

        top = GromacsTopFile('systems/bnz.top')
        gro = GromacsGroFile('systems/bnz.gro')
        system = top.createSystem()

        self.assertEqual(26, system.getNumParticles())
        self.assertEqual(1, len(top._moleculeTypes['BENX'].vsites2))

        context = Context(system, VerletIntegrator(1*femtosecond),
                          Platform.getPlatformByName('Reference'))
        context.setPositions(gro.positions)
        context.computeVirtualSites()
        ene = context.getState(getEnergy=True).getPotentialEnergy()
        # the energy output is from gromacs and it only prints out 6 sig digits.
        self.assertAlmostEqual(ene.value_in_unit(kilojoules_per_mole), 1.88855e+02, places=3)

if __name__ == '__main__':
    unittest.main()


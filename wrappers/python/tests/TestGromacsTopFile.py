import unittest
from validateConstraints import *
from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm.app.gromacstopfile import _defaultGromacsIncludeDir
import openmm.app.element as elem
from numpy.testing import assert_allclose

GROMACS_INCLUDE = _defaultGromacsIncludeDir()

@unittest.skipIf(not os.path.exists(GROMACS_INCLUDE), 'GROMACS is not installed')
class TestGromacsTopFile(unittest.TestCase):

    """Test the GromacsTopFile.createSystem() method."""

    def setUp(self):
        """Set up the tests by loading the input files."""

        # alanine dipeptide with explicit water
        self.top1 = GromacsTopFile('systems/explicit.top', unitCellDimensions=Vec3(6.223, 6.223, 6.223)*nanometers)

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
                          Platform.getPlatform('Reference'))
        context.setPositions(gro.positions)
        ene = context.getState(getEnergy=True, groups=2**1).getPotentialEnergy()
        self.assertAlmostEqual(ene.value_in_unit(kilojoules_per_mole), 341.6905133582857)

    def test_SMOG(self):
        """ Test to ensure that SMOG models can be run without problems """
        top = GromacsTopFile('systems/2ci2.pdb.top')
        gro = GromacsGroFile('systems/2ci2.pdb.gro')
        system = top.createSystem()

        context = Context(system, VerletIntegrator(1*femtosecond),
                          Platform.getPlatform('Reference'))
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
                          Platform.getPlatform('Reference'))
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

    def test_SwitchingFunction(self):
        """Test using a switching function."""
        top = GromacsTopFile('systems/ionic.top')
        for distance in (None, 0.8*nanometers):
            system = top.createSystem(nonbondedMethod=CutoffNonPeriodic, switchDistance=distance)
            for f in system.getForces():
                if isinstance(f, NonbondedForce) or isinstance(f, CustomNonbondedForce):
                    if distance is None:
                        self.assertFalse(f.getUseSwitchingFunction())
                    else:
                        self.assertTrue(f.getUseSwitchingFunction())
                        self.assertEqual(distance, f.getSwitchingDistance())

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

    def test_HydrogenMass(self):
        """Test that altering the mass of hydrogens works correctly."""

        topology = self.top1.topology
        hydrogenMass = 4*amu
        system1 = self.top1.createSystem()
        system2 = self.top1.createSystem(hydrogenMass=hydrogenMass)
        for atom in topology.atoms():
            if atom.element == elem.hydrogen:
                self.assertNotEqual(hydrogenMass, system1.getParticleMass(atom.index))
                if atom.residue.name == 'HOH':
                    self.assertEqual(system1.getParticleMass(atom.index), system2.getParticleMass(atom.index))
                else:
                    self.assertEqual(hydrogenMass, system2.getParticleMass(atom.index))
        totalMass1 = sum([system1.getParticleMass(i) for i in range(system1.getNumParticles())]).value_in_unit(amu)
        totalMass2 = sum([system2.getParticleMass(i) for i in range(system2.getNumParticles())]).value_in_unit(amu)
        self.assertAlmostEqual(totalMass1, totalMass2)

    def test_VirtualParticle(self):
        """Test virtual particle works correctly."""

        top = GromacsTopFile('systems/bnz.top')
        gro = GromacsGroFile('systems/bnz.gro')
        for atom in top.topology.atoms():
            if atom.name.startswith('C'):
                self.assertEqual(elem.carbon, atom.element)
            elif atom.name.startswith('H'):
                self.assertEqual(elem.hydrogen, atom.element)
            else:
                self.assertIsNone(atom.element)
        system = top.createSystem()

        self.assertEqual(26, system.getNumParticles())
        self.assertEqual(1, len(top._moleculeTypes['BENX'].vsites2))

        context = Context(system, VerletIntegrator(1*femtosecond),
                          Platform.getPlatform('Reference'))
        context.setPositions(gro.positions)
        context.computeVirtualSites()
        ene = context.getState(getEnergy=True).getPotentialEnergy()
        # the energy output is from gromacs and it only prints out 6 sig digits.
        self.assertAlmostEqual(ene.value_in_unit(kilojoules_per_mole), 1.88855e+02, places=3)

    def test_Vsite3Func1(self):
        """Test a three particle virtual site."""
        top = GromacsTopFile('systems/tip4pew.top')
        system = top.createSystem()
        self.assertEqual(3, system.getNumConstraints())
        self.assertTrue(system.isVirtualSite(3))
        vs = system.getVirtualSite(3)
        self.assertIsInstance(vs, ThreeParticleAverageSite)
        self.assertEqual(0, vs.getParticle(0))
        self.assertEqual(1, vs.getParticle(1))
        self.assertEqual(2, vs.getParticle(2))
        self.assertAlmostEqual(0.786646558, vs.getWeight(0))
        self.assertAlmostEqual(0.106676721, vs.getWeight(1))
        self.assertAlmostEqual(0.106676721, vs.getWeight(2))

    def test_Vsite3Func4(self):
        """Test a three particle virtual site."""
        top = GromacsTopFile('systems/tip5p.top')
        system = top.createSystem()
        self.assertEqual(3, system.getNumConstraints())
        for i in (3, 4):
            self.assertTrue(system.isVirtualSite(i))
            vs = system.getVirtualSite(i)
            self.assertIsInstance(vs, OutOfPlaneSite)
            self.assertEqual(0, vs.getParticle(0))
            self.assertEqual(1, vs.getParticle(1))
            self.assertEqual(2, vs.getParticle(2))
            self.assertAlmostEqual(-0.344908, vs.getWeight12())
            self.assertAlmostEqual(-0.344908, vs.getWeight13())
            wc = -6.4437903493
            if i == 4:
                wc = -wc
            self.assertAlmostEqual(wc, vs.getWeightCross())

    def test_GROMOS(self):
        """Test a system using the GROMOS 54a7 force field."""

        top = GromacsTopFile('systems/1ppt.top')
        gro = GromacsGroFile('systems/1ppt.gro')
        system = top.createSystem()
        for i, f in enumerate(system.getForces()):
            f.setForceGroup(i)
        context = Context(system, VerletIntegrator(1*femtosecond), Platform.getPlatform('Reference'))
        context.setPositions(gro.positions)
        energy = {}
        for i, f in enumerate(system.getForces()):
            energy[f.getName()] = context.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(kilojoules_per_mole)

        # Compare to energies computed with GROMACS.

        assert_allclose(1.12797e+03, energy['GROMOSBondForce'], rtol=1e-4)
        assert_allclose(5.59066e+02, energy['GROMOSAngleForce'], rtol=1e-4)
        assert_allclose(3.80152e+02, energy['PeriodicTorsionForce'], rtol=1e-4)
        assert_allclose(9.59178e+01, energy['HarmonicTorsionForce'], rtol=1e-4)
        assert_allclose(2.75307e+02, energy['LennardJonesExceptions'], rtol=1e-4)
        assert_allclose(-7.53704e+02, energy['LennardJonesForce'], rtol=1e-4)
        assert_allclose(-6.23055e+03+4.36880e+03, energy['NonbondedForce'], rtol=1e-4)
        total = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        assert_allclose(-1.77020e+02, total, rtol=1e-3)

if __name__ == '__main__':
    unittest.main()


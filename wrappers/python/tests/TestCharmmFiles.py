import unittest
from validateConstraints import *
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem

class TestCharmmFiles(unittest.TestCase):

    """Test the GromacsTopFile.createSystem() method."""
 
    def setUp(self):
        """Set up the tests by loading the input files."""

        # alanine tripeptide; no waters
        self.psf_c = CharmmPsfFile('systems/ala_ala_ala.psf')
        self.psf_x = CharmmPsfFile('systems/ala_ala_ala.xpsf')
        self.psf_v = CharmmPsfFile('systems/ala_ala_ala.vpsf')
        self.params = CharmmParameterSet(
                            'systems/charmm22.rtf', 'systems/charmm22.par')
        self.pdb = PDBFile('systems/ala_ala_ala.pdb')

    def test_NonbondedMethod(self):
        """Test both non-periodic methods for the systems"""

        methodMap = {NoCutoff:NonbondedForce.NoCutoff, 
                     CutoffNonPeriodic:NonbondedForce.CutoffNonPeriodic}
        for top in (self.psf_c, self.psf_x, self.psf_v):
            for method in methodMap:
                system = top.createSystem(self.params, nonbondedMethod=method)
                forces = system.getForces()
                self.assertTrue(any(isinstance(f, NonbondedForce) and 
                                    f.getNonbondedMethod()==methodMap[method] 
                                    for f in forces))

    def test_Cutoff(self):
        """Test to make sure the nonbondedCutoff parameter is passed correctly."""

        for top in (self.psf_c, self.psf_x, self.psf_v):
            for method in [CutoffNonPeriodic]:
                system = top.createSystem(self.params, nonbondedMethod=method,
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
            system = self.psf_c.createSystem(self.params, removeCMMotion=b)
            self.assertEqual(any(isinstance(f, CMMotionRemover) for f in system.getForces()), b)

    def test_ImplicitSolvent(self):
        """Test implicit solvent using the implicitSolvent parameter.

        """
        system = self.psf_v.createSystem(self.params, implicitSolvent=OBC2)
        self.assertTrue(any(isinstance(f, CustomGBForce) for f in system.getForces()))

    def test_ImplicitSolventParameters(self):
        """Test that solventDielectric and soluteDielectric are passed correctly.

        """
        system = self.psf_x.createSystem(self.params, implicitSolvent=GBn,
                                         solventDielectric=50.0, 
                                         soluteDielectric = 0.9)
        found_matching_solvent_dielectric=False
        found_matching_solute_dielectric=False
        for force in system.getForces():
            if isinstance(force, CustomGBForce):
                for i in range(force.getNumGlobalParameters()):
                    if force.getGlobalParameterName(i) == 'solventDielectric':
                        if force.getGlobalParameterDefaultValue(i) == 50.0:
                            found_matching_solvent_dielectric = True
                    elif force.getGlobalParameterName(i) == 'soluteDielectric':
                        if force.getGlobalParameterDefaultValue(i) == 0.9:
                            found_matching_solute_dielectric = True
            if isinstance(force, NonbondedForce):
                self.assertEqual(force.getReactionFieldDielectric(), 1.0)
        self.assertTrue(found_matching_solvent_dielectric and 
                        found_matching_solute_dielectric)

    def test_HydrogenMass(self):
        """Test that altering the mass of hydrogens works correctly."""
        
        topology = self.psf_v.topology
        hydrogenMass = 4*amu
        system1 = self.psf_v.createSystem(self.params)
        system2 = self.psf_v.createSystem(self.params, hydrogenMass=hydrogenMass)
        for atom in topology.atoms():
            if atom.element == elem.hydrogen:
                self.assertNotEqual(hydrogenMass, system1.getParticleMass(atom.index))
                self.assertEqual(hydrogenMass, system2.getParticleMass(atom.index))
        totalMass1 = sum([system1.getParticleMass(i) for i in range(system1.getNumParticles())]).value_in_unit(amu)
        totalMass2 = sum([system2.getParticleMass(i) for i in range(system2.getNumParticles())]).value_in_unit(amu)
        self.assertAlmostEqual(totalMass1, totalMass2)

    def test_NBFIX(self):
        """Tests CHARMM systems with NBFIX Lennard-Jones modifications"""
        import warnings
        warnings.filterwarnings('ignore', category=CharmmPSFWarning)
        psf = CharmmPsfFile('systems/ala3_solv.psf')
        crd = CharmmCrdFile('systems/ala3_solv.crd')
        params = CharmmParameterSet('systems/par_all36_prot.prm',
                                    'systems/toppar_water_ions.str')
        # Box dimensions (found from bounding box)
        psf.setBox(32.7119500*angstroms, 32.9959600*angstroms, 33.0071500*angstroms)

        # Turn off charges so we only test the Lennard-Jones energies
        for a in psf.atom_list:
            a.charge = 0.0

        # Now compute the full energy
        plat = Platform.getPlatformByName('Reference')
        system = psf.createSystem(params, nonbondedMethod=PME,
                                  nonbondedCutoff=8*angstroms)

        con = Context(system, VerletIntegrator(2*femtoseconds), plat)
        con.setPositions(crd.positions)

        state = con.getState(getEnergy=True, enforcePeriodicBox=True)
        ene = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
        self.assertAlmostEqual(ene, 15490.0033559, delta=0.05)

if __name__ == '__main__':
    unittest.main()


from collections import defaultdict
import unittest
import random

from validateModeller import *
from openmm.app import *
from openmm import *
from openmm.unit import *

if sys.version_info >= (3, 0):
    from io import StringIO
else:
    from cStringIO import StringIO

class TestModeller(unittest.TestCase):
    """ Test the Modeller class. """

    def setUp(self):
        # load the alanine dipeptide pdb file
        self.pdb = PDBFile('systems/alanine-dipeptide-explicit.pdb')
        self.topology_start = self.pdb.topology
        self.positions = self.pdb.positions
        self.forcefield = ForceField('amber10.xml', 'tip3p.xml')

        # load the T4-lysozyme-L99A receptor pdb file
        self.pdb2 = PDBFile('systems/lysozyme-implicit.pdb')
        self.topology_start2 = self.pdb2.topology
        self.positions2 = self.pdb2.positions

        # load the metallothionein pdb file
        self.pdb3 =  PDBFile('systems/1T2Y.pdb')
        self.topology_start3 = self.pdb3.topology
        self.positions3 = self.pdb3.positions

    def test_deleteWater(self):
        """ Test the deleteWater() method. """

        # build the chain dictionary
        chain_dict = {0:0}
        # 749 water chains are deleted
        chain_delta = -749

        # Build the residue and atom dictionaries for validate_preserved.
        # Also, count the number of deleted residues and atoms.
        residues_preserved = 0
        residue_delta = 0
        residue_dict = {}
        atoms_preserved = 0
        atom_delta = 0
        atom_dict = {}
        for residue in self.topology_start.residues():
            if residue.name!='HOH' and residue.name!='WAT':
                residue_dict[residue.index] = residues_preserved
                residues_preserved += 1
                for atom in residue.atoms():
                    atom_dict[atom.index] = atoms_preserved
                    atoms_preserved += 1
            else:
                residue_delta -= 1
                for atom in residue.atoms():
                    atom_delta -= 1

        modeller = Modeller(self.topology_start, self.positions)
        modeller.deleteWater()
        topology_after = modeller.getTopology()

        validate_preserved(self, self.topology_start, topology_after,
                           chain_dict, residue_dict, atom_dict)
        validate_deltas(self, self.topology_start, topology_after,
                         chain_delta, residue_delta, atom_delta)

    def test_delete(self):
        """ Test the delete() method. """

        modeller = Modeller(self.topology_start, self.positions)
        topology_before = modeller.getTopology()

        # Create the list of items to be deleted.
        # Start with the first 50 water chains
        chains = [chain for chain in topology_before.chains()]
        toDelete = chains[1:51]

        # Next add water residues 103->152 to the list of items to be deleted
        residues = [residue for residue in topology_before.residues()]
        toDelete.extend(residues[103:153])

        # Finally add water atoms 622->771 to the list of items to be deleted
        atoms = [atom for atom in topology_before.atoms()]
        toDelete.extend(atoms[622:772])

        modeller.delete(toDelete)
        topology_after = modeller.getTopology()

        # build the chain dictionary
        chain_dict = {0:0}
        for i in range(1,51):
            chain_dict[i+50] = i
        for i in range(51,101):
            chain_dict[i+100] = i
        for i in range(101, 600):
            chain_dict[i+150] = i

        # build the residue dictionary
        residue_dict = {}
        for i in range(3):
            residue_dict[i] = i
        for i in range(3,53):
            residue_dict[i+50] = i
        for i in range(53, 103):
            residue_dict[i+100] = i
        for i in range(103, 602):
            residue_dict[i+150] = i

        # build the atom dictionary
        atom_dict = {}
        for i in range(22):
            atom_dict[i] = i
        for i in range(22,172):
            atom_dict[i+150] = i
        for i in range(172,322):
            atom_dict[i+300] = i
        for i in range(322,1819):
            atom_dict[i+450] = i

        validate_preserved(self, topology_before, topology_after, chain_dict, residue_dict, atom_dict)

        chain_delta = -150
        residue_delta = -150
        atom_delta = -450

        validate_deltas(self, topology_before, topology_after, chain_delta, residue_delta, atom_delta)

    def test_add(self):
        """ Test the add() method. """

        # load the methanol-box pdb file
        pdb2 = PDBFile('systems/methanol-box.pdb')
        topology_toAdd = pdb2.topology
        positions_toAdd = pdb2.positions

        modeller = Modeller(self.topology_start, self.positions)
        modeller.deleteWater()
        topology_before = modeller.getTopology()

        modeller.add(topology_toAdd, positions_toAdd)
        topology_after = modeller.getTopology()

        # build the first chain dictionary for the first call of validate_preserved()
        chain_counter = 0
        chain_dict = {}
        for chain in topology_before.chains():
            chain_dict[chain.index] = chain_counter
            chain_counter += 1

        # build the residue and atom dictionaries for the first call of validate_preserved()
        residue_counter = 0
        residue_dict = {}
        atom_counter = 0
        atom_dict = {}
        for residue in topology_before.residues():
            residue_dict[residue.index] = residue_counter
            residue_counter += 1
            for atom in residue.atoms():
                atom_dict[atom.index] = atom_counter
                atom_counter += 1

        # Validate that the items from the before topology are preserved after addition of items.
        validate_preserved(self, topology_before, topology_after, chain_dict, residue_dict, atom_dict)

        # Next, we build another set of dictionaries to validate that the items added are
        # preserved.  Also, we calculate the number of chains, residues, and atoms added.

        # build the chain dictionary
        chain_delta = 0
        chain_dict = {}
        for chain in topology_toAdd.chains():
            chain_dict[chain.index] = chain_counter
            chain_counter += 1
            chain_delta += 1

        # build the residue and atom dictionaries for the second call of validate_preserved
        residue_delta = 0
        residue_dict = {}
        atom_delta = 0
        atom_dict = {}
        for residue in topology_toAdd.residues():
            residue_dict[residue.index] = residue_counter
            residue_counter += 1
            residue_delta += 1
            for atom in residue.atoms():
                atom_dict[atom.index] = atom_counter
                atom_counter += 1
                atom_delta += 1

        # validate that the items in the added topology are preserved
        validate_preserved(self, topology_toAdd, topology_after, chain_dict, residue_dict, atom_dict)
        # validate that the final topology has the correct number of items
        validate_deltas(self, topology_before, topology_after, chain_delta, residue_delta, atom_delta)

    def test_convertWater(self):
        """ Test the convertWater() method. """

        for model in ['tip3p', 'spce', 'tip4pew', 'tip5p']:
            if model == 'tip5p':
                firstmodel = 'tip4pew'
            else:
                firstmodel = 'tip5p'

            modeller = Modeller(self.topology_start, self.positions)
            modeller.convertWater(model=firstmodel)
            modeller.convertWater(model=model)
            topology_after = modeller.getTopology()

            for residue in topology_after.residues():
                if residue.name == "HOH":
                    oatom = [atom for atom in residue.atoms() if atom.element == element.oxygen]
                    hatoms = [atom for atom in residue.atoms() if atom.element == element.hydrogen]
                    matoms = [atom for atom in residue.atoms() if atom.name == 'M']
                    m1atoms = [atom for atom in residue.atoms() if atom.name == 'M1']
                    m2atoms = [atom for atom in residue.atoms() if atom.name == 'M2']
                    self.assertTrue(len(oatom)==1 and len(hatoms)==2)
                    if model=='tip3p' or model=='spce':
                        self.assertTrue(len(matoms)==0 and len(m1atoms)==0 and len(m2atoms)==0)
                    elif model=='tip4pew':
                        self.assertTrue(len(matoms)==1 and len(m1atoms)==0 and len(m2atoms)==0)
                    elif model=='tip5p':
                        self.assertTrue(len(matoms)==0 and len(m1atoms)==1 and len(m2atoms)==1)

                    # build the chain dictionary for validate_preserved
                    chain_counter = 0
                    chain_dict = {}
                    chain_delta = 0
                    for chain in self.topology_start.chains():
                        chain_dict[chain.index] = chain_counter
                        chain_counter += 1

                    # build the residue and atom dictionaries for validate_preserved
                    residue_counter = 0
                    residue_dict = {}
                    residue_delta = 0
                    atom_counter = 0
                    atom_dict = {}
                    atom_delta = 0
                    for residue in self.topology_start.residues():
                        residue_dict[residue.index] = residue_counter
                        residue_counter += 1
                        for atom in residue.atoms():
                            atom_dict[atom.index] = atom_counter
                            atom_counter += 1
                        if residue.name == 'HOH' and model == 'tip4pew':
                            atom_counter += 1
                            atom_delta += 1
                        if residue.name == 'HOH' and model == 'tip5p':
                            atom_counter += 2
                            atom_delta += 2

            validate_preserved(self, self.topology_start, topology_after,
                               chain_dict, residue_dict, atom_dict)

            validate_deltas(self, self.topology_start, topology_after,
                            chain_delta, residue_delta, atom_delta)

    def test_addSolventWaterModels(self):
        """ Test all addSolvent() method with all possible water models. """

        topology_start = self.pdb.topology
        topology_start.setUnitCellDimensions(Vec3(3.5, 3.5, 3.5)*nanometers)
        for model in ['tip3p', 'spce', 'tip4pew', 'tip5p', 'swm4ndp']:
            forcefield = ForceField('amber10.xml', model + '.xml')
            modeller = Modeller(topology_start, self.positions)
            # delete water to get the "before" topology
            modeller.deleteWater()
            topology_before = modeller.getTopology()
            # add the solvent to get the "after" topology
            modeller.addSolvent(forcefield, model=model)
            topology_after = modeller.getTopology()

            # First, check that everything that was there before has been preserved.

            # build the chain dictionary for validate_preserved
            chain_counter = 0
            chain_dict = {0:0}
            for chain in topology_before.chains():
                chain_dict[chain.index] = chain_counter
                chain_counter += 1

            # build the residue and atom dictionaries for validate_preserved
            residue_counter = 0
            residue_dict = {}
            atom_counter = 0
            atom_dict = {}
            for residue in topology_before.residues():
                residue_dict[residue.index] = residue_counter
                residue_counter += 1
                for atom in residue.atoms():
                    atom_dict[atom.index] = atom_counter
                    atom_counter += 1

            # validate that the items in the before topology remain after solvent is added
            validate_preserved(self, topology_before, topology_after, chain_dict, residue_dict, atom_dict)

            # Make sure water that was added was the correct model
            for residue in topology_after.residues():
                if residue.name == 'HOH':
                    oatom = [atom for atom in residue.atoms() if atom.element == element.oxygen]
                    hatoms = [atom for atom in residue.atoms() if atom.element == element.hydrogen]
                    matoms = [atom for atom in residue.atoms() if atom.name == 'M']
                    m1atoms = [atom for atom in residue.atoms() if atom.name == 'M1']
                    m2atoms = [atom for atom in residue.atoms() if atom.name == 'M2']
                    self.assertTrue(len(oatom)==1 and len(hatoms)==2)
                    if model=='tip3p' or model=='spce':
                        self.assertTrue(len(matoms)==0 and len(m1atoms)==0 and len(m2atoms)==0)
                    elif model=='tip4pew':
                        self.assertTrue(len(matoms)==1 and len(m1atoms)==0 and len(m2atoms)==0)
                    elif model=='tip5p':
                        self.assertTrue(len(matoms)==0 and len(m1atoms)==1 and len(m2atoms)==1)

    def test_addSolventPeriodicBox(self):
        """ Test the addSolvent() method; test that the five ways of passing in the periodic box all work. """

        # First way of passing in periodic box vectors:  set it in the original topology.
        topology_start = self.pdb.topology
        topology_start.setUnitCellDimensions(Vec3(3.5, 4.5, 5.5)*nanometers)
        modeller = Modeller(topology_start, self.positions)
        modeller.deleteWater()
        modeller.addSolvent(self.forcefield)
        topology_after = modeller.getTopology()
        dim3 = topology_after.getPeriodicBoxVectors()

        self.assertVecAlmostEqual(dim3[0]/nanometers, Vec3(3.5, 0, 0))
        self.assertVecAlmostEqual(dim3[1]/nanometers, Vec3(0, 4.5, 0))
        self.assertVecAlmostEqual(dim3[2]/nanometers, Vec3(0, 0, 5.5))

        # Second way of passing in the periodic box vectors: with the boxSize parameter to addSolvent()
        modeller = Modeller(topology_start, self.positions)
        modeller.deleteWater()
        modeller.addSolvent(self.forcefield, boxSize = Vec3(3.6, 4.6, 5.6)*nanometers)
        topology_after = modeller.getTopology()
        dim3 = topology_after.getPeriodicBoxVectors()

        self.assertVecAlmostEqual(dim3[0]/nanometers, Vec3(3.6, 0, 0))
        self.assertVecAlmostEqual(dim3[1]/nanometers, Vec3(0, 4.6, 0))
        self.assertVecAlmostEqual(dim3[2]/nanometers, Vec3(0, 0, 5.6))

        # Third way of passing in the periodic box vectors: with the boxVectors parameter to addSolvent()
        modeller = Modeller(topology_start, self.positions)
        modeller.deleteWater()
        modeller.addSolvent(self.forcefield, boxVectors = (Vec3(3.4, 0, 0), Vec3(0.5, 4.4, 0), Vec3(-1.0, -1.5, 5.4))*nanometers)
        topology_after = modeller.getTopology()
        dim3 = topology_after.getPeriodicBoxVectors()

        self.assertVecAlmostEqual(dim3[0]/nanometers, Vec3(3.4, 0, 0))
        self.assertVecAlmostEqual(dim3[1]/nanometers, Vec3(0.5, 4.4, 0))
        self.assertVecAlmostEqual(dim3[2]/nanometers, Vec3(-1.0, -1.5, 5.4))

        # Fourth way of passing in the periodic box vectors: pass a 'padding' value to addSolvent()
        modeller = Modeller(topology_start, self.positions)
        modeller.deleteWater()
        modeller.addSolvent(self.forcefield, padding = 0.9*nanometers)
        topology_after = modeller.getTopology()
        dim3 = topology_after.getPeriodicBoxVectors()

        self.assertVecAlmostEqual(dim3[0]/nanometers, Vec3(1.824363, 0, 0))
        self.assertVecAlmostEqual(dim3[1]/nanometers, Vec3(0, 1.824363, 0))
        self.assertVecAlmostEqual(dim3[2]/nanometers, Vec3(0, 0, 1.824363))

        # Fifth way: specify a number of molecules to add instead of a box size
        modeller = Modeller(topology_start, self.positions)
        modeller.deleteWater()
        numInitial = len(list(modeller.topology.residues()))
        modeller.addSolvent(self.forcefield, numAdded=1000)
        self.assertEqual(numInitial+1000, len(list(modeller.topology.residues())))

    def test_addSolventBoxShape(self):
        """Test the addSolvent() method; test the different box shapes."""
        modeller = Modeller(self.pdb.topology, self.positions)
        modeller.deleteWater()
        modeller.addSolvent(self.forcefield, padding=1.0*nanometers, boxShape='cube')
        cubeVectors = modeller.getTopology().getPeriodicBoxVectors()
        modeller = Modeller(self.pdb.topology, self.positions)
        modeller.deleteWater()
        modeller.addSolvent(self.forcefield, padding=1.0*nanometers, boxShape='dodecahedron')
        dodecVectors = modeller.getTopology().getPeriodicBoxVectors()
        modeller = Modeller(self.pdb.topology, self.positions)
        modeller.deleteWater()
        modeller.addSolvent(self.forcefield, padding=1.0*nanometers, boxShape='octahedron')
        octVectors = modeller.getTopology().getPeriodicBoxVectors()
        cubeVolume = cubeVectors[0][0]*cubeVectors[1][1]*cubeVectors[2][2]/(nanometers**3)
        dodecVolume = dodecVectors[0][0]*dodecVectors[1][1]*dodecVectors[2][2]/(nanometers**3)
        octVolume = octVectors[0][0]*octVectors[1][1]*octVectors[2][2]/(nanometers**3)
        self.assertAlmostEqual(0.707, dodecVolume/cubeVolume, places=3)
        self.assertAlmostEqual(0.770, octVolume/cubeVolume, places=3)

    def test_addSolventNeutralSolvent(self):
        """ Test the addSolvent() method; test adding ions to neutral solvent. """

        topology_start = self.pdb.topology
        topology_start.setUnitCellDimensions(Vec3(3.5, 3.5, 3.5)*nanometers)
        modeller = Modeller(topology_start, self.positions)
        modeller.deleteWater()
        modeller.addSolvent(self.forcefield, ionicStrength = 2.0*molar)
        topology_after = modeller.getTopology()

        water_count=0
        sodium_count=0
        chlorine_count=0
        for residue in topology_after.residues():
            if residue.name=='HOH':
                water_count += 1
            elif residue.name=='NA':
                sodium_count += 1
            elif residue.name=='CL':
                chlorine_count += 1

        total_added = water_count+sodium_count+chlorine_count
        self.assertEqual(total_added, 1364)
        expected_ion_fraction = 2.0*molar/(55.4*molar)
        expected_ions = math.floor(total_added*expected_ion_fraction+0.5)
        self.assertEqual(sodium_count, expected_ions)
        self.assertEqual(chlorine_count, expected_ions)

    def test_addSolventNegativeSolvent(self):
        """ Test the addSolvent() method; test adding ions to a negatively charged solvent. """

        topology_start = self.pdb.topology
        topology_start.setUnitCellDimensions(Vec3(3.5, 3.5, 3.5)*nanometers)

        for neutralize in (True, False):
            # set up modeller with no solvent
            modeller = Modeller(topology_start, self.positions)
            modeller.deleteWater()

            # add 5 Cl- ions to the original topology
            topology_toAdd = Topology()
            newChain = topology_toAdd.addChain()
            for i in range(5):
                topology_toAdd.addResidue('CL',  newChain)
            residues = [residue for residue in topology_toAdd.residues()]
            for i in range(5):
                topology_toAdd.addAtom('Cl',Element.getBySymbol('Cl'), residues[i])
            positions_toAdd = [Vec3(1.0,1.2,1.5), Vec3(1.7,1.0,1.4), Vec3(1.5,2.0,1.0),
                               Vec3(2.0,2.0,2.0), Vec3(2.0,1.5,1.0)]*nanometers
            modeller.add(topology_toAdd, positions_toAdd)
            modeller.addSolvent(self.forcefield, ionicStrength=1.0*molar, neutralize=neutralize)
            topology_after = modeller.getTopology()

            water_count = 0
            sodium_count = 0
            chlorine_count = 0
            for residue in topology_after.residues():
                if residue.name=='HOH':
                    water_count += 1
                elif residue.name=='NA':
                    sodium_count += 1
                elif residue.name=='CL':
                    chlorine_count += 1

            total_water_ions = water_count+sodium_count+chlorine_count
            expected_ion_fraction = 1.0*molar/(55.4*molar)
            expected_chlorine = math.floor((total_water_ions-10)*expected_ion_fraction+0.5)+5
            expected_sodium = expected_chlorine if neutralize else expected_chlorine-5
            self.assertEqual(sodium_count, expected_sodium)
            self.assertEqual(chlorine_count, expected_chlorine)

    def test_addSolventPositiveSolvent(self):
        """ Test the addSolvent() method; test adding ions to a positively charged solvent. """

        topology_start = self.pdb.topology
        topology_start.setUnitCellDimensions(Vec3(3.5, 3.5, 3.5)*nanometers)

        for neutralize in (True, False):
            # set up modeller with no solvent
            modeller = Modeller(topology_start, self.positions)
            modeller.deleteWater()

            # add 5 Na+ ions to the original topology
            topology_toAdd = Topology()
            newChain = topology_toAdd.addChain()
            for i in range(5):
                topology_toAdd.addResidue('NA', newChain)
            residues = [residue for residue in topology_toAdd.residues()]
            for i in range(5):
                 topology_toAdd.addAtom('Na',Element.getBySymbol('Na'), residues[i])
            positions_toAdd = [Vec3(1.0,1.2,1.5), Vec3(1.7,1.0,1.4), Vec3(1.5,2.0,1.0),
                               Vec3(2.0,2.0,2.0), Vec3(2.0,1.5,1.0)]*nanometers

            # positions_toAdd doesn't need to change
            modeller.add(topology_toAdd, positions_toAdd)
            modeller.addSolvent(self.forcefield, ionicStrength=1.0*molar, neutralize=neutralize)
            topology_after = modeller.getTopology()

            water_count = 0
            sodium_count = 0
            chlorine_count = 0
            for residue in topology_after.residues():
                if residue.name=='HOH':
                    water_count += 1
                elif residue.name=='NA':
                    sodium_count += 1
                elif residue.name=='CL':
                    chlorine_count += 1

            total_water_ions = water_count+sodium_count+chlorine_count
            expected_ion_fraction = 1.0*molar/(55.4*molar)
            expected_sodium = math.floor((total_water_ions-10)*expected_ion_fraction+0.5)+5
            expected_chlorine = expected_sodium if neutralize else expected_sodium-5
            self.assertEqual(sodium_count, expected_sodium)
            self.assertEqual(chlorine_count, expected_chlorine)

    def test_addSolventIons(self):
        """ Test the addSolvent() method with all possible choices for positive and negative ions. """

        topology_start = self.pdb.topology
        topology_start.setUnitCellDimensions(Vec3(3.5, 3.5, 3.5)*nanometers)

        # set up modeller with no solvent
        modeller = Modeller(topology_start, self.positions)
        modeller.deleteWater()
        topology_nowater = modeller.getTopology()
        positions_nowater = modeller.getPositions()

        expected_ion_fraction = 1.0*molar/(55.4*molar)

        for positiveIon in ['Cs+', 'K+', 'Li+', 'Na+', 'Rb+']:
            ionName = positiveIon[:-1].upper()
            modeller = Modeller(topology_nowater, positions_nowater)
            modeller.addSolvent(self.forcefield, positiveIon=positiveIon, ionicStrength=1.0*molar)
            topology_after = modeller.getTopology()

            water_count = 0
            positive_ion_count = 0
            chlorine_count = 0
            for residue in topology_after.residues():
                if residue.name=='HOH':
                    water_count += 1
                elif residue.name==ionName:
                    positive_ion_count += 1
                elif residue.name=='CL':
                    chlorine_count += 1

            total_added = water_count+positive_ion_count+chlorine_count
            self.assertEqual(total_added, 1364)
            expected_ions = math.floor(total_added*expected_ion_fraction+0.5)
            self.assertEqual(positive_ion_count, expected_ions)
            self.assertEqual(chlorine_count, expected_ions)

        for negativeIon in ['Cl-', 'Br-', 'F-', 'I-']:
            ionName = negativeIon[:-1].upper()
            modeller = Modeller(topology_nowater, positions_nowater)
            modeller.addSolvent(self.forcefield, negativeIon=negativeIon, ionicStrength=1.0*molar)

            topology_after = modeller.getTopology()

            water_count = 0
            sodium_count = 0
            negative_ion_count = 0
            for residue in topology_after.residues():
                if residue.name=='HOH':
                    water_count += 1
                elif residue.name=='NA':
                    sodium_count += 1
                elif residue.name==ionName:
                    negative_ion_count += 1

            total_added = water_count+sodium_count+negative_ion_count
            self.assertEqual(total_added, 1364)
            expected_ions = math.floor(total_added*expected_ion_fraction+0.5)
            self.assertEqual(positive_ion_count, expected_ions)
            self.assertEqual(chlorine_count, expected_ions)

    def test_addHydrogensPdb2(self):
        """ Test the addHydrogens() method on the T4-lysozyme-L99A pdb file. """

        # build the Modeller
        topology_start = self.topology_start2
        positions = self.positions2
        modeller = Modeller(topology_start, positions)

        # remove hydrogens from the topology
        toDelete = [atom for atom in topology_start.atoms() if atom.element==Element.getBySymbol('H')]
        modeller.delete(toDelete)

        # Create a variants list to force the one histidine to be of the right variation.
        residues = [residue for residue in topology_start.residues()]
        variants = [None]*len(residues)
        # For this protein, when you add hydrogens, the hydrogen is added to the delta nitrogen.
        # By setting variants[30] to 'HIE', we force the hydrogen onto the epsilon nitrogen, so
        # that it will match the topology in topology_start.
        variants[30] = 'HIE'

        # add the hydrogens back
        modeller.addHydrogens(self.forcefield, variants=variants)
        topology_after = modeller.getTopology()

        validate_equivalence(self, topology_start, topology_after)

    def test_addHydrogensPdb3(self):
        """ Test the addHydrogens() method on the metallothionein pdb file. """

        # build the Modeller
        topology_start = self.topology_start3
        positions = self.positions3
        modeller = Modeller(topology_start, positions)

        # remove hydrogens from the topology
        toDelete = [atom for atom in topology_start.atoms() if atom.element==Element.getBySymbol('H')]
        modeller.delete(toDelete)

        # add the hydrogens back
        modeller.addHydrogens(self.forcefield)
        topology_after = modeller.getTopology()

        validate_equivalence(self, topology_start, topology_after)

    def test_addHydrogensPdb3_keepPositions(self):
        """ Test addHydrogens() does not change existing Hs positions """

        # build the Modeller
        topology_start = self.topology_start3
        positions = self.positions3.value_in_unit(nanometers)
        modeller = Modeller(topology_start, positions)

        # Record original hydrogen positions
        oriH = [atom.index for atom in modeller.topology.atoms() if atom.element == element.hydrogen]
        oriH_pos = [positions[i] for i in oriH]

        # Remove hydrogens from last residue
        res_list = list(topology_start.residues())
        toDelete = [atom for atom in res_list[-1].atoms() if atom.element == element.hydrogen]
        modeller.delete(toDelete)

        n_deleted = len(toDelete)

        # Add hydrogen atoms back.
        modeller.addHydrogens(self.forcefield)
        topology_after = modeller.getTopology()

        # Fetch 'new' positions
        new_positions = modeller.positions.value_in_unit(nanometers)
        newH = [atom.index for atom in topology_after.atoms() if atom.element == element.hydrogen]
        newH_pos = [new_positions[i] for i in newH]

        # Did we add all Hs back in correctly?
        self.assertEqual(len(newH), len(oriH))

        # Are the old ones at the same position?
        # Negative control
        oriH_fixed = oriH_pos[:-1*n_deleted]
        newH_fixed = newH_pos[:-1*n_deleted]
        xyz_diff = any([norm(o-n) > 1e-6 for o, n in zip(oriH_fixed, newH_fixed)])
        self.assertEqual(xyz_diff, False)

        # Were the new ones optimized?
        # Positive control
        oriH_added = oriH_pos[-1*n_deleted:]
        newH_added = newH_pos[-1*n_deleted:]
        xyz_diff = all([norm(o-n) > 1e-6 for o, n in zip(oriH_added, newH_added)])
        self.assertEqual(xyz_diff, True)

    def test_addHydrogensASH(self):
        """ Test of addHydrogens() in which we force ASH to be a variant using the variants parameter. """

        # use the T4-lysozyme-L99A pdb file
        topology_start = self.topology_start2
        positions = self.positions2

        # build the Modeller
        modeller = Modeller(topology_start, positions)

        # remove hydrogens from the topology
        toDelete = [atom for atom in topology_start.atoms() if atom.element==Element.getBySymbol('H')]
        modeller.delete(toDelete)

        # Create a variants list to force the one histidine to be of the right variation.
        residues = [residue for residue in topology_start.residues()]
        variants = [None]*len(residues)
        # For this protein, when you add hydrogens, the hydrogen is added to the delta nitrogen.
        # By setting variants[30] to 'HIE', we force the hydrogen onto the epsilon nitrogen, so
        # that it will match the topology in topology_start.
        variants[30] = 'HIE'

        ASP_residue_list = [9,19,46,60,69,71,88,91,126,158]
        for residue_index in ASP_residue_list:
            variants[residue_index] = 'ASH'

        # add the hydrogens back, using the variants list we just built
        modeller.addHydrogens(self.forcefield, variants=variants)
        topology_ASH = modeller.getTopology()

        # There should be extra hydrogens on the ASP residues.  Assert that they exist,
        # then we delete them and validate that the topology matches what we started with.
        index_list_ASH = [176, 357, 761, 976, 1121, 1150, 1430, 1473, 2028, 2556]
        atoms = [atom for atom in topology_ASH.atoms()]
        toDelete2 = []
        for index in index_list_ASH:
            self.assertTrue(atoms[index].element.symbol=='H')
            toDelete2.append(atoms[index])
        modeller.delete(toDelete2)
        topology_ASP = modeller.getTopology()

        validate_equivalence(self, topology_ASP, topology_start)

    def test_addHydrogensCYX(self):
        """ Test of addHydrogens() in which we force CYX to be a variant using the variants parameter. """

        # use the metallothionein pdb file
        topology_start = self.topology_start3
        positions = self.positions3

        # build the Modeller
        modeller = Modeller(topology_start, positions)

        # remove hydrogens from the topology
        toDelete = [atom for atom in topology_start.atoms() if atom.element==Element.getBySymbol('H')]
        modeller.delete(toDelete)

        # Create a variants list to force the cysteins to be of the CYX variety.
        residues = [residue for residue in topology_start.residues()]
        variants = [None]*len(residues)
        CYS_residues = [2,4,10,12,16,18,21]
        for index in CYS_residues:
             variants[index] = 'CYX'

        # add the hydrogens
        modeller.addHydrogens(self.forcefield, variants=variants)
        topology_CYX = modeller.getTopology()

        # create a second modeller that we will attempt to match with topology_CYX
        modeller2 = Modeller(topology_start, positions)
        topology2 = modeller2.getTopology()

        # There should be extra hydrogens on the CYS residues.  Assert that they exist
        # on modeller2, then delete them and validate that the topologies match.

       # These are the indices of the hydrogens to delete from CYS to make CYX.
        index_list_CYS = [31, 49, 110, 135, 171, 193, 229]
        atoms = [atom for atom in topology2.atoms()]
        toDelete2 = []
        for index in index_list_CYS:
            self.assertTrue(atoms[index].element.symbol=='H')
            toDelete2.append(atoms[index])
        modeller2.delete(toDelete2)
        topology_after = modeller2.getTopology()

        validate_equivalence(self, topology_CYX, topology_after)

    def test_addHydrogensGLH(self):
        """ Test of addHydrogens() in which we force GLH to be a variant using the variants parameter. """

        # use the T4-lysozyme-L99A pdb file
        topology_start = self.topology_start2
        positions = self.positions2

        # build the Modeller
        modeller = Modeller(topology_start, positions)

        # remove hydrogens from the topology
        toDelete = [atom for atom in topology_start.atoms() if atom.element==Element.getBySymbol('H')]
        modeller.delete(toDelete)

        # Create a variants list to force the one histidine to be of the right variation.
        residues = [residue for residue in topology_start.residues()]
        variants = [None]*len(residues)
        # For this protein, when you add hydrogens, the hydrogen is added to the delta nitrogen.
        # By setting variants[30] to 'HIE', we force the hydrogen onto the epsilon nitrogen, so
        # that it will match the topology in topology_start.
        variants[30] = 'HIE'

        GLU_residue_list = [4,10,21,44,61,63,107,127]
        for residue_index in GLU_residue_list:
            variants[residue_index] = 'GLH'

        # add the hydrogens back, this time with the GLH variant in place of GLU
        modeller.addHydrogens(self.forcefield, variants=variants)
        topology_GLH = modeller.getTopology()

        # There should be extra hydrogens on the GLU residues.  Assert that they exist,
        # then we delete them and validate that the topology matches what we started with.
        index_list_GLH = [85, 192, 387, 731, 992, 1018, 1718, 2042]
        atoms = [atom for atom in topology_GLH.atoms()]
        toDelete2 = []
        for index in index_list_GLH:
            self.assertTrue(atoms[index].element.symbol=='H')
            toDelete2.append(atoms[index])
        modeller.delete(toDelete2)
        topology_GLU = modeller.getTopology()

        validate_equivalence(self, topology_GLU, topology_start)

    def test_addHydrogensLYN(self):
        """ Test of addHydrogens() in which we force LYN to be a variant using the variants parameter. """

        # use the T4-lysozyme-L99A pdb file
        topology_start = self.topology_start2
        positions = self.positions2

        # build the Modeller
        modeller = Modeller(topology_start, positions)

        # remove hydrogens from topology
        toDelete = [atom for atom in topology_start.atoms() if atom.element==Element.getBySymbol('H')]
        modeller.delete(toDelete)

        # Create a variants list to force the one histidine to be of the right variation.
        residues = [residue for residue in topology_start.residues()]
        variants = [None]*len(residues)
        # For this protein, when you add hydrogens, the hydrogen is added to the delta nitrogen.
        # By setting variants[30] to 'HIE', we force the hydrogen onto the epsilon nitrogen, so
        # that it will match the topology in topology_start.
        variants[30] = 'HIE'

        # Here we add the residues in which LYS is present to the variant list.  The final
        # LYS residue, 161, is not on the list because Amber force fields do not have an
        # entry for a terminal LYN residue.
        residue_list_LYS = [15,18,34,42,47,59,64,82,84,123,134,146]
        for residue_index in residue_list_LYS:
            variants[residue_index] = 'LYN'

        # add the hydrogens back, using the variants list we just built
        modeller.addHydrogens(self.forcefield, variants=variants)

        topology_LYN = modeller.getTopology()

        # create a second modeller that we will attempt to match with topology_LYN
        modeller2 = Modeller(topology_start, positions)

        # There should be extra hydrogens on the LYS residues.  Assert that they exist
        # on modeller2, then delete them and validate that the topologies match.

        # These are the indices of the hydrogens to delete from LYN to make LYS.
        index_list_LYN = [281,343,590,701,780,960,1034,1319,1360,1959,2135,2344]
        atoms = [atom for atom in topology_start.atoms()]
        toDelete2 = []
        for index in index_list_LYN:
             self.assertTrue(atoms[index].element.symbol=='H')
             toDelete2.append(atoms[index])
        modeller2.delete(toDelete2)
        topology_after = modeller2.getTopology()

        validate_equivalence(self, topology_LYN, topology_after)

    def test_addHydrogenspH4(self):
        """ Test of addHydrogens() with pH=4. """

        # use the T4-lysozyme-L99A pdb file
        topology_start = self.topology_start2
        positions = self.positions2

        # build the Modeller
        modeller = Modeller(topology_start, positions)

        # remove hydrogens from the topology
        toDelete = [atom for atom in topology_start.atoms() if atom.element==Element.getBySymbol('H')]
        modeller.delete(toDelete)

        # Create a variants list to force the one histidine to be of the right variation.
        residues = [residue for residue in topology_start.residues()]
        variants = [None]*len(residues)
        # For this protein, when you add hydrogens, the hydrogen is added to the delta nitrogen.
        # By setting variants[30] to 'HIE', we force the hydrogen onto the epsilon nitrogen, so
        # that it will match the topology in topology_start.
        variants[30] = 'HIE'

        # add the hydrogens back, this time at a lower pH
        modeller.addHydrogens(self.forcefield, variants=variants, pH=4.0)

        topology_ASH_GLH = modeller.getTopology()

        # There should be extra hydrogens on the ASP and GLU residues.  Assert that they exist,
        # then we delete them and validate that the topology matches what we started with.
        index_list_ASH = [177, 359, 765, 980, 1127, 1156, 1436, 1479, 2035, 2564]
        index_list_GLH = [85, 193, 389, 733, 996, 1022, 1726, 2051]
        atoms = [atom for atom in topology_ASH_GLH.atoms()]
        toDelete2 = []
        for index in index_list_ASH:
            self.assertTrue(atoms[index].element.symbol=='H')
            toDelete2.append(atoms[index])
        for index in index_list_GLH:
            self.assertTrue(atoms[index].element.symbol=='H')
            toDelete2.append(atoms[index])
        modeller.delete(toDelete2)
        topology_ASP_GLU = modeller.getTopology()

        validate_equivalence(self, topology_ASP_GLU, topology_start)

    def test_addHydrogenspH9(self):
        """ Test of addHydrogens() with pH=9. """

        # use the metallothionein pdb file
        topology_start = self.topology_start3
        positions = self.positions3

        # build the Modeller
        modeller = Modeller(topology_start, positions)

        # remove hydrogens from topology
        toDelete = [atom for atom in topology_start.atoms() if atom.element==Element.getBySymbol('H')]
        modeller.delete(toDelete)

        # add hydrogens with pH=9, so that the variation CYX will be chosen
        modeller.addHydrogens(self.forcefield, pH=9.0)
        topology_CYX = modeller.getTopology()

        # create a second modeller that we will attempt to match with topology_CYX
        modeller2 = Modeller(topology_start, positions)
        topology2 = modeller2.getTopology()

        # There should be extra hydrogens on the CYS residues.  Assert that they exist
        # on modeller2, then delete them and validate that the topologies match.

        # These are the indices of the hydrogens to delete from CYS to make CYX.
        index_list_CYS = [31, 49, 110, 135, 171, 193, 229]
        atoms = [atom for atom in topology2.atoms()]
        toDelete2 = []
        for index in index_list_CYS:
            self.assertTrue(atoms[index].element.symbol=='H')
            toDelete2.append(atoms[index])
        modeller2.delete(toDelete2)
        topology_after = modeller2.getTopology()

        validate_equivalence(self, topology_CYX, topology_after)

    def test_addHydrogenspH11(self):
        """ Test of addHydrogens() with pH=11. """

        # use the T4-lysozyme-L99A pdb file
        topology_start = self.topology_start2
        positions = self.positions2

        # build the Modeller
        modeller = Modeller(topology_start, positions)

        # remove hydrogens from topology
        toDelete = [atom for atom in topology_start.atoms() if atom.element==Element.getBySymbol('H')]
        modeller.delete(toDelete)

        # Create a variants list to force the one histidine to be of the right variation.
        residues = [residue for residue in topology_start.residues()]
        variants = [None]*len(residues)
        # For this protein, when you add hydrogens, the hydrogen is added to the delta nitrogen.
        # By setting variants[30] to 'HIE', we force the hydrogen onto the epsilon nitrogen, so
        # that it will match the topology in topology_start.
        variants[30] = 'HIE'

        # The Amber force fields do not have an entry for terminal LYN residues, so we need to
        # force residue 161 to be the LYS variant.
        variants[161] = 'LYS'

        # add the hydrogens back at pH = 11
        modeller.addHydrogens(self.forcefield, variants=variants, pH=11.0)
        topology_LYN = modeller.getTopology()

        # create a second modeller that we will attempt to match with topology_LYN
        modeller2 = Modeller(topology_start, positions)

        # There should be extra hydrogens on the LYS residues.  Assert that they exist
        # on modeller2, then delete them and validate that the topologies match.
        index_list_LYN = [281,343,590,701,780,960,1034,1319,1360,1959,2135,2344]
        atoms = [atom for atom in topology_start.atoms()]
        toDelete2 = []
        for index in index_list_LYN:
            self.assertTrue(atoms[index].element.symbol=='H')
            toDelete2.append(atoms[index])

        modeller2.delete(toDelete2)
        topology_after = modeller2.getTopology()

        validate_equivalence(self, topology_LYN, topology_after)

    def test_addHydrogensGlycam(self):
        """Test adding hydrogens for GLYCAM."""
        pdb = PDBFile('systems/glycopeptide.pdb')
        Modeller.loadHydrogenDefinitions('glycam-hydrogens.xml')
        modeller = Modeller(pdb.topology, pdb.positions)
        hydrogens = [a for a in modeller.topology.atoms() if a.element == element.hydrogen]
        modeller.delete(hydrogens)
        self.assertTrue(modeller.topology.getNumAtoms() < pdb.topology.getNumAtoms())
        modeller.addHydrogens()
        self.assertEqual(modeller.topology.getNumAtoms(), pdb.topology.getNumAtoms())
        for res1, res2 in zip(pdb.topology.residues(), modeller.topology.residues()):
            names1 = sorted([a.name for a in res1.atoms()])
            names2 = sorted([a.name for a in res2.atoms()])
            self.assertEqual(names1, names2)
        # Reset the loaded definitions so we don't affect other tests.
        Modeller._residueHydrogens = {}
        Modeller._hasLoadedStandardHydrogens = False

    def test_addSpecificHydrogens(self):
        """Test specifying exactly which hydrogens to add."""
        pdb = PDBFile('systems/glycopeptide.pdb')
        variants = [None]*pdb.topology.getNumResidues()
        for residue in pdb.topology.residues():
            if residue.name != 'ALA':
                var = []
                for atom1, atom2 in residue.bonds():
                    if atom1.element == element.hydrogen:
                        var.append((atom1.name, atom2.name))
                    elif atom2.element == element.hydrogen:
                        var.append((atom2.name, atom1.name))
                variants[residue.index] = var
        modeller = Modeller(pdb.topology, pdb.positions)
        hydrogens = [a for a in modeller.topology.atoms() if a.element == element.hydrogen and random.random() < 0.7]
        modeller.delete(hydrogens)
        self.assertTrue(modeller.topology.getNumAtoms() < pdb.topology.getNumAtoms())
        modeller.addHydrogens(variants=variants)
        self.assertEqual(modeller.topology.getNumAtoms(), pdb.topology.getNumAtoms())
        for res1, res2 in zip(pdb.topology.residues(), modeller.topology.residues()):
            names1 = sorted([a.name for a in res1.atoms()])
            names2 = sorted([a.name for a in res2.atoms()])
            self.assertEqual(names1, names2)

    def test_removeExtraHydrogens(self):
        """Test that addHydrogens() can remove hydrogens that shouldn't be there. """

        topology_start = self.topology_start3
        positions = self.positions3
        modeller = Modeller(topology_start, positions)

        # Add hydrogens, forcing residue 1 to be ASH.

        variants = [None]*25
        variants[1] = 'ASH'
        modeller.addHydrogens(self.forcefield, variants=variants)
        residue = list(modeller.topology.residues())[1]
        self.assertTrue(any(a.name == 'HD2' for a in residue.atoms()))

        # Now force it to be ASP and see if HD2 gets removed.

        variants[1] = 'ASP'
        modeller.addHydrogens(self.forcefield, variants=variants)
        residue = list(modeller.topology.residues())[1]
        self.assertFalse(any(a.name == 'HD2' for a in residue.atoms()))


    def test_addExtraParticles(self):
        """Test addExtraParticles()."""

        # Create a box of water.

        ff1 = ForceField('tip3p.xml')
        modeller = Modeller(Topology(), []*nanometers)
        modeller.addSolvent(ff1, 'tip3p', boxSize=Vec3(2,2,2)*nanometers)

        # Now convert the water to TIP4P.

        ff2 = ForceField('tip4pew.xml')
        modeller.addExtraParticles(ff2)
        for residue in modeller.topology.residues():
            atoms = list(residue.atoms())
            self.assertEqual(4, len(atoms))
            ep = [atom for atom in atoms if atom.element is None]
            self.assertEqual(1, len(ep))


    def test_addVirtualSites(self):
        """Test adding extra particles defined by virtual sites."""
        xml = """
            <ForceField>
             <AtomTypes>
              <Type name="C" class="C" element="C" mass="10"/>
              <Type name="N" class="N" element="N" mass="10"/>
              <Type name="O" class="O" element="O" mass="10"/>
              <Type name="V" class="V" mass="0.0"/>
             </AtomTypes>
             <Residues>
              <Residue name="Test">
               <Atom name="C" type="C"/>
               <Atom name="N" type="N"/>
               <Atom name="O" type="O"/>
               <Atom name="V1" type="V"/>
               <Atom name="V2" type="V"/>
               <Atom name="V3" type="V"/>
               <Atom name="V4" type="V"/>
               <VirtualSite type="average2" index="3" atom1="0" atom2="1" weight1="0.7" weight2="0.3"/>
               <VirtualSite type="average3" index="4" atom1="0" atom2="1" atom3="2" weight1="0.2" weight2="0.3" weight3="0.5"/>
               <VirtualSite type="outOfPlane" index="5" atom1="0" atom2="1" atom3="2" weight12="0.1" weight13="-0.2" weightCross="0.8"/>
               <VirtualSite type="localCoords" index="6" atom1="0" atom2="1" atom3="2" wo1="0.1" wo2="0.5" wo3="0.4" wx1="1" wx2="-0.6" wx3="-0.4" wy1="0.1" wy2="0.9" wy3="-1" p1="-0.5" p2="0.4" p3="1.1"/>
              </Residue>
             </Residues>
            </ForceField>"""
        ff = ForceField(StringIO(xml))

        # Create the three real atoms.

        topology = Topology()
        chain = topology.addChain()
        residue = topology.addResidue('Test', chain)
        topology.addAtom('C', element.carbon, residue)
        topology.addAtom('N', element.nitrogen, residue)
        topology.addAtom('V', element.oxygen, residue)


        # Add the virtual sites.

        modeller = Modeller(topology, [Vec3(0.1, 0.2, 0.3), Vec3(1.0, 0.9, 0.8), Vec3(1.5, 1.1, 0.7)]*nanometers)
        modeller.addExtraParticles(ff)
        top = modeller.topology
        pos = modeller.positions

        # Check that the correct particles were added.

        self.assertEqual(len(pos), 7)
        for atom, elem in zip(top.atoms(), [element.carbon, element.nitrogen, element.oxygen, None, None, None, None]):
            self.assertEqual(elem, atom.element)

        # Check that the positions were calculated correctly.

        system = ff.createSystem(top)
        integ = VerletIntegrator(1.0)
        context = Context(system, integ)
        context.setPositions(pos)
        context.computeVirtualSites()
        state = context.getState(getPositions=True)
        for p1, p2 in zip (pos, state.getPositions()):
            self.assertVecAlmostEqual(p1.value_in_unit(nanometers), p2.value_in_unit(nanometers), 1e-6)


    def test_multiSiteIon(self):
        """Test adding extra particles whose positions are determined based on bonds."""
        xml = """
            <ForceField>
             <AtomTypes>
              <Type name="Zn" class="Zn" element="Zn" mass="53.380"/>
              <Type name="DA" class="DA" mass="3.0"/>
             </AtomTypes>
             <Residues>
              <Residue name="ZN">
               <Atom name="ZN" type="Zn"/>
               <Atom name="D1" type="DA"/>
               <Atom name="D2" type="DA"/>
               <Atom name="D3" type="DA"/>
               <Atom name="D4" type="DA"/>
               <Bond from="0" to="2"/>
               <Bond from="0" to="1"/>
               <Bond from="0" to="3"/>
               <Bond from="0" to="4"/>
               <Bond from="1" to="2"/>
               <Bond from="1" to="3"/>
               <Bond from="1" to="4"/>
               <Bond from="2" to="4"/>
               <Bond from="2" to="3"/>
               <Bond from="3" to="4"/>
              </Residue>
             </Residues>
             <HarmonicBondForce>
              <Bond class1="DA" class2="Zn" length="0.09" k="535552.0"/>
              <Bond class1="DA" class2="DA" length="0.147" k="535552.0"/>
             </HarmonicBondForce>
            </ForceField>"""
        ff = ForceField(StringIO(xml))

        # Create two zinc atoms.

        topology = Topology()
        chain = topology.addChain()
        residue = topology.addResidue('ZN', chain)
        topology.addAtom('ZN', element.zinc, residue)
        residue = topology.addResidue('ZN', chain)
        topology.addAtom('ZN', element.zinc, residue)

        # Add the extra particles.

        modeller = Modeller(topology, [Vec3(0.5, 1.0, 1.5), Vec3(2.0, 2.0, 0.0)]*nanometers)
        modeller.addExtraParticles(ff)
        top = modeller.topology
        pos = modeller.positions

        # Check that the correct particles were added.

        self.assertEqual(len(pos), 10)
        for i, atom in enumerate(top.atoms()):
            self.assertEqual(element.zinc if i in (0,5) else None, atom.element)

        # Check that the positions in the first residue are reasonable.

        center = Vec3(0.5, 1.0, 1.5)*nanometers
        self.assertEqual(center, modeller.positions[0])
        for i in range(1, 5):
            for j in range(i):
                dist = norm(pos[i]-pos[j])
                expectedDist = 0.09 if j == 0 else 0.147
                self.assertTrue(dist > (expectedDist-0.01)*nanometers and dist < (expectedDist+0.01)*nanometers)


    def test_addMembrane(self):
        """Test adding a membrane to a realistic system."""

        mol = PDBxFile('systems/gpcr.cif')
        modeller = Modeller(mol.topology, mol.positions)
        ff = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

        # Add a membrane around the GPCR
        modeller.addMembrane(ff, minimumPadding=1.1*nanometers, ionicStrength=1*molar)

        # Make sure we added everything correctly
        resCount = defaultdict(int)
        for res in modeller.topology.residues():
            resCount[res.name] += 1

        self.assertEqual(16, resCount['ALA'])
        self.assertEqual(226, resCount['POP'])  # 2x128 - overlapping
        self.assertTrue(resCount['HOH'] > 1)

        deltaQ = resCount['CL'] - resCount['NA']
        self.assertEqual(deltaQ, 10)  # protein net q: +10

        # Check _addIons did the right thing.
        expected_ion_fraction = 1.0*molar/(55.4*molar)

        total_water = resCount['HOH']
        total_water_ions = resCount['HOH'] + resCount['CL'] + resCount['NA']

        # total_water_ions - protein charge
        expected_sodium = math.floor((total_water_ions-10)*expected_ion_fraction+0.5)
        expected_chlorine = expected_sodium + 10

        self.assertEqual(resCount['CL'], expected_chlorine)
        self.assertEqual(resCount['NA'], expected_sodium)

        # Check lipid numbering for repetitions
        lipidIdList = [(r.chain.id, r.id) for r in modeller.topology.residues()
                       if r.name == 'POP']
        self.assertEqual(len(lipidIdList), len(set(lipidIdList)))

        # Check dimensions to see if padding was respected
        originalSize = max(mol.positions) - min(mol.positions)
        newSize = modeller.topology.getUnitCellDimensions()
        for i in range(3):
            self.assertTrue(newSize[i] >= originalSize[i]+1.1*nanometers)

    def test_bondTypeAndOrderPreserved(self):
        """ Check that bond type and order are preserved across multiple operations. 

        Regression test for issue #4112 and similar behaviors.
        """

        # Given: an isolated molecule
        pdb = PDBFile("systems/alanine-dipeptide-implicit.pdb")
        topology, positions = pdb.topology, pdb.positions
        topology.setUnitCellDimensions(Vec3(3.5, 3.5, 3.5) * nanometers)
        # with some bonds carrying type and order information
        for bond in topology.bonds():
            if ((bond.atom1.element, bond.atom2.element) in [
                (element.carbon, element.oxygen), (element.oxygen, element.carbon)
            ]):
                bond.type = Double
                bond.order = 2.0
        modeller = Modeller(topology, positions)

        # When (1): add solvent
        forcefield = ForceField("amber10.xml", "tip3p.xml")
        modeller.addSolvent(forcefield, "tip3p")
        # sanity check: water was added
        self.assertTrue(any(r.name == "HOH" for r in modeller.topology.residues()))

        # When (2): convert water (no sites added)
        modeller.convertWater("spce")

        # When (3): convert water with addExtraParticles
        new_forcefield = ForceField('amber10.xml', 'tip4pew.xml')
        modeller.addExtraParticles(new_forcefield)
        # sanity check: extra sites were added
        self.assertEqual(
            set([len(list(res.atoms())) for res in modeller.topology.residues() if res.name == "HOH"]),
            {4}
        )

        # When (4): delete water (with deleteWater) and hydrogens (with delete)
        modeller.deleteWater()
        hydrogens = [a for a in modeller.topology.atoms() if a.element == element.hydrogen]
        modeller.delete(hydrogens)
        # sanity check: all gone
        self.assertFalse(any(a.element == element.hydrogen for a in modeller.topology.atoms()))
        self.assertFalse(any(r.name == "HOH" for r in modeller.topology.residues()))

        # When (5): add back hydrogens
        modeller.addHydrogens()
        # sanity check: hydrogens are back
        self.assertTrue(any(a.element == element.hydrogen for a in modeller.topology.atoms()))

        # Then (intermediate): bond info have been retained throughout the workflow
        self.assertIn((Double, 2.0), [(b.type, b.order) for b in modeller.topology.bonds()])

        # When (6): add a modeller (which also bears some bond info)
        to_add = PDBFile('systems/methanol-box.pdb')
        topology_to_add = to_add.topology
        positions_to_add = to_add.positions
        # add a dummy bond to the "to_add" system to check that it also is preserved
        atom0, atom1 = (atom for i, atom in enumerate(topology_to_add.atoms()) if i < 2)
        topology_to_add.addBond(atom0, atom1, Single, 1.0)
        modeller.add(topology_to_add, positions_to_add)

        # Then: bond info are retained for both the old and the new system
        all_bond_extra_data = [(b.type, b.order) for b in modeller.topology.bonds()]
        self.assertEqual(
            set(all_bond_extra_data),
            # None and Double from topology; other Nones and Single from topology_to_add
            {(None, None), (Single, 1.0), (Double, 2.0)}
        )

    def test_residueTemplates(self):
        """Test the residueTemplates argument to Modeller methods"""

        # Create a Topology and ForceField involving residues that match multiple templates.

        topology = Topology()
        chain = topology.addChain()
        residue1 = topology.addResidue('Fe', chain)
        topology.addAtom('Fe', element.iron, residue1)
        residue2 = topology.addResidue('Fe', chain)
        topology.addAtom('Fe', element.iron, residue2)
        positions = [Vec3(0, 0, 0), Vec3(1, 0, 0)]*nanometers
        ff = ForceField('amber14/tip3pfb.xml', 'amber14/lipid17.xml')
        residueTemplates = {residue1: 'FE2', residue2: 'FE'}

        # Test addSolvent().

        modeller = Modeller(topology, positions)
        with self.assertRaises(Exception):
            modeller.addSolvent(ff, padding=1*nanometer)
        modeller.addSolvent(ff, padding=1*nanometer, residueTemplates=residueTemplates)

        # Test addHydrogens().

        modeller = Modeller(topology, positions)
        with self.assertRaises(Exception):
            modeller.addHydrogens(ff)
        modeller.addHydrogens(ff, residueTemplates=residueTemplates)

        # Test addExtraParticles().

        modeller = Modeller(topology, positions)
        with self.assertRaises(Exception):
            modeller.addExtraParticles(ff)
        modeller.addExtraParticles(ff, residueTemplates=residueTemplates)

        # Test addMembrane().

        modeller = Modeller(topology, positions)
        with self.assertRaises(Exception):
            modeller.addMembrane(ff)
        modeller.addMembrane(ff, residueTemplates=residueTemplates)

    def assertVecAlmostEqual(self, p1, p2, tol=1e-7):
        scale = max(1.0, norm(p1),)
        for i in range(3):
            diff = abs(p1[i]-p2[i])/scale
            self.assertTrue(diff < tol)

if __name__ == '__main__':
    unittest.main()

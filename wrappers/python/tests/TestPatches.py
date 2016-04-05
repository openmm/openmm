import unittest
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.forcefield as forcefield
import math
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO
import os

class TestForceField(unittest.TestCase):
    """Test ForceFields that use patches."""

    def testParsePatch(self):
        """Test parsing a <Patch> tag."""

        xml = """
<ForceField>
 <AtomTypes>
  <Type name="A type" class="A class" element="O" mass="15.99943"/>
  <Type name="B type" class="B class" element="H" mass="1.007947"/>
 </AtomTypes>
 <Patches>
  <Patch name="Test">
    <AddAtom name="A" type="A type"/>
    <ChangeAtom name="B" type="B type"/>
    <RemoveAtom name="C"/>
    <AddBond atomName1="A" atomName2="B"/>
    <RemoveBond atomName1="B" atomName2="C"/>
    <AddExternalBond atomName="A"/>
    <RemoveExternalBond atomName="C"/>
    <ApplyToResidue name="RES"/>
  </Patch>
 </Patches>
</ForceField>"""
        ff = ForceField(StringIO(xml))
        self.assertEqual(1, len(ff._patches))
        patch = ff._patches['Test']
        self.assertEqual(1, len(patch.addedAtoms))
        self.assertEqual(1, len(patch.changedAtoms))
        self.assertEqual(1, len(patch.deletedAtoms))
        self.assertEqual(1, len(patch.addedBonds))
        self.assertEqual(1, len(patch.deletedBonds))
        self.assertEqual(1, len(patch.addedExternalBonds))
        self.assertEqual(1, len(patch.deletedExternalBonds))
        self.assertEqual(1, len(ff._templatePatches))
        self.assertEqual(1, len(ff._templatePatches['RES']))
        self.assertEqual('A', patch.addedAtoms[0].name.name)
        self.assertEqual(0, patch.addedAtoms[0].name.residue)
        self.assertEqual('A type', patch.addedAtoms[0].type)
        self.assertEqual('B', patch.changedAtoms[0].name.name)
        self.assertEqual(0, patch.changedAtoms[0].name.residue)
        self.assertEqual('B type', patch.changedAtoms[0].type)
        self.assertEqual('C', patch.deletedAtoms[0].name)
        self.assertEqual(0, patch.deletedAtoms[0].residue)
        self.assertEqual('A', patch.addedBonds[0][0].name)
        self.assertEqual(0, patch.addedBonds[0][0].residue)
        self.assertEqual('B', patch.addedBonds[0][1].name)
        self.assertEqual(0, patch.addedBonds[0][1].residue)
        self.assertEqual('B', patch.deletedBonds[0][0].name)
        self.assertEqual(0, patch.deletedBonds[0][0].residue)
        self.assertEqual('C', patch.deletedBonds[0][1].name)
        self.assertEqual(0, patch.deletedBonds[0][1].residue)
        self.assertEqual('A', patch.addedExternalBonds[0].name)
        self.assertEqual(0, patch.addedExternalBonds[0].residue)
        self.assertEqual('C', patch.deletedExternalBonds[0].name)
        self.assertEqual(0, patch.deletedExternalBonds[0].residue)
        self.assertEqual('Test', ff._templatePatches['RES'][0][0])
        self.assertEqual(0, ff._templatePatches['RES'][0][1])

if __name__ == '__main__':
    unittest.main()

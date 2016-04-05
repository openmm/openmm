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
        self.assertEqual(1, len(patch.addedAtoms[0]))
        self.assertEqual(1, len(patch.changedAtoms))
        self.assertEqual(1, len(patch.changedAtoms[0]))
        self.assertEqual(1, len(patch.deletedAtoms))
        self.assertEqual(1, len(patch.addedBonds))
        self.assertEqual(1, len(patch.deletedBonds))
        self.assertEqual(1, len(patch.addedExternalBonds))
        self.assertEqual(1, len(patch.deletedExternalBonds))
        self.assertEqual(1, len(ff._templatePatches))
        self.assertEqual(1, len(ff._templatePatches['RES']))
        self.assertEqual('A', patch.addedAtoms[0][0].name)
        self.assertEqual('A type', patch.addedAtoms[0][0].type)
        self.assertEqual('B', patch.changedAtoms[0][0].name)
        self.assertEqual('B type', patch.changedAtoms[0][0].type)
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

    def testParseMultiresiduePatch(self):
        """Test parsing a <Patch> tag that affects two residues."""

        xml = """
<ForceField>
 <AtomTypes>
  <Type name="A type" class="A class" element="O" mass="15.99943"/>
  <Type name="B type" class="B class" element="H" mass="1.007947"/>
 </AtomTypes>
 <Patches>
  <Patch name="Test" residues="2">
    <AddAtom name="1:A" type="A type"/>
    <ChangeAtom name="2:B" type="B type"/>
    <AddBond atomName1="1:A" atomName2="2:B"/>
    <ApplyToResidue name="1:RESA"/>
    <ApplyToResidue name="2:RESB"/>
  </Patch>
 </Patches>
</ForceField>"""
        ff = ForceField(StringIO(xml))
        self.assertEqual(1, len(ff._patches))
        patch = ff._patches['Test']
        self.assertEqual(2, len(patch.addedAtoms))
        self.assertEqual(1, len(patch.addedAtoms[0]))
        self.assertEqual(0, len(patch.addedAtoms[1]))
        self.assertEqual(2, len(patch.changedAtoms))
        self.assertEqual(0, len(patch.changedAtoms[0]))
        self.assertEqual(1, len(patch.changedAtoms[1]))
        self.assertEqual(1, len(patch.addedBonds))
        self.assertEqual(2, len(ff._templatePatches))
        self.assertEqual(1, len(ff._templatePatches['RESA']))
        self.assertEqual(1, len(ff._templatePatches['RESB']))
        self.assertEqual('A', patch.addedAtoms[0][0].name)
        self.assertEqual('A type', patch.addedAtoms[0][0].type)
        self.assertEqual('B', patch.changedAtoms[1][0].name)
        self.assertEqual('B type', patch.changedAtoms[1][0].type)
        self.assertEqual('A', patch.addedBonds[0][0].name)
        self.assertEqual(0, patch.addedBonds[0][0].residue)
        self.assertEqual('B', patch.addedBonds[0][1].name)
        self.assertEqual(1, patch.addedBonds[0][1].residue)
        self.assertEqual('Test', ff._templatePatches['RESA'][0][0])
        self.assertEqual(0, ff._templatePatches['RESA'][0][1])
        self.assertEqual('Test', ff._templatePatches['RESB'][0][0])
        self.assertEqual(1, ff._templatePatches['RESB'][0][1])

    def testApplyPatch(self):
        """Test applying a patch to a template."""

        xml = """
<ForceField>
 <AtomTypes>
  <Type name="A type" class="A class" element="O" mass="15.99943"/>
  <Type name="B type" class="B class" element="H" mass="1.007947"/>
  <Type name="C type" class="C class" element="H" mass="1.007947"/>
  <Type name="D type" class="D class" element="C" mass="12.010000"/>
 </AtomTypes>
 <Residues>
  <Residue name="RES">
   <Atom name="A" type="A type"/>
   <Atom name="B" type="B type"/>
   <Atom name="C" type="C type"/>
   <Bond atomName1="A" atomName2="B"/>
   <Bond atomName1="B" atomName2="C"/>
   <ExternalBond atomName="C"/>
   <VirtualSite type="average2" siteName="C" atomName1="B" atomName2="C" weight1="0.6" weight2="0.4"/>
  </Residue>
 </Residues>
 <Patches>
  <Patch name="Test">
    <AddAtom name="D" type="D type"/>
    <ChangeAtom name="B" type="A type"/>
    <RemoveAtom name="A"/>
    <AddBond atomName1="B" atomName2="D"/>
    <RemoveBond atomName1="A" atomName2="B"/>
    <AddExternalBond atomName="D"/>
    <RemoveExternalBond atomName="C"/>
    <ApplyToResidue name="RES"/>
  </Patch>
 </Patches>
</ForceField>"""
        ff = ForceField(StringIO(xml))
        self.assertEqual(1, len(ff._patches))
        patch = ff._patches['Test']
        template = ff._templates['RES']
        newTemplates = patch.createPatchedTemplates([template])
        self.assertEqual(1, len(newTemplates))
        t = newTemplates[0]
        self.assertEqual(3, len(t.atoms))
        self.assertTrue(any(a.name == 'B' and a.type == 'A type' for a in t.atoms))
        self.assertTrue(any(a.name == 'C' and a.type == 'C type' for a in t.atoms))
        self.assertTrue(any(a.name == 'D' and a.type == 'D type' for a in t.atoms))
        indexMap = dict([(a.name, i) for i, a in enumerate(t.atoms)])
        self.assertEqual(2, len(t.bonds))
        self.assertTrue((indexMap['B'], indexMap['C']) in t.bonds)
        self.assertTrue((indexMap['B'], indexMap['D']) in t.bonds)
        self.assertEqual(1, len(t.externalBonds))
        self.assertTrue(indexMap['D'] in t.externalBonds)
        self.assertEqual(1, len(t.virtualSites))
        v = t.virtualSites[0]
        self.assertEqual('average2', v.type)
        self.assertEqual(0.6, v.weights[0])
        self.assertEqual(0.4, v.weights[1])
        self.assertEqual(indexMap['C'], v.index)
        self.assertEqual(indexMap['B'], v.atoms[0])
        self.assertEqual(indexMap['C'], v.atoms[1])

if __name__ == '__main__':
    unittest.main()

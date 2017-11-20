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

class TestPatches(unittest.TestCase):
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
        patch = list(ff._templatePatches['RES'])[0]
        self.assertEqual('Test', patch[0])
        self.assertEqual(0, patch[1])

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
        patchA = list(ff._templatePatches['RESA'])[0]
        self.assertEqual('Test', patchA[0])
        self.assertEqual(0, patchA[1])
        patchB = list(ff._templatePatches['RESB'])[0]
        self.assertEqual('Test', patchB[0])
        self.assertEqual(1, patchB[1])

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

    def testAlaAlaAla(self):
        """Test constructing a System that involves two patches."""

        xml = """
<ForceField>
 <AtomTypes>
  <Type name="N" class="N" element="N" mass="14.00672"/>
  <Type name="H" class="H" element="H" mass="1.007947"/>
  <Type name="CT" class="CT" element="C" mass="12.01078"/>
  <Type name="H1" class="H1" element="H" mass="1.007947"/>
  <Type name="HC" class="HC" element="H" mass="1.007947"/>
  <Type name="C" class="C" element="C" mass="12.01078"/>
  <Type name="O" class="O" element="O" mass="15.99943"/>
  <Type name="O2" class="O2" element="O" mass="15.99943"/>
  <Type name="N3" class="N3" element="N" mass="14.00672"/>
 </AtomTypes>
 <Residues>
  <Residue name="ALA">
   <Atom name="N" type="N"/>
   <Atom name="H" type="H"/>
   <Atom name="CA" type="CT"/>
   <Atom name="HA" type="H1"/>
   <Atom name="CB" type="CT"/>
   <Atom name="HB1" type="HC"/>
   <Atom name="HB2" type="HC"/>
   <Atom name="HB3" type="HC"/>
   <Atom name="C" type="C"/>
   <Atom name="O" type="O"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="8"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="8" to="9"/>
   <ExternalBond from="0"/>
   <ExternalBond from="8"/>
   <AllowPatch name="CTER"/>
   <AllowPatch name="NTER"/>
  </Residue>
 </Residues>
 <Patches>
  <Patch name="CTER">
    <AddAtom name="OXT" type="O2"/>
    <ChangeAtom name="O" type="O2"/>
    <AddBond atomName1="C" atomName2="OXT"/>
    <RemoveExternalBond atomName="C"/>
  </Patch>
  <Patch name="NTER">
    <RemoveAtom name="H"/>
    <AddAtom name="H1" type="H"/>
    <AddAtom name="H2" type="H"/>
    <AddAtom name="H3" type="H"/>
    <ChangeAtom name="N" type="N3"/>
    <RemoveBond atomName1="N" atomName2="H"/>
    <AddBond atomName1="N" atomName2="H1"/>
    <AddBond atomName1="N" atomName2="H2"/>
    <AddBond atomName1="N" atomName2="H3"/>
    <RemoveExternalBond atomName="N"/>
  </Patch>
 </Patches>
 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
  <Atom type="N" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="H" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="CT" charge="0.0337" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="H1" charge="0.0823" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="HC" charge="0.0603" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="C" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="O" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="O2" charge="-0.8055" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="N3" charge="0.1414" sigma="0.324999852378" epsilon="0.71128"/>
 </NonbondedForce>
</ForceField>"""
        ff = ForceField(StringIO(xml))
        pdb = PDBFile(os.path.join('systems', 'ala_ala_ala.pdb'))
        system = ff.createSystem(pdb.topology)
        nb = system.getForce(0)
        expectedCharges = [0.1414, 0.2719, 0.2719, 0.2719, 0.0337, 0.0823, 0.0337, 0.0603, 0.0603, 0.0603, 0.5973, -0.5679,
                           -0.4157, 0.2719, 0.0337, 0.0823, 0.0337, 0.0603, 0.0603, 0.0603, 0.5973, -0.5679,
                           0.5973, -0.8055, -0.8055, -0.4157, 0.2719, 0.0337, 0.0823, 0.0337, 0.0603, 0.0603, 0.0603]
        for i in range(system.getNumParticles()):
            self.assertEqual(expectedCharges[i], nb.getParticleParameters(i)[0].value_in_unit(elementary_charge))

    def testDisulfidePatch(self):
        pdb = PDBFile(os.path.join('systems', 'bpti.pdb'))
        ff = ForceField('amber99sb.xml')
        system1 = ff.createSystem(pdb.topology)
        
        # Override the CYX template so it will no longer match.
        
        xml = """
<ForceField>
 <Residues>
  <Residue name="CYX" override="1">
  </Residue>
 </Residues>
</ForceField>"""
        ff.loadFile(StringIO(xml))
        try:
            ff.createSystem(pdb.topology)
            failed = False
        except:
            failed = True
        self.assertTrue(failed)
        
        # Now add a patch for matching disulfides.
        
        xml = """
<ForceField>
 <Patches>
  <Patch name="Disulfide" residues="2">
    <RemoveAtom name="1:HG"/>
    <RemoveAtom name="2:HG"/>
    <ChangeAtom name="1:CA" type="83"/>
    <ChangeAtom name="2:CA" type="83"/>
    <ChangeAtom name="1:HA" type="84"/>
    <ChangeAtom name="2:HA" type="84"/>
    <ChangeAtom name="1:CB" type="85"/>
    <ChangeAtom name="2:CB" type="85"/>
    <ChangeAtom name="1:HB2" type="86"/>
    <ChangeAtom name="2:HB2" type="86"/>
    <ChangeAtom name="1:HB3" type="86"/>
    <ChangeAtom name="2:HB3" type="86"/>
    <ChangeAtom name="1:SG" type="87"/>
    <ChangeAtom name="2:SG" type="87"/>
    <AddBond atomName1="1:SG" atomName2="2:SG"/>
    <ApplyToResidue name="1:CYS"/>
    <ApplyToResidue name="2:CYS"/>
  </Patch>
 </Patches>
</ForceField>"""
        ff.loadFile(StringIO(xml))
        system2 = ff.createSystem(pdb.topology)
        self.assertEqual(XmlSerializer.serialize(system1), XmlSerializer.serialize(system2))

if __name__ == '__main__':
    unittest.main()

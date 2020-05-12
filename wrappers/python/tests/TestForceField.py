import unittest
from validateConstraints import *
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem
import simtk.openmm.app.forcefield as forcefield
import math
import textwrap
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO
import os
import warnings

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

        # alanine dipeptide with implicit water
        self.pdb2 = PDBFile('systems/alanine-dipeptide-implicit.pdb')
        self.forcefield2 = ForceField('amber99sb.xml', 'amber99_obc.xml')


    def test_NonbondedMethod(self):
        """Test all six options for the nonbondedMethod parameter."""

        methodMap = {NoCutoff:NonbondedForce.NoCutoff,
                     CutoffNonPeriodic:NonbondedForce.CutoffNonPeriodic,
                     CutoffPeriodic:NonbondedForce.CutoffPeriodic,
                     Ewald:NonbondedForce.Ewald,
                     PME:NonbondedForce.PME,
                     LJPME:NonbondedForce.LJPME}
        for method in methodMap:
            system = self.forcefield1.createSystem(self.pdb1.topology,
                                                  nonbondedMethod=method)
            forces = system.getForces()
            self.assertTrue(any(isinstance(f, NonbondedForce) and
                                f.getNonbondedMethod()==methodMap[method]
                                for f in forces))

    def test_DispersionCorrection(self):
        """Test to make sure that the dispersion/long-range correction is set properly."""
        top = Topology()
        chain = top.addChain()

        for lrc in (True, False):
            xml = textwrap.dedent(
                """
                <ForceField>
                 <LennardJonesForce lj14scale="0.3" useDispersionCorrection="{lrc}">
                  <Atom type="A" sigma="1" epsilon="0.1"/>
                  <Atom type="B" sigma="2" epsilon="0.2"/>
                  <NBFixPair type1="A" type2="B" sigma="2.5" epsilon="1.1"/>
                 </LennardJonesForce>
                 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5" useDispersionCorrection="{lrc2}">
                  <Atom type="A" sigma="0.315" epsilon="0.635"/>
                 </NonbondedForce>
                </ForceField>
                """
            )
            ff = ForceField(StringIO(xml.format(lrc=lrc, lrc2=lrc)))
            system = ff.createSystem(top)
            checked_nonbonded = False
            checked_custom = False
            for force in system.getForces():
                if isinstance(force, NonbondedForce):
                    self.assertEqual(force.getUseDispersionCorrection(), lrc)
                    checked_nonbonded = True
                elif isinstance(force, CustomNonbondedForce):
                    self.assertEqual(force.getUseLongRangeCorrection(), lrc)
                    checked_custom = True
            self.assertTrue(checked_nonbonded and checked_custom)

            # check that the keyword argument overwrites xml input
            lrc_kwarg = not lrc
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                system2 = ff.createSystem(top, useDispersionCorrection=lrc_kwarg)
                self.assertTrue(len(w) == 2)
                assert "conflict" in str(w[-1].message).lower()
            checked_nonbonded = False
            checked_custom = False
            for force in system2.getForces():
                if isinstance(force, NonbondedForce):
                    self.assertEqual(force.getUseDispersionCorrection(), lrc_kwarg)
                    checked_nonbonded = True
                elif isinstance(force, CustomNonbondedForce):
                    self.assertEqual(force.getUseLongRangeCorrection(), lrc_kwarg)
                    checked_custom = True
            self.assertTrue(checked_nonbonded and checked_custom)

            # check that no warning is generated when useDispersionCorrection is not in the xml file
            xml = textwrap.dedent(
                """
                <ForceField>
                 <LennardJonesForce lj14scale="0.3">
                  <Atom type="A" sigma="1" epsilon="0.1"/>
                  <Atom type="B" sigma="2" epsilon="0.2"/>
                  <NBFixPair type1="A" type2="B" sigma="2.5" epsilon="1.1"/>
                 </LennardJonesForce>
                 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
                  <Atom type="A" sigma="0.315" epsilon="0.635"/>
                 </NonbondedForce>
                </ForceField>
                """
            )
            ff = ForceField(StringIO(xml))
            system = ff.createSystem(top)
            for lrc_kwarg in [True, False]:
                with warnings.catch_warnings():
                    warnings.simplefilter("error")
                    system2 = ff.createSystem(top, useDispersionCorrection=lrc_kwarg)

    def test_Cutoff(self):
        """Test to make sure the nonbondedCutoff parameter is passed correctly."""

        for method in [CutoffNonPeriodic, CutoffPeriodic, Ewald, PME, LJPME]:
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

    def test_SwitchingDistance(self):
        """Test that the switchDistance parameter is processed correctly."""

        for switchDistance in [None, 0.9*nanometers]:
            system = self.forcefield1.createSystem(self.pdb1.topology,
                                                   nonbondedMethod=PME,
                                                   switchDistance=switchDistance)
            for force in system.getForces():
                if isinstance(force, NonbondedForce):
                    if switchDistance is None:
                        self.assertFalse(force.getUseSwitchingFunction())
                    else:
                        self.assertTrue(force.getUseSwitchingFunction())
                        self.assertEqual(switchDistance, force.getSwitchingDistance())

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

    def test_flexibleConstraints(self):
        """ Test the flexibleConstraints keyword """
        topology = self.pdb1.topology
        system1 = self.forcefield1.createSystem(topology, constraints=HAngles,
                                                rigidWater=True)
        system2 = self.forcefield1.createSystem(topology, constraints=HAngles,
                                                rigidWater=True, flexibleConstraints=True)
        system3 = self.forcefield1.createSystem(topology, constraints=None, rigidWater=False)
        validateConstraints(self, topology, system1, HAngles, True)
        # validateConstraints fails for system2 since by definition atom pairs can be in both bond
        # and constraint lists. So just check that the number of constraints is the same for both
        # system1 and system2
        self.assertEqual(system1.getNumConstraints(), system2.getNumConstraints())
        for force in system1.getForces():
            if isinstance(force, HarmonicBondForce):
                bf1 = force
            elif isinstance(force, HarmonicAngleForce):
                af1 = force
        for force in system2.getForces():
            if isinstance(force, HarmonicBondForce):
                bf2 = force
            elif isinstance(force, HarmonicAngleForce):
                af2 = force
        for force in system3.getForces():
            if isinstance(force, HarmonicAngleForce):
                af3 = force
        # Make sure we picked up extra bond terms with flexibleConstraints
        self.assertGreater(bf2.getNumBonds(), bf1.getNumBonds())
        # Make sure flexibleConstraints yields just as many angles as no constraints
        self.assertEqual(af2.getNumAngles(), af3.getNumAngles())

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
        with open('systems/lysozyme-implicit-forces.xml') as input:
            state2 = XmlSerializer.deserialize(input.read())
        numDifferences = 0
        for f1, f2, in zip(state1.getForces().value_in_unit(kilojoules_per_mole/nanometer), state2.getForces().value_in_unit(kilojoules_per_mole/nanometer)):
            diff = norm(f1-f2)
            if diff > 0.1 and diff/norm(f1) > 1e-3:
                numDifferences += 1
        self.assertTrue(numDifferences < system.getNumParticles()/20) # Tolerate occasional differences from numerical error

    def test_ProgrammaticForceField(self):
        """Test building a ForceField programmatically."""

        # Build the ForceField for TIP3P programmatically.
        ff = ForceField()
        ff.registerAtomType({'name':'tip3p-O', 'class':'OW', 'mass':15.99943*daltons, 'element':elem.oxygen})
        ff.registerAtomType({'name':'tip3p-H', 'class':'HW', 'mass':1.007947*daltons, 'element':elem.hydrogen})
        residue = ForceField._TemplateData('HOH')
        residue.atoms.append(ForceField._TemplateAtomData('O', 'tip3p-O', elem.oxygen))
        residue.atoms.append(ForceField._TemplateAtomData('H1', 'tip3p-H', elem.hydrogen))
        residue.atoms.append(ForceField._TemplateAtomData('H2', 'tip3p-H', elem.hydrogen))
        residue.addBond(0, 1)
        residue.addBond(0, 2)
        ff.registerResidueTemplate(residue)
        bonds = forcefield.HarmonicBondGenerator(ff)
        bonds.registerBond({'class1':'OW', 'class2':'HW', 'length':0.09572*nanometers, 'k':462750.4*kilojoules_per_mole/nanometer})
        ff.registerGenerator(bonds)
        angles = forcefield.HarmonicAngleGenerator(ff)
        angles.registerAngle({'class1':'HW', 'class2':'OW', 'class3':'HW', 'angle':1.82421813418*radians, 'k':836.8*kilojoules_per_mole/radian})
        ff.registerGenerator(angles)
        nonbonded = forcefield.NonbondedGenerator(ff, 0.833333, 0.5, True)
        nonbonded.registerAtom({'type':'tip3p-O', 'charge':-0.834, 'sigma':0.31507524065751241*nanometers, 'epsilon':0.635968*kilojoules_per_mole})
        nonbonded.registerAtom({'type':'tip3p-H', 'charge':0.417, 'sigma':1*nanometers, 'epsilon':0*kilojoules_per_mole})
        ff.registerGenerator(nonbonded)

        # Build a water box.
        modeller = Modeller(Topology(), [])
        modeller.addSolvent(ff, boxSize=Vec3(3, 3, 3)*nanometers)

        # Create a system using the programmatic force field as well as one from an XML file.
        system1 = ff.createSystem(modeller.topology)
        ff2 = ForceField('tip3p.xml')
        system2 = ff2.createSystem(modeller.topology)
        self.assertEqual(XmlSerializer.serialize(system1), XmlSerializer.serialize(system2))

    def test_PeriodicBoxVectors(self):
        """Test setting the periodic box vectors."""

        vectors = (Vec3(5, 0, 0), Vec3(-1.5, 4.5, 0), Vec3(0.4, 0.8, 7.5))*nanometers
        self.pdb1.topology.setPeriodicBoxVectors(vectors)
        self.assertEqual(Vec3(5, 4.5, 7.5)*nanometers, self.pdb1.topology.getUnitCellDimensions())
        system = self.forcefield1.createSystem(self.pdb1.topology)
        for i in range(3):
            self.assertEqual(vectors[i], self.pdb1.topology.getPeriodicBoxVectors()[i])
            self.assertEqual(vectors[i], system.getDefaultPeriodicBoxVectors()[i])

    def test_ResidueAttributes(self):
        """Test a ForceField that gets per-particle parameters from residue attributes."""

        xml = """
<ForceField>
 <AtomTypes>
  <Type name="tip3p-O" class="OW" element="O" mass="15.99943"/>
  <Type name="tip3p-H" class="HW" element="H" mass="1.007947"/>
 </AtomTypes>
 <Residues>
  <Residue name="HOH">
   <Atom name="O" type="tip3p-O" charge="-0.834"/>
   <Atom name="H1" type="tip3p-H" charge="0.417"/>
   <Atom name="H2" type="tip3p-H" charge="0.417"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
  </Residue>
 </Residues>
 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
  <UseAttributeFromResidue name="charge"/>
  <Atom type="tip3p-O" sigma="0.315" epsilon="0.635"/>
  <Atom type="tip3p-H" sigma="1" epsilon="0"/>
 </NonbondedForce>
</ForceField>"""
        ff = ForceField(StringIO(xml))

        # Build a water box.
        modeller = Modeller(Topology(), [])
        modeller.addSolvent(ff, boxSize=Vec3(3, 3, 3)*nanometers)

        # Create a system and make sure all nonbonded parameters are correct.
        system = ff.createSystem(modeller.topology)
        nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
        atoms = list(modeller.topology.atoms())
        for i in range(len(atoms)):
            params = nonbonded.getParticleParameters(i)
            if atoms[i].element == elem.oxygen:
                self.assertEqual(params[0], -0.834*elementary_charge)
                self.assertEqual(params[1], 0.315*nanometers)
                self.assertEqual(params[2], 0.635*kilojoule_per_mole)
            else:
                self.assertEqual(params[0], 0.417*elementary_charge)
                self.assertEqual(params[1], 1.0*nanometers)
                self.assertEqual(params[2], 0.0*kilojoule_per_mole)

    def test_residueTemplateGenerator(self):
        """Test the ability to add residue template generators to parameterize unmatched residues."""
        def simpleTemplateGenerator(forcefield, residue):
            """\
            Simple residue template generator.
            This implementation uses the programmatic API to define residue templates.

            NOTE: We presume we have already loaded the force definitions into ForceField.
            """
            # Generate a unique prefix name for generating parameters.
            from uuid import uuid4
            template_name = uuid4()
            # Create residue template.
            from simtk.openmm.app.forcefield import _createResidueTemplate
            template = _createResidueTemplate(residue) # use helper function
            template.name = template_name # replace template name
            for (template_atom, residue_atom) in zip(template.atoms, residue.atoms()):
                template_atom.type = 'XXX' # replace atom type
            # Register the template.
            forcefield.registerResidueTemplate(template)

            # Signal that we have successfully parameterized the residue.
            return True

        # Define forcefield parameters used by simpleTemplateGenerator.
        # NOTE: This parameter definition file will currently only work for residues that either have
        # no external bonds or external bonds to other residues parameterized by the simpleTemplateGenerator.
        simple_ffxml_contents = """
<ForceField>
 <AtomTypes>
  <Type name="XXX" class="XXX" element="C" mass="12.0"/>
 </AtomTypes>
 <HarmonicBondForce>
  <Bond type1="XXX" type2="XXX" length="0.1409" k="392459.2"/>
 </HarmonicBondForce>
 <HarmonicAngleForce>
  <Angle type1="XXX" type2="XXX" type3="XXX" angle="2.09439510239" k="527.184"/>
 </HarmonicAngleForce>
 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
  <Atom type="XXX" charge="0.000" sigma="0.315" epsilon="0.635"/>
 </NonbondedForce>
</ForceField>"""

        #
        # Test where we generate parameters for only a ligand.
        #

        # Load the PDB file.
        pdb = PDBFile(os.path.join('systems', 'T4-lysozyme-L99A-p-xylene-implicit.pdb'))
        # Create a ForceField object.
        forcefield = ForceField('amber99sb.xml', 'tip3p.xml', StringIO(simple_ffxml_contents))
        # Add the residue template generator.
        forcefield.registerTemplateGenerator(simpleTemplateGenerator)
        # Parameterize system.
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
        # TODO: Test energies are finite?

        #
        # Test for a few systems where we generate all parameters.
        #

        tests = [
            { 'pdb_filename' : 'alanine-dipeptide-implicit.pdb', 'nonbondedMethod' : NoCutoff },
            { 'pdb_filename' : 'lysozyme-implicit.pdb', 'nonbondedMethod' : NoCutoff },
            { 'pdb_filename' : 'alanine-dipeptide-explicit.pdb', 'nonbondedMethod' : CutoffPeriodic },
            ]

        # Test all systems with separate ForceField objects.
        for test in tests:
            # Load the PDB file.
            pdb = PDBFile(os.path.join('systems', test['pdb_filename']))
            # Create a ForceField object.
            forcefield = ForceField(StringIO(simple_ffxml_contents))
            # Add the residue template generator.
            forcefield.registerTemplateGenerator(simpleTemplateGenerator)
            # Parameterize system.
            system = forcefield.createSystem(pdb.topology, nonbondedMethod=test['nonbondedMethod'])
            # TODO: Test energies are finite?

        # Now test all systems with a single ForceField object.
        # Create a ForceField object.
        forcefield = ForceField(StringIO(simple_ffxml_contents))
        # Add the residue template generator.
        forcefield.registerTemplateGenerator(simpleTemplateGenerator)
        for test in tests:
            # Load the PDB file.
            pdb = PDBFile(os.path.join('systems', test['pdb_filename']))
            # Parameterize system.
            system = forcefield.createSystem(pdb.topology, nonbondedMethod=test['nonbondedMethod'])
            # TODO: Test energies are finite?

    def test_getUnmatchedResidues(self):
        """Test retrieval of list of residues for which no templates are available."""

        # Load the PDB file.
        pdb = PDBFile(os.path.join('systems', 'T4-lysozyme-L99A-p-xylene-implicit.pdb'))
        # Create a ForceField object.
        forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
        # Get list of unmatched residues.
        unmatched_residues = forcefield.getUnmatchedResidues(pdb.topology)
        # Check results.
        self.assertEqual(len(unmatched_residues), 1)
        self.assertEqual(unmatched_residues[0].name, 'TMP')
        self.assertEqual(unmatched_residues[0].id, '163')

        # Load the PDB file.
        pdb = PDBFile(os.path.join('systems', 'ala_ala_ala.pdb'))
        # Create a ForceField object.
        forcefield = ForceField('tip3p.xml')
        # Get list of unmatched residues.
        unmatched_residues = forcefield.getUnmatchedResidues(pdb.topology)
        # Check results.
        self.assertEqual(len(unmatched_residues), 3)
        self.assertEqual(unmatched_residues[0].name, 'ALA')
        self.assertEqual(unmatched_residues[0].chain.id, 'X')
        self.assertEqual(unmatched_residues[0].id, '1')

    def test_ggenerateTemplatesForUnmatchedResidues(self):
        """Test generation of blank forcefield residue templates for unmatched residues."""
        #
        # Test where we generate parameters for only a ligand.
        #

        # Load the PDB file.
        pdb = PDBFile(os.path.join('systems', 'nacl-water.pdb'))
        # Create a ForceField object.
        forcefield = ForceField('tip3p.xml')
        # Get list of unmatched residues.
        unmatched_residues = forcefield.getUnmatchedResidues(pdb.topology)
        [templates, residues] = forcefield.generateTemplatesForUnmatchedResidues(pdb.topology)
        # Check results.
        self.assertEqual(len(unmatched_residues), 24)
        self.assertEqual(len(residues), 2)
        self.assertEqual(len(templates), 2)
        unique_names = set([ residue.name for residue in residues ])
        self.assertTrue('HOH' not in unique_names)
        self.assertTrue('NA' in unique_names)
        self.assertTrue('CL' in unique_names)
        template_names = set([ template.name for template in templates ])
        self.assertTrue('HOH' not in template_names)
        self.assertTrue('NA' in template_names)
        self.assertTrue('CL' in template_names)

        # Define forcefield parameters using returned templates.
        # NOTE: This parameter definition file will currently only work for residues that either have
        # no external bonds or external bonds to other residues parameterized by the simpleTemplateGenerator.
        simple_ffxml_contents = """
<ForceField>
 <AtomTypes>
  <Type name="XXX" class="XXX" element="C" mass="12.0"/>
 </AtomTypes>
 <HarmonicBondForce>
  <Bond type1="XXX" type2="XXX" length="0.1409" k="392459.2"/>
 </HarmonicBondForce>
 <HarmonicAngleForce>
  <Angle type1="XXX" type2="XXX" type3="XXX" angle="2.09439510239" k="527.184"/>
 </HarmonicAngleForce>
 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
  <Atom type="XXX" charge="0.000" sigma="0.315" epsilon="0.635"/>
 </NonbondedForce>
</ForceField>"""

        #
        # Test the pre-geenration of missing residue template for a ligand.
        #

        # Load the PDB file.
        pdb = PDBFile(os.path.join('systems', 'T4-lysozyme-L99A-p-xylene-implicit.pdb'))
        # Create a ForceField object.
        forcefield = ForceField('amber99sb.xml', 'tip3p.xml', StringIO(simple_ffxml_contents))
        # Get list of unique unmatched residues.
        [templates, residues] = forcefield.generateTemplatesForUnmatchedResidues(pdb.topology)
        # Add residue templates to forcefield.
        for template in templates:
            # Replace atom types.
            for atom in template.atoms:
                atom.type = 'XXX'
            # Register the template.
            forcefield.registerResidueTemplate(template)
        # Parameterize system.
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
        # TODO: Test energies are finite?

    def test_getMatchingTemplates(self):
        """Test retrieval of list of templates that match residues in a topology."""

        # Load the PDB file.
        pdb = PDBFile(os.path.join('systems', 'ala_ala_ala.pdb'))
        # Create a ForceField object.
        forcefield = ForceField('amber99sb.xml')
        # Get list of matching residue templates.
        templates = forcefield.getMatchingTemplates(pdb.topology)
        # Check results.
        residues = [ residue for residue in pdb.topology.residues() ]
        self.assertEqual(len(templates), len(residues))
        self.assertEqual(templates[0].name, 'NALA')
        self.assertEqual(templates[1].name, 'ALA')
        self.assertEqual(templates[2].name, 'CALA')

    def test_Wildcard(self):
        """Test that PeriodicTorsionForces using wildcard ('') for atom types / classes in the ffxml are correctly registered"""

        # Use wildcards in types
        xml = """
<ForceField>
 <AtomTypes>
  <Type name="C" class="C" element="C" mass="12.010000"/>
  <Type name="O" class="O" element="O" mass="16.000000"/>
 </AtomTypes>
 <PeriodicTorsionForce>
  <Proper type1="" type2="C" type3="C" type4="" periodicity1="2" phase1="3.141593" k1="15.167000"/>
  <Improper type1="C" type2="" type3="" type4="O" periodicity1="2" phase1="3.141593" k1="43.932000"/>
 </PeriodicTorsionForce>
</ForceField>"""

        ff = ForceField(StringIO(xml))

        self.assertEqual(len(ff._forces[0].proper), 1)
        self.assertEqual(len(ff._forces[0].improper), 1)

       # Use wildcards in classes
        xml = """
<ForceField>
 <AtomTypes>
  <Type name="C" class="C" element="C" mass="12.010000"/>
  <Type name="O" class="O" element="O" mass="16.000000"/>
 </AtomTypes>
 <PeriodicTorsionForce>
  <Proper class1="" class2="C" class3="C" class4="" periodicity1="2" phase1="3.141593" k1="15.167000"/>
  <Improper class1="C" class2="" class3="" class4="O" periodicity1="2" phase1="3.141593" k1="43.932000"/>
 </PeriodicTorsionForce>
</ForceField>"""

        ff = ForceField(StringIO(xml))

        self.assertEqual(len(ff._forces[0].proper), 1)
        self.assertEqual(len(ff._forces[0].improper), 1)

    def test_ScalingFactorCombining(self):
        """ Tests that FFs can be combined if their scaling factors are very close """
        forcefield = ForceField('amber99sb.xml', os.path.join('systems', 'test_amber_ff.xml'))
        # This would raise an exception if it didn't work

    def test_MultipleFilesandForceTags(self):
        """Test that the order of listing of multiple ffxmls does not matter.
           Tests that one generator per force type is created and that the ffxml
           defining atom types does not have to be listed first"""

        ffxml = """<ForceField>
 <Residues>
  <Residue name="ACE-Test">
   <Atom name="HH31" type="710"/>
   <Atom name="CH3" type="711"/>
   <Atom name="HH32" type="710"/>
   <Atom name="HH33" type="710"/>
   <Atom name="C" type="712"/>
   <Atom name="O" type="713"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="1" to="3"/>
   <Bond from="1" to="4"/>
   <Bond from="4" to="5"/>
   <ExternalBond from="4"/>
  </Residue>
 </Residues>
 <PeriodicTorsionForce>
  <Proper class1="C" class2="C" class3="C" class4="C" periodicity1="2" phase1="3.14159265359" k1="10.46"/>
  <Improper class1="C" class2="C" class3="C" class4="C" periodicity1="2" phase1="3.14159265359" k1="43.932"/>
 </PeriodicTorsionForce>
</ForceField>"""

        ff1 = ForceField(StringIO(ffxml), 'amber99sbildn.xml')
        ff2 = ForceField('amber99sbildn.xml', StringIO(ffxml))

        self.assertEqual(len(ff1._forces), 4)
        self.assertEqual(len(ff2._forces), 4)

        pertorsion1 = ff1._forces[0]
        pertorsion2 = ff2._forces[2]

        self.assertEqual(len(pertorsion1.proper), 110)
        self.assertEqual(len(pertorsion1.improper), 42)
        self.assertEqual(len(pertorsion2.proper), 110)
        self.assertEqual(len(pertorsion2.improper), 42)

    def test_ResidueTemplateUserChoice(self):
        """Test createSystem does not allow multiple matching templates, unless
           user has specified which template to use via residueTemplates arg"""
        ffxml = """<ForceField>
 <AtomTypes>
  <Type name="Fe2+" class="Fe2+" element="Fe" mass="55.85"/>
  <Type name="Fe3+" class="Fe3+" element="Fe" mass="55.85"/>
 </AtomTypes>
 <Residues>
  <Residue name="FE2">
   <Atom name="FE2" type="Fe2+" charge="2.0"/>
  </Residue>
  <Residue name="FE">
   <Atom name="FE" type="Fe3+" charge="3.0"/>
  </Residue>
 </Residues>
 <NonbondedForce coulomb14scale="0.833333333333" lj14scale="0.5">
  <UseAttributeFromResidue name="charge"/>
  <Atom type="Fe2+" sigma="0.227535532613" epsilon="0.0150312292"/>
  <Atom type="Fe3+" sigma="0.192790482606" epsilon="0.00046095128"/>
 </NonbondedForce>
</ForceField>"""

        pdb_string = "ATOM      1 FE    FE A   1      20.956  27.448 -29.067  1.00  0.00          Fe"
        ff = ForceField(StringIO(ffxml))
        pdb = PDBFile(StringIO(pdb_string))

        self.assertRaises(Exception, lambda: ff.createSystem(pdb.topology))
        sys = ff.createSystem(pdb.topology, residueTemplates={list(pdb.topology.residues())[0] : 'FE2'})
        # confirm charge
        self.assertEqual(sys.getForce(0).getParticleParameters(0)[0]._value, 2.0)
        sys = ff.createSystem(pdb.topology, residueTemplates={list(pdb.topology.residues())[0] : 'FE'})
        # confirm charge
        self.assertEqual(sys.getForce(0).getParticleParameters(0)[0]._value, 3.0)

    def test_ResidueOverriding(self):
        """Test residue overriding via override tag in the XML"""

        ffxml1 = """<ForceField>
 <AtomTypes>
  <Type name="Fe2+_tip3p_HFE" class="Fe2+_tip3p_HFE" element="Fe" mass="55.85"/>
 </AtomTypes>
 <Residues>
  <Residue name="FE2">
   <Atom name="FE2" type="Fe2+_tip3p_HFE" charge="2.0"/>
  </Residue>
 </Residues>
 <NonbondedForce coulomb14scale="0.833333333333" lj14scale="0.5">
  <UseAttributeFromResidue name="charge"/>
  <Atom type="Fe2+_tip3p_HFE" sigma="0.227535532613" epsilon="0.0150312292"/>
 </NonbondedForce>
</ForceField>"""

        ffxml2 = """<ForceField>
 <AtomTypes>
  <Type name="Fe2+_tip3p_standard" class="Fe2+_tip3p_standard" element="Fe" mass="55.85"/>
 </AtomTypes>
 <Residues>
  <Residue name="FE2">
   <Atom name="FE2" type="Fe2+_tip3p_standard" charge="2.0"/>
  </Residue>
 </Residues>
 <NonbondedForce coulomb14scale="0.833333333333" lj14scale="0.5">
  <UseAttributeFromResidue name="charge"/>
  <Atom type="Fe2+_tip3p_standard" sigma="0.241077193129" epsilon="0.03940482832"/>
 </NonbondedForce>
</ForceField>"""

        ffxml3 = """<ForceField>
 <AtomTypes>
  <Type name="Fe2+_tip3p_standard" class="Fe2+_tip3p_standard" element="Fe" mass="55.85"/>
 </AtomTypes>
 <Residues>
  <Residue name="FE2" override="1">
   <Atom name="FE2" type="Fe2+_tip3p_standard" charge="2.0"/>
  </Residue>
 </Residues>
 <NonbondedForce coulomb14scale="0.833333333333" lj14scale="0.5">
  <UseAttributeFromResidue name="charge"/>
  <Atom type="Fe2+_tip3p_standard" sigma="0.241077193129" epsilon="0.03940482832"/>
 </NonbondedForce>
</ForceField>"""

        pdb_string = "ATOM      1 FE    FE A   1      20.956  27.448 -29.067  1.00  0.00          Fe"
        pdb = PDBFile(StringIO(pdb_string))

        self.assertRaises(Exception, lambda: ForceField(StringIO(ffxml1), StringIO(ffxml2)))
        ff = ForceField(StringIO(ffxml1), StringIO(ffxml3))
        self.assertEqual(ff._templates['FE2'].atoms[0].type, 'Fe2+_tip3p_standard')
        ff.createSystem(pdb.topology)

    def test_LennardJonesGenerator(self):
        """ Test the LennardJones generator"""
        warnings.filterwarnings('ignore', category=CharmmPSFWarning)
        psf = CharmmPsfFile('systems/ions.psf')
        pdb = PDBFile('systems/ions.pdb')
        params = CharmmParameterSet('systems/toppar_water_ions.str'
                                    )

        # Box dimensions (found from bounding box)
        psf.setBox(12.009*angstroms,   12.338*angstroms,   11.510*angstroms)

        # Turn off charges so we only test the Lennard-Jones energies
        for a in psf.atom_list:
            a.charge = 0.0

        # Now compute the full energy
        plat = Platform.getPlatformByName('Reference')
        system = psf.createSystem(params, nonbondedMethod=PME,
                                  nonbondedCutoff=5*angstroms)

        con = Context(system, VerletIntegrator(2*femtoseconds), plat)
        con.setPositions(pdb.positions)

        # Now set up system from ffxml.
        xml = """
<ForceField>
 <AtomTypes>
  <Type name="SOD" class="SOD" element="Na" mass="22.98977"/>
  <Type name="CLA" class="CLA" element="Cl" mass="35.45"/>
 </AtomTypes>
 <Residues>
  <Residue name="CLA">
   <Atom name="CLA" type="CLA"/>
  </Residue>
  <Residue name="SOD">
   <Atom name="SOD" type="SOD"/>
  </Residue>
 </Residues>
 <LennardJonesForce lj14scale="1.0" useDispersionCorrection="False">
  <Atom type="CLA" sigma="0.404468018036" epsilon="0.6276"/>
  <Atom type="SOD" sigma="0.251367073323" epsilon="0.1962296"/>
  <NBFixPair type1="CLA" type2="SOD" sigma="0.33239431" epsilon="0.350933"/>
 </LennardJonesForce>
</ForceField> """
        ff = ForceField(StringIO(xml))
        system2 = ff.createSystem(pdb.topology, nonbondedMethod=PME,
                                  nonbondedCutoff=5*angstroms)
        con2 = Context(system2, VerletIntegrator(2*femtoseconds), plat)
        con2.setPositions(pdb.positions)

        state = con.getState(getEnergy=True, enforcePeriodicBox=True)
        ene = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
        state2 = con2.getState(getEnergy=True, enforcePeriodicBox=True)
        ene2 = state2.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
        self.assertAlmostEqual(ene, ene2)

    def test_NBFix(self):
        """Test using LennardJonesGenerator to implement NBFix terms."""
        # Create a chain of five atoms.

        top = Topology()
        chain = top.addChain()
        res = top.addResidue('RES', chain)
        top.addAtom('A', elem.oxygen, res)
        top.addAtom('B', elem.carbon, res)
        top.addAtom('C', elem.carbon, res)
        top.addAtom('D', elem.carbon, res)
        top.addAtom('E', elem.nitrogen, res)
        atoms = list(top.atoms())
        top.addBond(atoms[0], atoms[1])
        top.addBond(atoms[1], atoms[2])
        top.addBond(atoms[2], atoms[3])
        top.addBond(atoms[3], atoms[4])

        # Create the force field and system.

        xml = """
<ForceField>
 <AtomTypes>
  <Type name="A" class="A" element="O" mass="1"/>
  <Type name="B" class="B" element="C" mass="1"/>
  <Type name="C" class="C" element="C" mass="1"/>
  <Type name="D" class="D" element="C" mass="1"/>
  <Type name="E" class="E" element="N" mass="1"/>
 </AtomTypes>
 <Residues>
  <Residue name="RES">
   <Atom name="A" type="A"/>
   <Atom name="B" type="B"/>
   <Atom name="C" type="C"/>
   <Atom name="D" type="D"/>
   <Atom name="E" type="E"/>
   <Bond atomName1="A" atomName2="B"/>
   <Bond atomName1="B" atomName2="C"/>
   <Bond atomName1="C" atomName2="D"/>
   <Bond atomName1="D" atomName2="E"/>
  </Residue>
 </Residues>
 <LennardJonesForce lj14scale="0.3">
  <Atom type="A" sigma="1" epsilon="0.1"/>
  <Atom type="B" sigma="2" epsilon="0.2"/>
  <Atom type="C" sigma="3" epsilon="0.3"/>
  <Atom type="D" sigma="4" epsilon="0.4"/>
  <Atom type="E" sigma="4" epsilon="0.4"/>
  <NBFixPair type1="A" type2="D" sigma="2.5" epsilon="1.1"/>
  <NBFixPair type1="A" type2="E" sigma="3.5" epsilon="1.5"/>
 </LennardJonesForce>
</ForceField> """
        ff = ForceField(StringIO(xml))
        system = ff.createSystem(top)

        # Check that it produces the correct energy.

        integrator = VerletIntegrator(0.001)
        context = Context(system, integrator, Platform.getPlatform(0))
        positions = [Vec3(i, 0, 0) for i in range(5)]*nanometers
        context.setPositions(positions)
        def ljEnergy(sigma, epsilon, r):
            return 4*epsilon*((sigma/r)**12-(sigma/r)**6)
        expected = 0.3*ljEnergy(2.5, 1.1, 3) + 0.3*ljEnergy(3.0, sqrt(0.08), 3) + ljEnergy(3.5, 1.5, 4)
        self.assertAlmostEqual(expected, context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoules_per_mole))

    def test_IgnoreExternalBonds(self):
        """Test the ignoreExternalBonds option"""

        modeller = Modeller(self.pdb2.topology, self.pdb2.positions)
        modeller.delete([next(modeller.topology.residues())])
        self.assertRaises(Exception, lambda: self.forcefield2.createSystem(modeller.topology))
        system = self.forcefield2.createSystem(modeller.topology, ignoreExternalBonds=True)
        templates = self.forcefield2.getMatchingTemplates(modeller.topology, ignoreExternalBonds=True)
        self.assertEqual(2, len(templates))
        self.assertEqual('ALA', templates[0].name)
        self.assertEqual('NME', templates[1].name)

    def test_Includes(self):
        """Test using a ForceField that includes other files."""
        forcefield = ForceField(os.path.join('systems', 'ff_with_includes.xml'))
        self.assertTrue(len(forcefield._atomTypes) > 10)
        self.assertTrue('spce-O' in forcefield._atomTypes)
        self.assertTrue('HOH' in forcefield._templates)

    def test_ImpropersOrdering(self):
        """Test correctness of the ordering of atom indexes in improper torsions
        and the torsion.ordering parameter.
        """

        xml = """
<ForceField>
 <PeriodicTorsionForce ordering="amber">
  <Improper class1="C" class2="" class3="O2" class4="O2" periodicity1="2" phase1="3.14159265359" k1="43.932"/>
 </PeriodicTorsionForce>
</ForceField>
"""
        pdb = PDBFile('systems/impropers_ordering_tetrapeptide.pdb')
        # ff1 uses default ordering of impropers, ff2 uses "amber" for the one
        # problematic improper
        ff1 = ForceField('amber99sbildn.xml')
        ff2 = ForceField(StringIO(xml), 'amber99sbildn.xml')

        system1 = ff1.createSystem(pdb.topology)
        system2 = ff2.createSystem(pdb.topology)

        imp1 = system1.getForce(2).getTorsionParameters(158)
        imp2 = system2.getForce(0).getTorsionParameters(158)

        system1_indexes = [imp1[0], imp1[1], imp1[2], imp1[3]]
        system2_indexes = [imp2[0], imp2[1], imp2[2], imp2[3]]

        self.assertEqual(system1_indexes, [51, 56, 54, 55])
        self.assertEqual(system2_indexes, [51, 55, 54, 56])

    def test_ImpropersOrdering_smirnoff(self):
        """Test correctness of the ordering of atom indexes in improper torsions
        and the torsion.ordering parameter when using the 'smirnoff' mode.
        """

        # SMIRNOFF parameters for formaldehyde
        xml = """
<ForceField>
  <AtomTypes>
    <Type name="[H]C(=O)[H]$C1#0" element="C" mass="12.01078" class="[H]C(=O)[H]$C1#0"/>
    <Type name="[H]C(=O)[H]$O1#1" element="O" mass="15.99943" class="[H]C(=O)[H]$O1#1"/>
    <Type name="[H]C(=O)[H]$H1#2" element="H" mass="1.007947" class="[H]C(=O)[H]$H1#2"/>
    <Type name="[H]C(=O)[H]$H2#3" element="H" mass="1.007947" class="[H]C(=O)[H]$H2#3"/>
  </AtomTypes>
  <PeriodicTorsionForce ordering="smirnoff">
    <Improper class1="[H]C(=O)[H]$C1#0" class2="[H]C(=O)[H]$O1#1" class3="[H]C(=O)[H]$H1#2" class4="[H]C(=O)[H]$H2#3" periodicity1="2" phase1="3.141592653589793" k1="1.5341333333333336"/>
    <Improper class1="[H]C(=O)[H]$C1#0" class2="[H]C(=O)[H]$H1#2" class3="[H]C(=O)[H]$H2#3" class4="[H]C(=O)[H]$O1#1" periodicity1="2" phase1="3.141592653589793" k1="1.5341333333333336"/>
    <Improper class1="[H]C(=O)[H]$C1#0" class2="[H]C(=O)[H]$H2#3" class3="[H]C(=O)[H]$O1#1" class4="[H]C(=O)[H]$H1#2" periodicity1="2" phase1="3.141592653589793" k1="1.5341333333333336"/>
  </PeriodicTorsionForce>
  <Residues>
    <Residue name="[H]C(=O)[H]">
      <Atom name="C1" type="[H]C(=O)[H]$C1#0" charge="0.5632799863815308"/>
      <Atom name="O1" type="[H]C(=O)[H]$O1#1" charge="-0.514739990234375"/>
      <Atom name="H1" type="[H]C(=O)[H]$H1#2" charge="-0.02426999807357788"/>
      <Atom name="H2" type="[H]C(=O)[H]$H2#3" charge="-0.02426999807357788"/>
      <Bond atomName1="C1" atomName2="O1"/>
      <Bond atomName1="C1" atomName2="H1"/>
      <Bond atomName1="C1" atomName2="H2"/>
    </Residue>
  </Residues>
</ForceField>
"""
        pdb = PDBFile('systems/formaldehyde.pdb')
        # ff1 uses default ordering of impropers, ff2 uses "amber" for the one
        # problematic improper
        ff = ForceField(StringIO(xml))

        system = ff.createSystem(pdb.topology)

        # Check that impropers are applied in the correct three-fold trefoil pattern
        forces = { force.__class__.__name__ : force for force in system.getForces() }
        force = forces['PeriodicTorsionForce']
        created_torsions = set()
        for index in range(force.getNumTorsions()):
            i,j,k,l,_,_,_ = force.getTorsionParameters(index)
            created_torsions.add((i,j,k,l))
        expected_torsions = set([(0,3,1,2), (0,1,2,3), (0,2,3,1)])
        self.assertEqual(expected_torsions, created_torsions)

    def test_Disulfides(self):
        """Test that various force fields handle disulfides correctly."""
        pdb = PDBFile('systems/bpti.pdb')
        for ff in ['amber99sb.xml', 'amber14-all.xml', 'charmm36.xml', 'amberfb15.xml', 'amoeba2013.xml']:
            forcefield = ForceField(ff)
            system = forcefield.createSystem(pdb.topology)


class AmoebaTestForceField(unittest.TestCase):
    """Test the ForceField.createSystem() method with the AMOEBA forcefield."""

    def setUp(self):
        """Set up the tests by loading the input pdb files and force field
        xml files.

        """

        self.pdb1 = PDBFile('systems/amoeba-ion-in-water.pdb')
        self.forcefield1 = ForceField('amoeba2013.xml')
        self.topology1 = self.pdb1.topology


    def test_NonbondedMethod(self):
        """Test both options for the nonbondedMethod parameter."""

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

    def test_RigidWater(self):
        """Test that AMOEBA creates rigid water with the correct geometry."""

        system = self.forcefield1.createSystem(self.pdb1.topology, rigidWater=True)
        constraints = dict()
        for i in range(system.getNumConstraints()):
            p1,p2,dist = system.getConstraintParameters(i)
            if p1 < 3:
                constraints[(min(p1,p2), max(p1,p2))] = dist.value_in_unit(nanometers)
        hoDist = 0.09572
        hohAngle = 108.50*math.pi/180.0
        hohDist = math.sqrt(2*hoDist**2 - 2*hoDist**2*math.cos(hohAngle))
        self.assertAlmostEqual(constraints[(0,1)], hoDist)
        self.assertAlmostEqual(constraints[(0,2)], hoDist)
        self.assertAlmostEqual(constraints[(1,2)], hohDist)

    def test_Forces(self):
        """Compute forces and compare them to ones generated with a previous version of OpenMM to ensure they haven't changed."""

        pdb = PDBFile('systems/alanine-dipeptide-implicit.pdb')
        forcefield = ForceField('amoeba2013.xml', 'amoeba2013_gk.xml')
        system = forcefield.createSystem(pdb.topology, polarization='direct')
        integrator = VerletIntegrator(0.001)
        context = Context(system, integrator)
        context.setPositions(pdb.positions)
        state1 = context.getState(getForces=True)
        with open('systems/alanine-dipeptide-amoeba-forces.xml') as input:
            state2 = XmlSerializer.deserialize(input.read())
        for f1, f2, in zip(state1.getForces().value_in_unit(kilojoules_per_mole/nanometer), state2.getForces().value_in_unit(kilojoules_per_mole/nanometer)):
            diff = norm(f1-f2)
            self.assertTrue(diff < 0.1 or diff/norm(f1) < 1e-3)


if __name__ == '__main__':
    unittest.main()

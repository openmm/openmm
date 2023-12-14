.. default-domain:: py

.. _creating-force-fields:

Creating Force Fields
#####################

OpenMM uses a simple XML file format to describe force fields.  It includes many
common force fields, but you can also create your own.  A force field can use
all the standard OpenMM force classes, as well as the very flexible custom force
classes.  You can even extend the ForceField class to add support for completely
new forces, such as ones defined in plugins.  This makes it a powerful tool for
force field development.

Basic Concepts
**************

Let’s start by considering how OpenMM defines a force field.  There are a small
number of basic concepts to understand.

Atom Types and Atom Classes
===========================

Force field parameters are assigned to atoms based on their “atom types”.  Atom
types should be the most specific identification of an atom that will ever be
needed.  Two atoms should have the same type only if the force field will always
treat them identically in every way.

Multiple atom types can be grouped together into “atom classes”.  In general,
two types should be in the same class if the force field usually (but not
necessarily always) treats them identically.  For example, the :math:`\alpha`\ -carbon of an
alanine residue will probably have a different atom type than the :math:`\alpha`\ -carbon of a
leucine residue, but both of them will probably have the same atom class.

All force field parameters can be specified either by atom type or atom class.
Classes exist as a convenience to make force field definitions more compact.  If
necessary, you could define everything in terms of atom types, but when many
types all share the same parameters, it is convenient to only have to specify
them once.

Residue Templates
=================

Types are assigned to atoms by matching residues to templates.  A template
specifies a list of atoms, the type of each one, and the bonds between them.
For each residue in the PDB file, the force field searches its list of templates
for one that has an identical set of atoms with identical bonds between them.
When matching templates, neither the order of the atoms nor their names matter;
it only cares about their elements and the set of bonds between them.  (The PDB
file reader does care about names, of course, since it needs to figure out which
atom each line of the file corresponds to.)

Forces
======

Once a force field has defined its atom types and residue templates, it must
define its force field parameters.  This generally involves one block of XML for
each Force object that will be added to the System.  The details are different
for each Force, but it generally consists of a set of rules for adding
interactions based on bonds and atom types or classes.  For example, when adding
a HarmonicBondForce, the force field will loop over every pair of bonded atoms,
check their types and classes, and see if they match any of its rules.  If so,
it will call :code:`addBond()` on the HarmonicBondForce.  If none of them
match, it simply ignores that pair and continues.

Writing the XML File
********************

The root element of the XML file must be a :code:`<ForceField>` tag:

.. code-block:: xml

    <ForceField>
    ...
    </ForceField>

The :code:`<ForceField>` tag contains the following children:

* An :code:`<AtomTypes>` tag containing the atom type definitions
* A :code:`<Residues>` tag containing the residue template definitions
* Optionally a :code:`<Patches>` tag containing patch definitions
* Zero or more tags defining specific forces


The order of these tags does not matter.  They are described in detail below.

<AtomTypes>
===========

The atom type definitions look like this:

.. code-block:: xml

    <AtomTypes>
     <Type name="0" class="N" element="N" mass="14.00672"/>
     <Type name="1" class="H" element="H" mass="1.007947"/>
     <Type name="2" class="CT" element="C" mass="12.01078"/>
     ...
    </AtomTypes>

There is one :code:`<Type>` tag for each atom type.  It specifies the name
of the type, the name of the class it belongs to, the symbol for its element,
and its mass in amu.  The names are arbitrary strings: they need not be numbers,
as in this example.  The only requirement is that all types have unique names.
The classes are also arbitrary strings, and in general will not be unique.  Two
types belong to the same class if they list the same value for the
:code:`class` attribute.

<Residues>
==========

The residue template definitions look like this:

.. code-block:: xml

    <Residues>
     <Residue name="ACE">
      <Atom name="HH31" type="710"/>
      <Atom name="CH3" type="711"/>
      <Atom name="HH32" type="710"/>
      <Atom name="HH33" type="710"/>
      <Atom name="C" type="712"/>
      <Atom name="O" type="713"/>
      <Bond atomName1="HH31" atomName2="CH3"/>
      <Bond atomName1="CH3" atomName2="HH32"/>
      <Bond atomName1="CH3" atomName2="HH33"/>
      <Bond atomName1="CH3" atomName2="C"/>
      <Bond atomName1="C" atomName2="O"/>
      <ExternalBond atomName="C"/>
     </Residue>
     <Residue name="ALA">
      ...
     </Residue>
     ...
    </Residues>

There is one :code:`<Residue>` tag for each residue template.  That in turn
contains the following tags:

* An :code:`<Atom>` tag for each atom in the residue.  This specifies the
  name of the atom and its atom type.
* A :code:`<Bond>` tag for each pair of atoms that are bonded to each
  other.  The :code:`atomName1` and :code:`atomName2` attributes are the names
  of the two bonded atoms.  (Some older force fields use the alternate tags
  :code:`to` and :code:`from` to specify the atoms by index instead of name.
  This is still supported for backward compatibility, but specifying atoms by
  name is recommended, since it makes the residue definition much easier to
  understand.)
* An :code:`<ExternalBond>` tag for each atom that will be bonded to an
  atom of a different residue.  :code:`atomName` is the name of the atom.
  (Alternatively, the deprecated :code:`from` tag may indicate the atom by
  index instead of name.)


The :code:`<Residue>` tag may also contain :code:`<VirtualSite>` tags,
as in the following example:

.. code-block:: xml

    <Residue name="HOH">
     <Atom name="O" type="tip4pew-O"/>
     <Atom name="H1" type="tip4pew-H"/>
     <Atom name="H2" type="tip4pew-H"/>
     <Atom name="M" type="tip4pew-M"/>
     <VirtualSite type="average3" siteName="M" atomName1="O" atomName2="H1" atomName3="H2"
         weight1="0.786646558" weight2="0.106676721" weight3="0.106676721"/>
     <Bond atomName1="O" atomName2="H1"/>
     <Bond atomName1="O" atomName2="H2"/>
    </Residue>

Each :code:`<VirtualSite>` tag indicates an atom in the residue that should
be represented with a virtual site.  The :code:`type` attribute may equal
:code:`"average2"`\ , :code:`"average3"`\ , :code:`"outOfPlane"`\ , or
:code:`"localCoords"`\ , which correspond to the TwoParticleAverageSite, ThreeParticleAverageSite,
OutOfPlaneSite, and LocalCoordinatesSite classes respectively.  See Section
:numref:`virtual-sites` for descriptions of how the virtual site classes work.  The :code:`siteName`
attribute gives the name of the atom to represent with a virtual site.  The atoms
it is calculated based on are specified by :code:`atomName1`\ , :code:`atomName2`\ , etc.
(Some old force fields use the deprecated tags :code:`index`, :code:`atom1`,
:code:`atom2`, etc. to refer to them by index instead of name.)

The remaining attributes are specific to the virtual site class, and specify the
parameters for calculating the site position.  For a TwoParticleAverageSite,
they are :code:`weight1` and :code:`weight2`\ .  For a
ThreeParticleAverageSite, they are :code:`weight1`\ , :code:`weight2`\ , and
\ :code:`weight3`\ . For an OutOfPlaneSite, they are :code:`weight12`\ ,
:code:`weight13`\ , and :code:`weightCross`\ . For a LocalCoordinatesSite, they
are :code:`p1`\ , :code:`p2`\ , and :code:`p3` (giving the x, y, and z coordinates
of the site position in the local coordinate system), and :code:`wo1`\ ,
:code:`wx1`\ , :code:`wy1`\ , :code:`wo2`\ , :code:`wx2`\ , :code:`wy2`\ , ...
(giving the weights for computing the origin, x axis, and y axis).

<Patches>
=========

A "patch" is a set of rules for modifying a residue template (or possibly multiple
templates at once).  For example a terminal amino acid is slightly different from
one in the middle of a chain.  A force field could of course define multiple
templates for each amino acid (standard, N-terminal, C-terminal, and monomer),
but since the modifications are the same for nearly all amino acids, it is simpler
to include only the "standard" templates, along with a set of patches for
modifying terminal residues.

Here is an example of a patch definition:

.. code-block:: xml

    <Patches>
     <Patch name="NTER">
      <RemoveAtom name="H"/>
      <RemoveBond atomName1="N" atomName2="H"/>
      <AddAtom name="H1" type="H"/>
      <AddAtom name="H2" type="H"/>
      <AddAtom name="H3" type="H"/>
      <AddBond atomName1="N" atomName2="H1"/>
      <AddBond atomName1="N" atomName2="H2"/>
      <AddBond atomName1="N" atomName2="H3"/>
      <RemoveExternalBond atomName="N"/>
      <ChangeAtom name="N" type="N3"/>
     </Patch>
    </Patches>

There is one :code:`<Patch>` tag for each patch definition.  That in turn may
contain any of the following tags:

 * An :code:`<AddAtom>` tag indicates that an atom should be added to the
   template.  It specifies the name of the atom and its atom type.
 * A :code:`<ChangeAtom>` tag indicates that the type of an atom already present
   in the template should be altered.  It specifies the name of the atom and its
   new atom type.
 * A :code:`<RemoveAtom>` tag indicates that an atom should be removed from the
   template.  It specifies the name of the atom to remove.
 * An :code:`<AddBond>` tag indicates that a bond should be added to the
   template.  It specifies the names of the two bonded atoms.
 * A :code:`<RemoveBond>` tag indicates that a bond already present in the
   template should be removed.  It specifies the names of the two bonded atoms.
 * An :code:`<AddExternalBond>` tag indicates that a new external bond should be
   added to the template.  It specifies the name of the bonded atom.
 * A :code:`<RemoveExternalBond>` tag indicates that an external bond aleady
   present in the template should be removed.  It specifies the name of the
   bonded atom.

In addition to defining the patches, you also must identify which residue
templates each patch can be applied to.  This can be done in two ways.  The more
common one is to have each template identify the patches that can be applied to
it.  This is done with an :code:`<AllowPatch>` tag:

.. code-block:: xml

    <Residue name="ALA">
     <AllowPatch name="CTER"/>
     <AllowPatch name="NTER"/>
     ...
    </Residue>

Alternatively, the patch can indicate which residues it may be applied to.  This
is done with an :code:`<ApplyToResidue>` tag:

.. code-block:: xml

    <Patch name="NTER">
     <ApplyToResidue name="ALA"/>
     <ApplyToResidue name="ARG"/>
     ...
    </Patch>

A patch can alter multiple templates at once.  This is useful for creating bonds
between molecules, and allows the atom types in one residue to depend on the
identity of the other residue it is bonded to.  To create a multi-residue patch,
added a :code:`residues` attribute to the :code:`<Patch>` tag specifying how many
residues that patch covers.  Then whenever you refer to an atom, prefix its name
with the index of the residue it belongs to:

.. code-block:: xml

  <Patch name="Disulfide" residues="2">
    <RemoveAtom name="1:HG"/>
    <RemoveAtom name="2:HG"/>
    <AddBond atomName1="1:SG" atomName2="2:SG"/>
    <ApplyToResidue name="1:CYS"/>
    <ApplyToResidue name="2:CYS"/>
  </Patch>

In this example, the patch modifies two residues of the same type, but that need
not always be true.  Each :code:`<ApplyToResidue>` tag therefore indicates which
one of the residue templates it modifies may be of the specified type.  Similarly,
if a residue template includes an :code:`<AcceptPatch>` tag for a multi-residue
patch, it must specify the name of the patch, followed by the index of the residue
within that patch:

.. code-block:: xml

    <AllowPatch name="Disulfide:1"/>


Missing residue templates
=========================

.. CAUTION::
   These features are experimental, and their API is subject to change.

You can use the :meth:`getUnmatchedResidues()` method to get a list of residues
in the provided :code:`topology` object that do not currently have a matching
residue template defined in the :class:`ForceField`.
::

    pdb = PDBFile('input.pdb')
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    unmatched_residues = forcefield.getUnmatchedResidues(topology)

This is useful for identifying issues with prepared systems, debugging issues
with residue template definitions, or identifying which additional residues need
to be parameterized.

As a convenience for parameterizing new residues, you can also get a list of
residues and empty residue templates using :meth:`generateTemplatesForUnmatchedResidues`
::

    pdb = PDBFile('input.pdb')
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    [templates, residues] = forcefield.generateTemplatesForUnmatchedResidues(topology)
    # Se the atom types
    for template in templates:
        for atom in template.atoms:
            atom.type = ... # set the atom types here
        # Register the template with the forcefield.
        forcefield.registerResidueTemplate(template)

If you find that templates seem to be incorrectly matched, another useful
function :meth:`getMatchingTemplates()` can help you identify which templates
are being matched:
::

    pdb = PDBFile('input.pdb')
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    templates = forcefield.getMatchingTemplates(topology)
    for (residue, template) in zip(pdb.topology.residues(), templates):
        print("Residue %d %s matched template %s" % (residue.id, residue.name, template.name))

<HarmonicBondForce>
===================

To add a HarmonicBondForce to the System, include a tag that looks like this:

.. code-block:: xml

    <HarmonicBondForce>
     <Bond class1="C" class2="C" length="0.1525" k="259408.0"/>
     <Bond class1="C" class2="CA" length="0.1409" k="392459.2"/>
     <Bond class1="C" class2="CB" length="0.1419" k="374049.6"/>
     ...
    </HarmonicBondForce>

Every :code:`<Bond>` tag defines a rule for creating harmonic bond
interactions between atoms.  Each tag may identify the atoms either by type
(using the attributes :code:`type1` and :code:`type2`\ ) or by class
(using the attributes :code:`class1` and :code:`class2`\ ).  For every
pair of bonded atoms, the force field searches for a rule whose atom types or
atom classes match the two atoms.  If it finds one, it calls
:code:`addBond()` on the HarmonicBondForce with the specified parameters.
Otherwise, it ignores that pair and continues.  :code:`length` is the
equilibrium bond length in nm, and :code:`k` is the spring constant in
kJ/mol/nm\ :sup:`2`\ .

<HarmonicAngleForce>
====================

To add a HarmonicAngleForce to the System, include a tag that looks like this:

.. code-block:: xml

    <HarmonicAngleForce>
     <Angle class1="C" class2="C" class3="O" angle="2.094" k="669.44"/>
     <Angle class1="C" class2="C" class3="OH" angle="2.094" k="669.44"/>
     <Angle class1="CA" class2="C" class3="CA" angle="2.094" k="527.184"/>
     ...
    </HarmonicAngleForce>

Every :code:`<Angle>` tag defines a rule for creating harmonic angle
interactions between triplets of atoms.  Each tag may identify the atoms either
by type (using the attributes :code:`type1`\ , :code:`type2`\ , ...) or by
class (using the attributes :code:`class1`\ , :code:`class2`\ , ...).  The
force field identifies every set of three atoms in the system where the first is
bonded to the second, and the second to the third.  For each one, it searches
for a rule whose atom types or atom classes match the three atoms.  If it finds
one, it calls :code:`addAngle()` on the HarmonicAngleForce with the
specified parameters.  Otherwise, it ignores that set and continues.
:code:`angle` is the equilibrium angle in radians, and :code:`k` is the
spring constant in kJ/mol/radian\ :sup:`2`\ .

<PeriodicTorsionForce>
======================

To add a PeriodicTorsionForce to the System, include a tag that looks like this:

.. code-block:: xml

    <PeriodicTorsionForce>
     <Proper class1="HC" class2="CT" class3="CT" class4="CT" periodicity1="3" phase1="0.0"
         k1="0.66944"/>
     <Proper class1="HC" class2="CT" class3="CT" class4="HC" periodicity1="3" phase1="0.0"
         k1="0.6276"/>
     ...
     <Improper class1="N" class2="C" class3="CT" class4="O" periodicity1="2"
         phase1="3.14159265359" k1="4.6024"/>
     <Improper class1="N" class2="C" class3="CT" class4="H" periodicity1="2"
         phase1="3.14159265359" k1="4.6024"/>
     ...
    </PeriodicTorsionForce>

Every child tag defines a rule for creating periodic torsion interactions
between sets of four atoms.  Each tag may identify the atoms either by type
(using the attributes :code:`type1`\ , :code:`type2`\ , ...) or by class
(using the attributes :code:`class1`\ , :code:`class2`\ , ...).

The force field recognizes two different types of torsions: proper and improper.
A proper torsion involves four atoms that are bonded in sequence: 1 to 2, 2 to
3, and 3 to 4.  An improper torsion involves a central atom and three others
that are bonded to it: atoms 2, 3, and 4 are all bonded to atom 1.  The force
field begins by identifying every set of atoms in the system of each of these
types. For each one, it searches for a rule whose atom types or atom classes
match the four atoms.  If it finds one, it calls :code:`addTorsion()` on the
PeriodicTorsionForce with the specified parameters.  Otherwise, it ignores that
set and continues.  :code:`periodicity1` is the periodicity of the torsion,
\ :code:`phase1` is the phase offset in radians, and :code:`k1` is the
force constant in kJ/mol.

Each torsion definition can specify multiple periodic torsion terms to add to
its atoms.  To add a second one, just add three more attributes:
:code:`periodicity2`\ , :code:`phase2`\ , and :code:`k2`\ .  You can have as
many terms as you want.  Here is an example of a rule that adds three torsion
terms to its atoms:

.. code-block:: xml

    <Proper class1="CT" class2="CT" class3="CT" class4="CT"
        periodicity1="3" phase1="0.0" k1="0.75312"
        periodicity2="2" phase2="3.14159265359" k2="1.046"
        periodicity3="1" phase3="3.14159265359" k3="0.8368"/>

You can also use wildcards when defining torsions.  To do this, simply leave the
type or class name for an atom empty.  That will cause it to match any atom.
For example, the following definition will match any sequence of atoms where the
second atom has class OS and the third has class P:

.. code-block:: xml

    <Proper class1="" class2="OS" class3="P" class4="" periodicity1="3" phase1="0.0" k1="1.046"/>

The :code:`<PeriodicTorsionForce>` tag also supports an optional
:code:`ordering` attribute to provide better compatibility with the way
impropers are assigned in different simulation packages:

 * :code:`ordering="default"` specifies the default behavior if the attribute
   is omitted.
 * :code:`ordering="amber"` produces behavior that replicates the behavior of
   AmberTools LEaP
 * :code:`ordering="charmm"` produces behavior more consistent with CHARMM
 * :code:`ordering="smirnoff"` allows multiple impropers to be added using
   exact matching to replicate the beheavior of `SMIRNOFF <https://open-forcefield-toolkit.readthedocs.io/en/latest/users/smirnoff.html>`_
   impropers

Different :code:`<PeriodicTorsionForce>` tags can specify different :code:`ordering`
values to be used for the sub-elements appearing within their tags.

<RBTorsionForce>
================

To add an RBTorsionForce to the System, include a tag that looks like this:

.. code-block:: xml

    <RBTorsionForce>
     <Proper class1="CT" class2="CT" class3="OS" class4="CT" c0="2.439272" c1="4.807416"
         c2="-0.8368" c3="-6.409888" c4="0" c5="0" />
     <Proper class1="C" class2="N" class3="CT" class4="C" c0="10.46" c1="-3.34720"
         c2="-7.1128" c3="0" c4="0" c5="0" />
     ...
     <Improper class1="N" class2="C" class3="CT" class4="O" c0="0.8368" c1="0"
         c2="-2.76144" c3="0" c4="3.3472" c5="0" />
     <Improper class1="N" class2="C" class3="CT" class4="H" c0="29.288" c1="-8.368"
         c2="-20.92" c3="0" c4="0" c5="0" />
     ...
    </RBTorsionForce>

Every child tag defines a rule for creating Ryckaert-Bellemans torsion
interactions between sets of four atoms.  Each tag may identify the atoms either
by type (using the attributes :code:`type1`\ , :code:`type2`\ , ...) or by
class (using the attributes :code:`class1`\ , :code:`class2`\ , ...).

The force field recognizes two different types of torsions: proper and improper.
A proper torsion involves four atoms that are bonded in sequence: 1 to 2, 2 to
3, and 3 to 4.  An improper torsion involves a central atom and three others
that are bonded to it: atoms 2, 3, and 4 are all bonded to atom 1.  The force
field begins by identifying every set of atoms in the system of each of these
types. For each one, it searches for a rule whose atom types or atom classes
match the four atoms.  If it finds one, it calls :code:`addTorsion()` on the
RBTorsionForce with the specified parameters.  Otherwise, it ignores that set
and continues.  The attributes :code:`c0` through :code:`c5` are the
coefficients of the terms in the Ryckaert-Bellemans force expression.

You can also use wildcards when defining torsions.  To do this, simply leave the
type or class name for an atom empty.  That will cause it to match any atom.
For example, the following definition will match any sequence of atoms where the
second atom has class OS and the third has class P:

.. code-block:: xml

    <Proper class1="" class2="OS" class3="P" class4="" c0="2.439272" c1="4.807416"
        c2="-0.8368" c3="-6.409888" c4="0" c5="0" />

<CMAPTorsionForce>
==================

To add a CMAPTorsionForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CMAPTorsionForce>
     <Map>
      0.0 0.809 0.951 0.309
      -0.587 -1.0 -0.587 0.309
      0.951 0.809 0.0 -0.809
      -0.951 -0.309 0.587 1.0
     </Map>
     <Torsion map="0" class1="CT" class2="CT" class3="C" class4="N" class5="CT"/>
     <Torsion map="0" class1="N" class2="CT" class3="C" class4="N" class5="CT"/>
     ...
    </CMAPTorsionForce>

Each :code:`<Map>` tag defines an energy correction map.  Its content is the
list of energy values in kJ/mole, listed in the correct order for
CMAPTorsionForce’s :code:`addMap()` method and separated by white space.
See the API documentation for details.  The size of the map is determined from
the number of energy values.

Each :code:`<Torsion>` tag defines a rule for creating CMAP torsion
interactions between sets of five atoms.  The tag may identify the atoms either
by type (using the attributes :code:`type1`\ , :code:`type2`\ , ...) or by
class (using the attributes :code:`class1`\ , :code:`class2`\ , ...).  The
force field identifies every set of five atoms that are bonded in sequence: 1 to
2, 2 to 3, 3 to 4, and 4 to 5.  For each one, it searches for a rule whose atom
types or atom classes match the five atoms.  If it finds one, it calls
:code:`addTorsion()` on the CMAPTorsionForce with the specified parameters.
Otherwise, it ignores that set and continues.  The first torsion is defined by
the sequence of atoms 1-2-3-4, and the second one by atoms 2-3-4-5.
:code:`map` is the index of the map to use, starting from 0, in the order they
are listed in the file.

You can also use wildcards when defining torsions.  To do this, simply leave the
type or class name for an atom empty.  That will cause it to match any atom.
For example, the following definition will match any sequence of five atoms
where the middle three have classes CT, C, and N respectively:

.. code-block:: xml

    <Torsion map="0" class1="" class2="CT" class3="C" class4="N" class5=""/>

<NonbondedForce>
================

To add a NonbondedForce to the System, include a tag that looks like this:

.. code-block:: xml

    <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
     <Atom type="0" charge="-0.4157" sigma="0.32499" epsilon="0.71128"/>
     <Atom type="1" charge="0.2719" sigma="0.10690" epsilon="0.06568"/>
     <Atom type="2" charge="0.0337" sigma="0.33996" epsilon="0.45772"/>
     ...
    </NonbondedForce>

The :code:`<NonbondedForce>` tag has two attributes
:code:`coulomb14scale` and :code:`lj14scale` that specify the scale
factors between pairs of atoms separated by three bonds.  After setting the
nonbonded parameters for all atoms, the force field calls
:code:`createExceptionsFromBonds()` on the NonbondedForce, passing in these
scale factors as arguments.

Each :code:`<Atom>` tag specifies the nonbonded parameters for one atom type
(specified with the :code:`type` attribute) or atom class (specified with
the :code:`class` attribute).  It is fine to mix these two methods, having
some tags specify a type and others specify a class.  However you do it, you
must make sure that a unique set of parameters is defined for every atom type.
:code:`charge` is measured in units of the proton charge, :code:`sigma`
is in nm, and :code:`epsilon` is in kJ/mole.

<GBSAOBCForce>
==============

To add a GBSAOBCForce to the System, include a tag that looks like this:

.. code-block:: xml

    <GBSAOBCForce>
     <Atom type="0" charge="-0.4157" radius="0.1706" scale="0.79"/>
     <Atom type="1" charge="0.2719" radius="0.115" scale="0.85"/>
     <Atom type="2" charge="0.0337" radius="0.19" scale="0.72"/>
     ...
    </GBSAOBCForce>

Each :code:`<Atom>` tag specifies the OBC parameters for one atom type
(specified with the :code:`type` attribute) or atom class (specified with
the :code:`class` attribute).  It is fine to mix these two methods, having
some tags specify a type and others specify a class.  However you do it, you
must make sure that a unique set of parameters is defined for every atom type.
:code:`charge` is measured in units of the proton charge, :code:`radius`
is the GBSA radius in nm, and :code:`scale` is the OBC scaling factor.

<CustomBondForce>
=================

To add a CustomBondForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CustomBondForce energy="scale*k*(r-r0)^2">
     <GlobalParameter name="scale" defaultValue="0.5"/>
     <PerBondParameter name="k"/>
     <PerBondParameter name="r0"/>
     <Bond class1="OW" class2="HW" r0="0.09572" k="462750.4"/>
     <Bond class1="HW" class2="HW" r0="0.15136" k="462750.4"/>
     <Bond class1="C" class2="C" r0="0.1525" k="259408.0"/>
     ...
    </CustomBondForce>

The energy expression for the CustomBondForce is specified by the
:code:`energy` attribute.  This is a mathematical expression that gives the
energy of each bond as a function of its length *r*\ .  It also may depend on
an arbitrary list of global or per-bond parameters.  Use a
:code:`<GlobalParameter>` tag to define a global parameter, and a
:code:`<PerBondParameter>` tag to define a per-bond parameter.

Every :code:`<Bond>` tag defines a rule for creating custom bond
interactions between atoms.  Each tag may identify the atoms either by type
(using the attributes :code:`type1` and :code:`type2`\ ) or by class
(using the attributes :code:`class1` and :code:`class2`\ ).  For every
pair of bonded atoms, the force field searches for a rule whose atom types or
atom classes match the two atoms.  If it finds one, it calls
:code:`addBond()` on the CustomBondForce.  Otherwise, it ignores that pair and
continues.  The remaining attributes are the values to use for the per-bond
parameters.  All per-bond parameters must be specified for every
:code:`<Bond>` tag, and the attribute name must match the name of the
parameter.  For instance, if there is a per-bond parameter with the name “k”,
then every :code:`<Bond>` tag must include an attribute called :code:`k`\ .

<CustomAngleForce>
==================

To add a CustomAngleForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CustomAngleForce energy="scale*k*(theta-theta0)^2">
     <GlobalParameter name="scale" defaultValue="0.5"/>
     <PerAngleParameter name="k"/>
     <PerAngleParameter name=" theta0"/>
     <Angle class1="HW" class2="OW" class3="HW" theta0="1.824218" k="836.8"/>
     <Angle class1="HW" class2="HW" class3="OW" theta0="2.229483" k="0.0"/>
     <Angle class1="C" class2="C" class3="O" theta0="2.094395" k="669.44"/>
     ...
    </CustomAngleForce>

The energy expression for the CustomAngleForce is specified by the
:code:`energy` attribute.  This is a mathematical expression that gives the
energy of each angle as a function of the angle *theta*\ .  It also may depend
on an arbitrary list of global or per-angle parameters.  Use a
:code:`<GlobalParameter>` tag to define a global parameter, and a
:code:`<PerAngleParameter>` tag to define a per-angle parameter.

Every :code:`<Angle>` tag defines a rule for creating custom angle
interactions between triplets of atoms.  Each tag may identify the atoms either
by type (using the attributes :code:`type1`\ , :code:`type2`\ , ...) or by
class (using the attributes :code:`class1`\ , :code:`class2`\ , ...).  The
force field identifies every set of three atoms in the system where the first is
bonded to the second, and the second to the third.  For each one, it searches
for a rule whose atom types or atom classes match the three atoms.  If it finds
one, it calls :code:`addAngle()` on the CustomAngleForce.  Otherwise, it
ignores that set and continues. The remaining attributes are the values to use
for the per-angle parameters. All per-angle parameters must be specified for
every :code:`<Angle>` tag, and the attribute name must match the name of the
parameter.  For instance, if there is a per-angle parameter with the name “k”,
then every :code:`<Angle>` tag must include an attribute called :code:`k`\ .

<CustomTorsionForce>
====================

To add a CustomTorsionForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CustomTorsionForce energy="scale*k*(1+cos(per*theta-phase))">
     <GlobalParameter name="scale" defaultValue="1"/>
     <PerTorsionParameter name="k"/>
     <PerTorsionParameter name="per"/>
     <PerTorsionParameter name="phase"/>
     <Proper class1="HC" class2="CT" class3="CT" class4="CT" per="3" phase="0.0" k="0.66944"/>
     <Proper class1="HC" class2="CT" class3="CT" class4="HC" per="3" phase="0.0" k="0.6276"/>
     ...
     <Improper class1="N" class2="C" class3="CT" class4="O" per="2" phase="3.14159265359"
         k="4.6024"/>
     <Improper class1="N" class2="C" class3="CT" class4="H" per="2" phase="3.14159265359"
         k="4.6024"/>
     ...
    </CustomTorsionForce>

The energy expression for the CustomTorsionForce is specified by the
:code:`energy` attribute.  This is a mathematical expression that gives the
energy of each torsion as a function of the angle *theta*\ .  It also may
depend on an arbitrary list of global or per-torsion parameters.  Use a
:code:`<GlobalParameter>` tag to define a global parameter, and a
:code:`<PerTorsionParameter>` tag to define a per-torsion parameter.

Every child tag defines a rule for creating custom torsion interactions between
sets of four atoms.  Each tag may identify the atoms either by type (using the
attributes :code:`type1`\ , :code:`type2`\ , ...) or by class (using the
attributes :code:`class1`\ , :code:`class2`\ , ...).

The force field recognizes two different types of torsions: proper and improper.
A proper torsion involves four atoms that are bonded in sequence: 1 to 2, 2 to
3, and 3 to 4.  An improper torsion involves a central atom and three others
that are bonded to it: atoms 2, 3, and 4 are all bonded to atom 1.  The force
field begins by identifying every set of atoms in the system of each of these
types. For each one, it searches for a rule whose atom types or atom classes
match the four atoms.  If it finds one, it calls :code:`addTorsion()` on the
CustomTorsionForce with the specified parameters.  Otherwise, it ignores that
set and continues. The remaining attributes are the values to use for the per-
torsion parameters.  Every :code:`<Torsion>` tag must include one attribute
for every per-torsion parameter, and the attribute name must match the name of
the parameter.

You can also use wildcards when defining torsions.  To do this, simply leave the
type or class name for an atom empty.  That will cause it to match any atom.
For example, the following definition will match any sequence of atoms where the
second atom has class OS and the third has class P:

.. code-block:: xml

    <Proper class1="" class2="OS" class3="P" class4="" per="3" phase="0.0" k="0.66944"/>

<CustomNonbondedForce>
======================

To add a CustomNonbondedForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CustomNonbondedForce energy="scale*epsilon1*epsilon2*((sigma1+sigma2)/r)^12" bondCutoff="3">
     <GlobalParameter name="scale" defaultValue="1"/>
     <PerParticleParameter name="sigma"/>
     <PerParticleParameter name="epsilon"/>
     <Atom type="0" sigma="0.3249" epsilon="0.7112"/>
     <Atom type="1" sigma="0.1069" epsilon="0.0656"/>
     <Atom type="2" sigma="0.3399" epsilon="0.4577"/>
     ...
    </CustomNonbondedForce>

The energy expression for the CustomNonbondedForce is specified by the
:code:`energy` attribute.  This is a mathematical expression that gives the
energy of each pairwise interaction as a function of the distance *r*\ .  It
also may depend on an arbitrary list of global or per-particle parameters.  Use
a :code:`<GlobalParameter>` tag to define a global parameter, and a
:code:`<PerParticleParameter>` tag to define a per-particle parameter.

The expression may also depend on computed values, each defined with a
:code:`<ComputedValue>` tag.  The tag should have two attributes, :code:`name`
with the name of the computed value, and :code:`expression` with the expression
used to compute it.

Exclusions are created automatically based on the :code:`bondCutoff` attribute.
After setting the nonbonded parameters for all atoms, the force field calls
:code:`createExclusionsFromBonds()` on the CustomNonbondedForce, passing in this
value as its argument.  To avoid creating exclusions, set :code:`bondCutoff` to 0.

Each :code:`<Atom>` tag specifies the parameters for one atom type
(specified with the :code:`type` attribute) or atom class (specified with
the :code:`class` attribute).  It is fine to mix these two methods, having
some tags specify a type and others specify a class.  However you do it, you
must make sure that a unique set of parameters is defined for every atom type.
The remaining attributes are the values to use for the per-atom parameters. All
per-atom parameters must be specified for every :code:`<Atom>` tag, and the
attribute name must match the name of the parameter.  For instance, if there is
a per-atom parameter with the name “radius”, then every :code:`<Atom>` tag
must include an attribute called :code:`radius`\ .

CustomNonbondedForce also allows you to define tabulated functions.  See Section
:numref:`tabulated-functions` for details.

<CustomGBForce>
===============

To add a CustomGBForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CustomGBForce>
     <GlobalParameter name="solventDielectric" defaultValue="78.3"/>
     <GlobalParameter name="soluteDielectric" defaultValue="1"/>
     <PerParticleParameter name="charge"/>
     <PerParticleParameter name="radius"/>
     <PerParticleParameter name="scale"/>
     <ComputedValue name="I" type="ParticlePairNoExclusions">
        step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(1/U^2-1/L^2)*(r-sr2*sr2/r)+0.5*log(L/U)/r+C);
        U=r+sr2; C=2*(1/or1-1/L)*step(sr2-r-or1); L=max(or1, D); D=abs(r-sr2); sr2 =
        scale2*or2; or1 = radius1-0.009; or2 = radius2-0.009
     </ComputedValue>
     <ComputedValue name="B" type="SingleParticle">
      1/(1/or-tanh(1*psi-0.8*psi^2+4.85*psi^3)/radius); psi=I*or; or=radius-0.009
     </ComputedValue>
     <EnergyTerm type="SingleParticle">
      28.3919551*(radius+0.14)^2*(radius/B)^6-0.5*138.935456*
              (1/soluteDielectric-1/solventDielectric)*charge^2/B
     </EnergyTerm>
     <EnergyTerm type="ParticlePair">
      -138.935456*(1/soluteDielectric-1/solventDielectric)*charge1*charge2/f;
              f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))
     </EnergyTerm>
     <Atom type="0" charge="-0.4157" radius="0.1706" scale="0.79"/>
     <Atom type="1" charge="0.2719" radius="0.115" scale="0.85"/>
     <Atom type="2" charge="0.0337" radius="0.19" scale="0.72"/>
     ...
    </CustomGBForce>

The above (rather complicated) example defines a generalized Born model that is
equivalent to GBSAOBCForce.  The definition consists of a set of computed values
(defined by :code:`<ComputedValue>` tags) and energy terms (defined by
:code:`<EnergyTerm>` tags), each of which is evaluated according to a
mathematical expression.  See the API documentation for details.

The expressions may depend on an arbitrary list of global or per-atom
parameters.  Use a :code:`<GlobalParameter>` tag to define a global
parameter, and a :code:`<PerAtomParameter>` tag to define a per-atom
parameter.

Each :code:`<Atom>` tag specifies the parameters for one atom type
(specified with the :code:`type` attribute) or atom class (specified with
the :code:`class` attribute).  It is fine to mix these two methods, having
some tags specify a type and others specify a class.  However you do it, you
must make sure that a unique set of parameters is defined for every atom type.
The remaining attributes are the values to use for the per-atom parameters. All
per-atom parameters must be specified for every :code:`<Atom>` tag, and the
attribute name must match the name of the parameter.  For instance, if there is
a per-atom parameter with the name “radius”, then every :code:`<Atom>` tag
must include an attribute called :code:`radius`\ .

CustomGBForce also allows you to define tabulated functions.  See Section
:numref:`tabulated-functions` for details.

<CustomHbondForce>
=========================

To add a CustomHbondForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CustomHbondForce particlesPerDonor="3" particlesPerAcceptor="2" bondCutoff="2"
        energy="scale*k*(distance(a1,d1)-r0)^2*(angle(a1,d1,d2)-theta0)^2">
     <GlobalParameter name="scale" defaultValue="1"/>
     <PerDonorParameter name="theta0"/>
     <PerAcceptorParameter name="k"/>
     <PerAcceptorParameter name="r0"/>
     <Donor class1="H" class2="N" class3="C" theta0="2.1"/>
     <Acceptor class1="O" class2="C" k="115.0" r0="0.2"/>
     ...
    </CustomHbondForce>

The energy expression for the CustomHbondForce is specified by the
:code:`energy` attribute.  This is a mathematical expression that gives the
energy of each donor-acceptor interaction as a function of various particle coordinates,
distances, and angles.  See the API documentation for details.  :code:`particlesPerDonor`
specifies the number of particles that make up a donor group, and :code:`particlesPerAcceptor`
specifies the number of particles that make up an acceptor group.

The expression may depend on an arbitrary list of global, per-donor, or
per-acceptor parameters.  Use a :code:`<GlobalParameter>` tag to define a global
parameter, a :code:`<PerDonorParameter>` tag to define a per-donor parameter,
and a :code:`<PerAcceptorParameter>` tag to define a per-acceptor parameter.

Exclusions are created automatically based on the :code:`bondCutoff` attribute.
If any atom of a donor is within the specified distance (measured in bonds) of
any atom of an acceptor, an exclusion is added to prevent them from interacting
with each other.  If a donor and an acceptor share any atom in common, that is a
bond distance of 0, so they are always excluded.

Every :code:`<Donor>` or :code:`<Acceptor>` tag defines a rule for creating donor
or acceptor groups.  The number of atoms specified in each one must match the
value of :code:`particlesPerDonor` or :code:`particlesPerAcceptor` specified in the
parent tag. Each tag may identify the atoms either by type (using the attributes
:code:`type1`\ , :code:`type2`\ , ...) or by class (using the attributes
:code:`class1`\ , :code:`class2`\ , ...).  The force field considers every atom
in the system (if the number of atoms is 1), every pair of bonded atoms (if the number
of atoms is 2), or every set of three atoms where the first is bonded to the second
and the second to the third (if the number of atoms is 3).  For each one, it searches
for a rule whose atom types or atom classes match the atoms.  If it finds one,
it calls :code:`addDonor()` or :code:`addAcceptor()` on the CustomHbondForce.
Otherwise, it ignores that set and continues. The remaining attributes are the
values to use for the per-donor and per-acceptor parameters. All parameters must
be specified for every tag, and the attribute name must match the name of the
parameter.  For instance, if there is a per-donor parameter with the name “k”,
then every :code:`<Donor>` tag must include an attribute called :code:`k`\ .

CustomHbondForce also allows you to define tabulated functions.  See Section
:numref:`tabulated-functions` for details.

<CustomManyParticleForce>
=========================

To add a CustomManyParticleForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CustomManyParticleForce particlesPerSet="3" permutationMode="UniqueCentralParticle"
        bondCutoff="3" energy="scale*(distance(p1,p2)-r1)*(distance(p1,p3)-r1)">
     <GlobalParameter name="scale" defaultValue="1"/>
     <PerParticleParameter name="r"/>
     <TypeFilter index="0" types="1,2"/>
     <Atom type="0" r="0.31" filterType="0"/>
     <Atom type="1" r="0.25" filterType="0"/>
     <Atom type="2" r="0.33" filterType="1"/>
     ...
    </CustomManyParticleForce>

The energy expression for the CustomManyParticleForce is specified by the
:code:`energy` attribute.  This is a mathematical expression that gives the
energy of each interaction as a function of various particle coordinates,
distances, and angles.  See the API documentation for details.  :code:`particlesPerSet`
specifies the number of particles involved in the interaction and
:code:`permutationMode` specifies the permutation mode.

The expression may depend on an arbitrary list of global or per-atom
parameters.  Use a :code:`<GlobalParameter>` tag to define a global
parameter, and a :code:`<PerAtomParameter>` tag to define a per-atom
parameter.

Exclusions are created automatically based on the :code:`bondCutoff` attribute.
After setting the nonbonded parameters for all atoms, the force field calls
:code:`createExclusionsFromBonds()` on the CustomManyParticleForce, passing in this
value as its argument.  To avoid creating exclusions, set :code:`bondCutoff` to 0.

Type filters may be specified with a :code:`<TypeFilter>` tag.  The :code:`index`
attribute specifies the index of the particle to apply the filter to, and
:code:`types` is a comma separated list of allowed types.

Each :code:`<Atom>` tag specifies the parameters for one atom type
(specified with the :code:`type` attribute) or atom class (specified with
the :code:`class` attribute).  It is fine to mix these two methods, having
some tags specify a type and others specify a class.  However you do it, you
must make sure that a unique set of parameters is defined for every atom type.
In addition, each :code:`<Atom>` tag must include the :code:`filterType`
attribute, which specifies the atom type for use in type filters.
The remaining attributes are the values to use for the per-atom parameters. All
per-atom parameters must be specified for every :code:`<Atom>` tag, and the
attribute name must match the name of the parameter.  For instance, if there is
a per-atom parameter with the name “radius”, then every :code:`<Atom>` tag
must include an attribute called :code:`radius`\ .

CustomManyParticleForce also allows you to define tabulated functions.  See Section
:numref:`tabulated-functions` for details.

<LennardJonesForce>
===================

The :code:`<LennardJonesForce>` tag provides an alternative to :code:`<NonbondedForce>`
for implementing Lennard-Jones nonbonded interactions.  It instead implements
them with a CustomNonbondedForce, and if necessary a CustomBondForce for 1-4
interactions.  The advantage of using this tag is that it provides more flexibility
in how the sigma and epsilon parameters are determined.  The disadvantage is
that it tends to run slightly slower, so when the extra flexibility is not
needed, it is better to use :code:`<NonbondedForce>` intead.

To use it, include a tag that looks like this:

.. code-block:: xml

    <LennardJonesForce lj14scale="1.0" useDispersionCorrection="True">
     <Atom epsilon="0.192464" sigma="0.040001" type="H"/>
     <Atom epsilon="0.192464" sigma="0.040001" type="HC"/>
     <Atom epsilon="0.092048" sigma="0.235197" type="HA"/>
     ...
     <NBFixPair epsilon="0.134306" sigma="0.29845" type1="CRL1" type2="HAL2"/>
     <NBFixPair epsilon="0.150205" sigma="0.23876" type1="HAL2" type2="HGA1"/>
     ...
    </LennardJonesForce>

The :code:`<LennardJonesForce>` tag has two attributes: :code:`lj14scale` specifies
the scale factor between pairs of atoms separated by three bonds, and
:code:`useDispersionCorrection` specifies whether to include a long range
dispersion correction.

Each :code:`<Atom>` tag specifies the nonbonded parameters for one atom type
(specified with the :code:`type` attribute) or atom class (specified with
the :code:`class` attribute).  It is fine to mix these two methods, having
some tags specify a type and others specify a class.  However you do it, you
must make sure that a unique set of parameters is defined for every atom type.
:code:`sigma` is in nm, and :code:`epsilon` is in kJ/mole.

An :code:`<Atom>` tag can optionally include two other attributes: :code:`sigma14`
and :code:`epsilon14`.  If they are present, they are used in place of the
standard sigma and epsilon parameters when computing 1-4 interactions.

Each :code:`<NBFixPair>` tag specifies a pair of atom types (:code:`type1` and
:code:`type2`) or classes (:code:`class1` and :code:`class2`) whose interaction
should be computed differently.  Instead of using the standard Lorentz-Berthelot
combining rules to determine sigma and epsilon based on the parameters for the
individual atoms, the tag specifies different values to use.

Because this tag only computes Lennard-Jones interactions, it usually is used together
with a :code:`<NonbondedForce>` to compute the Coulomb interactions.  In that
case, the NonbondedForce should specify :code:`epsilon="0"` for every atom so it
does not also compute the Lennard-Jones interactions.

Writing Custom Expressions
==========================

The custom forces described in this chapter involve user defined algebraic
expressions.  These expressions are specified as character strings, and may
involve a variety of standard operators and mathematical functions.

The following operators are supported: + (add), - (subtract), * (multiply), /
(divide), and ^ (power).  Parentheses “(“and “)” may be used for grouping.

The following standard functions are supported: sqrt, exp, log, sin, cos, sec,
csc, tan, cot, asin, acos, atan, sinh, cosh, tanh, erf, erfc, min, max, abs,
floor, ceil, step, delta, select. step(x) = 0 if x < 0, 1 otherwise.
delta(x) = 1 if x is 0, 0 otherwise.  select(x,y,z) = z if x = 0, y otherwise.
Some custom forces allow additional functions to be defined from tabulated values.

Numbers may be given in either decimal or exponential form.  All of the
following are valid numbers: 5, -3.1, 1e6, and 3.12e-2.

The variables that may appear in expressions are specified in the API
documentation for each force class.  In addition, an expression may be followed
by definitions for intermediate values that appear in the expression.  A
semicolon “;” is used as a delimiter between value definitions.  For example,
the expression
::

    a^2+a*b+b^2; a=a1+a2; b=b1+b2

is exactly equivalent to
::

    (a1+a2)^2+(a1+a2)*(b1+b2)+(b1+b2)^2

The definition of an intermediate value may itself involve other intermediate
values.  All uses of a value must appear *before* that value’s definition.

.. _tabulated-functions:

Tabulated Functions
===================

Some forces, such as CustomNonbondedForce and CustomGBForce, allow you to define
tabulated functions.  To define a function, include a :code:`<Function>` tag inside the
:code:`<CustomNonbondedForce>` or :code:`<CustomGBForce>` tag:

.. code-block:: xml

    <Function name="myfn" type="Continuous1D" min="-5" max="5">
    0.983674857694 -0.980096396266 -0.975743130031 -0.970451936613 -0.964027580076
    -0.956237458128 -0.946806012846 -0.935409070603 -0.921668554406 -0.905148253645
    -0.885351648202 -0.861723159313 -0.833654607012 -0.800499021761 -0.761594155956
    -0.716297870199 -0.664036770268 -0.604367777117 -0.537049566998 -0.46211715726
    -0.379948962255 -0.291312612452 -0.197375320225 -0.099667994625 0.0
    0.099667994625 0.197375320225 0.291312612452 0.379948962255 0.46211715726
    0.537049566998 0.604367777117 0.664036770268 0.716297870199 0.761594155956
    0.800499021761 0.833654607012 0.861723159313 0.885351648202 0.905148253645
    0.921668554406 0.935409070603 0.946806012846 0.956237458128 0.964027580076
    0.970451936613 0.975743130031 0.980096396266 0.983674857694 0.986614298151
    0.989027402201
    </Function>

The tag’s attributes define the name of the function, the type of function, and
the range of values for which it is defined.  The required set of attributed
depends on the function type:

.. tabularcolumns:: |l|L|

============  =======================================================
Type          Required Attributes
============  =======================================================
Continuous1D  min, max
Continuous2D  xmin, ymin, xmax, ymax, xsize, ysize
Continuous3D  xmin, ymin, zmin, xmax, ymax, zmax, xsize, ysize, zsize
Discrete1D
Discrete2D    xsize, ysize
Discrete3D    xsize, ysize, zsize
============  =======================================================


The "min" and "max" attributes define the range of the independent variables for
a continuous function.  The "size" attributes define the size of the table along
each axis.  The tabulated values are listed inside the body of the tag, with
successive values separated by white space.  See the API documentation for more
details.


Residue Template Parameters
===========================

In forces that use an :code:`<Atom>` tag to define parameters for atom types or
classes, there is an alternate mechanism you can also use: defining those
parameter values in the residue template.  This is useful for situations that
come up in certain force fields.  For example, :code:`NonbondedForce` and
:code:`GBSAOBCForce` each have a :code:`charge` attribute.  If you only have to
define the charge of each atom type once, that is more convenient and avoids
potential bugs.  Also, many force fields have a different charge for each atom
type, but Lennard-Jones parameters that are the same for all types in a class.
It would be preferable not to have to repeat those parameter values many times
over.

When writing a residue template, you can add arbitrary additional attributes
to each :code:`<Atom>` tag.  For example, you might include a :code:`charge`
attribute as follows:

.. code-block:: xml

   <Atom name="CA" type="53" charge="0.0381"/>

When writing the tag for a force, you can then include a
:code:`<UseAttributeFromResidue>` tag inside it.  This indicates that a
specified attribute should be taken from the residue template.  Finally, you
simply omit that attribute in the force's own :code:`<Atom>` tags.  For example:

.. code-block:: xml

    <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
     <UseAttributeFromResidue name="charge"/>
     <Atom class="CX" sigma="0.339966950842" epsilon="0.4577296"/>
     ...
    </NonbondedForce>

Notice that the :code:`charge` attribute is missing, and that the parameters
are specified by class, not by type.  This means that sigma and epsilon only
need to be specified once for each class.  The atom charges, which are different
for each type, are taken from the residue template instead.


Including Other Files
=====================

Sometimes it is useful to split a force field definition into multiple files,
but still be able to use the force field by specifying only a single file.  You
can accomplish this with the :code:`<Include>` tag.  For example:

.. code-block:: xml

    <ForceField>
     <Include file="definitions.xml"/>
     ...
    </ForceField>

The :code:`file` attribute gives the path of the file to include.  It may be
relative either to the directory containing the parent XML file (the one with
the :code:`<Include>` tag) or the OpenMM data directory (the one containing
built in force fields).


Using Multiple Files
********************

If multiple XML files are specified when a ForceField is created, their
definitions must be combined into a single force field.  This process involves
a few steps.

The first step is to load the atom type definitions from all files.  Because this
happens before any residue template or force definitions are loaded, it is
possible for the residue templates and forces in one file to refer to atom types
or classes defined in a different file.

Next the residue and patch definitions are loaded, and finally the forces.
Special care is needed when two files each define a force of the same type.

For standard forces like HarmonicBondForce or NonbondedForce, the definitions
are combined into a single definition.  For example, if two files each contain
a :code:`<HarmonicBondForce>` tag, only a single HarmonicBondForce will be added
to the system, containing bonds based on the definitions from both files.  An
important consequence is that only a single bond will be added for any pair of
atoms.  If each file contains a :code:`<Bond>` tag that matches the same pair of
atoms, and if those tags specify different parameters, no guarantee is made about
which parameters will end up being used.  It is generally best to avoid doing
this.

Combining of standard forces is especially important for forces that involve
per-atom parameters, such as NonbondedForce or GBSAOBCForce.  These forces
require parameter values to be defined for every atom type.  Because the
definitions are combined, it does not matter which file each type is defined in.
For example, files that define explicit water models generally define a small
number of atom types, as well as nonbonded parameters for those types.  They are
combined with the parameters from the main force field to create a single
NonbondedForce.  In contrast, files that define implicit solvent models do not
define any new atom types, but provide parameters for the atom types that were
defined in the main force field file.

In contrast to standard forces, custom forces do not get combined since two
objects of the same class may still represent different potential functions.
For example, if there are two :code:`<CustomBondForce>` tags, two separate
CustomBondForces are added to the System, each containing bonds based only on
the definitions in that tag.

A consequence is that custom forces requiring per-atom parameters, such as
CustomNonbondedForce or CustomGBForce, cannot be split between files.  A single
tag must define parameters for all atom types.  There is no way for a second
file to define parameters for additional atom types.


Extending ForceField
********************

The ForceField class is designed to be modular and extensible.  This means you
can add support for entirely new force types, such as ones implemented with
plugins.

Adding new force types
======================

For every force class, there is a “generator” class that parses the
corresponding XML tag, then creates Force objects and adds them to the System.
ForceField maintains a map of tag names to generator classes.  When a ForceField
is created, it scans through the XML files, looks up the generator class for
each tag, and asks that class to create a generator object based on it.  Then,
when you call :code:`createSystem()`\ ,  it loops over each of its generators
and asks each one to create its Force object.  Adding a new Force type therefore
is simply a matter of creating a new generator class and adding it to
ForceField’s map.

The generator class must define two methods.  First, it needs a static method
with the following signature to parse the XML tag and create the generator:
::

    @staticmethod
    def parseElement(element, forcefield):

:code:`element` is the XML tag (an xml.etree.ElementTree.Element object) and
:code:`forcefield` is the ForceField being created.  This method should
create a generator and add it to the ForceField:
::

    generator = MyForceGenerator()
    forcefield._forces.append(generator)

It then should parse the information contained in the XML tag and configure the
generator based on it.

Second, it must define a method with the following signature:
::

    def createForce(self, system, data, nonbondedMethod, nonbondedCutoff, args):

When :code:`createSystem()` is called on the ForceField, it first creates
the System object, then loops over each of its generators and calls
:code:`createForce()` on each one.  This method should create the Force object
and add it to the System.  :code:`data` is a ForceField._SystemData object
containing information about the System being created (atom types, bonds,
angles, etc.), :code:`system` is the System object, and the remaining
arguments are values that were passed to :code:`createSystem()`\ .  To get a
better idea of how this works, look at the existing generator classes in
forcefield.py.

The generator class may optionally also define a method with the following
signature:
::

    def postprocessSystem(self, system, data, args):

If this method exists, it will be called after all Forces have been created.
This gives generators a chance to make additional changes to the System.

Finally, you need to register your class by adding it to ForceField’s map:
::

    forcefield.parsers['MyForce'] = MyForceGenerator.parseElement

The key is the XML tag name, and the value is the static method to use for
parsing it.

Now you can simply create a ForceField object as usual.  If an XML file contains
a :code:`<MyForce>` tag, it will be recognized and processed correctly.

Adding residue template generators
==================================

.. CAUTION::
   This feature is experimental, and its API is subject to change.

Typically, when :class:`ForceField` encounters a residue it does not have a template for,
it simply raises an :code:`Exception`, since it does not know how to assign atom types for
the unknown residue.

However, :class:`ForceField` has an API for registering *residue template generators* that are
called when a residue without an existing template is encountered.  These generators
may create new residue templates that match existing atom types and parameters, or can
even create new atom types and new parameters that are added to :class:`ForceField`. This
functionality can be useful for adding residue template generators that are able to
parameterize small molecules that are not represented in a protein or nucleic acid
forcefield, for example, or for creating new residue templates for post-translationally
modified residues, covalently-bound ligands, or unnatural amino acids or bases.

To register a new residue template generator named :code:`generator`, simply call the
:meth:`registerTemplateGenerator` method on an existing :class:`ForceField` object:
::

    forcefield.registerTemplateGenerator(generator)

This :code:`generator` function must conform to the following API:
::

    def generator(forcefield, residue):
        """
        Parameters
        ----------
        forcefield : openmm.app.ForceField
            The ForceField object to which residue templates and/or parameters are to be added.
        residue : openmm.app.Topology.Residue
            The residue topology for which a template is to be generated.

        Returns
        -------
        success : bool
            If the generator is able to successfully parameterize the residue, `True` is returned.
            If the generator cannot parameterize the residue, it should return `False` and not
            modify `forcefield`.

        The generator should either register a residue template directly with
        `forcefield.registerResidueTemplate(template)` or it should call `forcefield.loadFile(file)`
        to load residue definitions from an ffxml file.

        It can also use the `ForceField` programmatic API to add additional atom types (via
        `forcefield.registerAtomType(parameters)`) or additional parameters.

        """

The :class:`ForceField` object will be modified by the residue template generator as residues without previously
defined templates are encountered.  Because these templates are added to the :class:`ForceField` as new residue
types are encountered, subsequent residues will be parameterized using the same residue templates without
calling the :code:`generator` again.

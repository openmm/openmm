.. default-domain:: py

.. _model-building-and-editing:

Model Building and Editing
##########################

Sometimes you have a PDB file that needs some work before you can simulate it.
Maybe it doesn’t contain hydrogen atoms (which is common for structures
determined by X-ray crystallography), so you need to add them.  Or perhaps you
want to simulate the system in explicit water, but the PDB file doesn’t contain
water molecules.  Or maybe it does contain water molecules, but they contain the
wrong number of interaction sites for the water model you want to use.  OpenMM’s
Modeller class can fix problems such as these.

To use it, create a :class:`Modeller` object, providing the initial :class:`Topology` and atom
positions.  You then can invoke various modelling functions on it.  Each one
modifies the system in some way, creating a new :class:`Topology` and list of positions.
When you are all done, you can retrieve them from the :class:`Modeller` and use them as
the starting point for your simulation:

.. samepage::
    ::

        ...
        pdb = PDBFile('input.pdb')
        modeller = Modeller(pdb.topology, pdb.positions)
        # ... Call some modelling functions here ...
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

    .. caption::

        :autonumber:`Example,Modeller outline`

Now let’s consider the particular functions you can call.

Adding Hydrogens
****************

Call the :meth:`addHydrogens` function to add missing hydrogen atoms:
::

    modeller.addHydrogens(forcefield)

The force field is needed to determine the positions for the hydrogen atoms.  If
the system already contains some hydrogens but is missing others, that is fine.
The Modeller will recognize the existing ones and figure out which ones need to
be added.

Some residues can exist in different protonation states depending on the pH and
on details of the local environment.  By default it assumes pH 7, but you can
specify a different value:
::

    modeller.addHydrogens(forcefield, pH=5.0)

For each residue, it selects the protonation state that is most common at the
specified pH.  In the case of Cysteine residues, it also checks whether the
residue participates in a disulfide bond when selecting the state to use.
Histidine has two different protonation states that are equally likely at
neutral pH.  It therefore selects which one to use based on which will form a
better hydrogen bond.

If you want more control, it is possible to specify exactly which protonation
state to use for particular residues.  For details, consult the API
documentation for the Modeller class.

By default, :class:`Modeller` loads information about hydrogens in standard
amino acids and nucleic acids.  You can call :meth:`loadHydrogenDefinitions` to
load definitions for other types of molecules.  In particular, if your system
contains carbohydrates that you plan to simulate with the GLYCAM force field, call
::

    Modeller.loadHydrogenDefinitions('glycam-hydrogens.xml')

All subsequent calls to :meth:`addHydrogens` will make use of the newly loaded
definitions.

Adding Solvent
**************

Call :meth:`addSolvent` to create a box of solvent (water and ions) around the model:
::

    modeller.addSolvent(forcefield)

This constructs a box of water around the solute, ensuring that no water
molecule comes closer to any solute atom than the sum of their van der Waals
radii.  It also determines the charge of the solute, and adds enough positive or
negative ions to make the system neutral.

When called as shown above, :meth:`addSolvent` expects that periodic box dimensions were
specified in the PDB file, and it uses them as the size for the water box.  If
your PDB file does not specify a box size, or if you want to use a different
size, you can specify one:
::

    modeller.addSolvent(forcefield, boxSize=Vec3(5.0, 3.5, 3.5)*nanometers)

This requests a 5 nm by 3.5 nm by 3.5 nm box.  For a non-rectangular box, you
can specify the three box vectors defining the unit cell:
::

    modeller.addSolvent(forcefield, boxVectors=(avec, bvec, cvec))

Another option is to specify a padding distance:
::

    modeller.addSolvent(forcefield, padding=1.0*nanometers)

This determines the largest size of the solute along any axis (x, y, or z).  It
then creates a cubic box of width (solute size)+2*(padding).  The above line
guarantees that no part of the solute comes closer than 1 nm to any edge of the
box.

Finally, you can specify the exact number of solvent molecules (including both
water and ions) to add.  This is useful when you want to solvate several different
conformations of the same molecule while guaranteeing they all have the same
amount of solvent:
::

    modeller.addSolvent(forcefield, numAdded=5000)

By default, :meth:`addSolvent` creates TIP3P water molecules, but it also supports other
water models:
::

    modeller.addSolvent(forcefield, model='tip5p')

Allowed values for the :code:`model` option are ``'tip3p'``, ``'spce'``,
``'tip4pew'``, ``'tip5p'``, and ``'swm4ndp'``.  Be sure to include the single quotes
around the value.

Another option is to add extra ion pairs to give a desired total ionic strength.
For example:
::

    modeller.addSolvent(forcefield, ionicStrength=0.1*molar)

This solvates the system with a salt solution whose ionic strength is 0.1 molar.
Note that when computing the ionic strength, it does *not* consider the ions
that were added to neutralize the solute.  It assumes those are bound to the
solute and do not contribute to the bulk ionic strength.

By default, Na\ :sup:`+` and Cl\ :sup:`-` ions are used, but you can specify
different ones using the :code:`positiveIon` and :code:`negativeIon`
options.  For example, this creates a potassium chloride solution:
::

    modeller.addSolvent(forcefield, ionicStrength=0.1*molar, positiveIon='K+')

Allowed values for :code:`positiveIon` are ``'Cs+'``, ``'K+'``, ``'Li+'``, ``'Na+'``, and
``'Rb+'``.  Allowed values for :code:`negativeIon` are ``'Cl-'``, ``'Br-'``, ``'F-'``, and
``'I-'``.  Be sure to include the single quotes around the value.  Also be aware
some force fields do not include parameters for all of these ion types, so you
need to use types that are supported by your chosen force field.

Adding a Membrane
*****************

If you want to simulate a membrane protein, you may need to create a membrane as
well.  You can do this by calling :meth:`addMembrane`.  Call it *instead* of
:meth:`addSolvent`, not in addition to it.  This one method adds the membrane,
solvent, and ions all at once, making sure the lipid head groups are properly
solvated.  For example, this creates a POPC membrane, ensuring at least 1 nm of
padding on all sides:
::

    modeller.addMembrane(forcefield, lipidType='POPC', minimumPadding=1*nanometer)

The membrane is added in the XY plane, and the existing protein is assumed to already be oriented
and positioned correctly.  When possible, it is recommended to start with a model
from the `Orientations of Proteins in Membranes`_ (OPM) database.  Otherwise, it
is up to you to select the protein position yourself.

Because this method also adds solvent, it takes many of the same arguments as
:meth:`addSolvent`.  See the API documentation for details.

.. _`Orientations of Proteins in Membranes`: http://opm.phar.umich.edu

.. _adding-or-removing-extra-particles:

Adding or Removing Extra Particles
**********************************

“Extra particles” are particles that do not represent ordinary atoms.  This
includes the virtual interaction sites used in many water models, Drude
particles, etc.  If you are using a force field that involves extra particles,
you must add them to the :class:`Topology`.  To do this, call:
::

    modeller.addExtraParticles(forcefield)

This looks at the force field to determine what extra particles are needed, then
modifies each residue to include them.  This function can remove extra particles
as well as adding them.

Removing Water
**************

Call deleteWater to remove all water molecules from the system:
::

    modeller.deleteWater()

This is useful, for example, if you want to simulate it with implicit solvent.
Be aware, though, that this only removes water molecules, not ions or other
small molecules that might be considered “solvent”.

.. _saving-the-results:

Saving The Results
******************

Once you have finished editing your model, you can immediately use the resulting
:class:`Topology` object and atom positions as the input to a :class:`Simulation`.  If you plan to
simulate it many times, though, it is usually better to save the result to a new
PDB file, then use that as the input for the simulations.  This avoids the cost
of repeating the modelling operations at the start of every simulation, and also
ensures that all your simulations are really starting from exactly the same
structure.

The following example loads a PDB file, adds missing hydrogens, builds a solvent
box around it, performs an energy minimization, and saves the result to a new
PDB file.

.. samepage::
    ::

        from openmm.app import *
        from openmm import *
        from openmm.unit import *

        print('Loading...')
        pdb = PDBFile('input.pdb')
        forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
        modeller = Modeller(pdb.topology, pdb.positions)
        print('Adding hydrogens...')
        modeller.addHydrogens(forcefield)
        print('Adding solvent...')
        modeller.addSolvent(forcefield, model='tip3p', padding=1*nanometer)
        print('Minimizing...')
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
        integrator = VerletIntegrator(0.001*picoseconds)
        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)
        simulation.minimizeEnergy(maxIterations=100)
        print('Saving...')
        positions = simulation.context.getState(positions=True).getPositions()
        PDBFile.writeFile(simulation.topology, positions, open('output.pdb', 'w'))
        print('Done')

    .. caption::

        :autonumber:`Example,Modeller complete`


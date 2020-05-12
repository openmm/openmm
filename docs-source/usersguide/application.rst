.. include:: header.rst

.. default-domain:: py

.. _the-openmm-application-layer-introduction:

The OpenMM Application Layer: Getting Started
#############################################

Introduction
************

The first thing to understand about the OpenMM “application layer” is that it is
not exactly an application in the traditional sense: there is no program called
“OpenMM” that you run.  Rather, it is a collection of libraries written in the
Python programming language.  Those libraries can easily be chained together to
create Python programs that run simulations.  But don’t worry!  You don’t need
to know anything about Python programming (or programming at all) to use it.
Nearly all molecular simulation applications ask you to write some sort of
“script” that specifies the details of the simulation to run.  With OpenMM, that
script happens to be written in Python.  But it is no harder to write than those
for most other applications, and this guide will teach you everything you need
to know.  There is even a graphical interface that can write the script for you
based on a simple set of options (see Section :ref:`the-script-builder-application`),
so you never need to type a single line of code!

On the other hand, if you don’t mind doing a little programming, this approach
gives you enormous power and flexibility.  Your script has complete access to
the entire OpenMM application programming interface (API), as well as the full
power of the Python language and libraries.  You have complete control over
every detail of the simulation, from defining the molecular system to analyzing
the results.


.. _installing-openmm:

Installing OpenMM
*****************

OpenMM is installed using the Conda package manager (http://conda.pydata.org).
Conda is included as part of the Anaconda Python distribution, which you can
download from http://docs.continuum.io/anaconda/install.  This is a Python
distribution specifically designed for scientific applications, with many of the
most popular mathematical and scientific packages preinstalled.  Alternatively
you can use Miniconda (available from http://conda.pydata.org/miniconda.html),
which includes only Python itself, plus the Conda package manager.  That offers
a much smaller initial download, with the ability to then install only the
packages you want.

(A third option is to compile OpenMM from source.  This provides more flexibility,
but it is much more work, and there is rarely a need for anyone but advanced users
to compile from source.  Detailed instruction are in Chapter :ref:`compiling-openmm-from-source-code`.)

\1. Begin by installing the most recent 64 bit, Python 3.x version of either
Anaconda or Miniconda.

\2. (Optional) If you want to run OpenMM on a GPU, install CUDA and/or OpenCL.

  * If you have an Nvidia GPU, download CUDA from
    https://developer.nvidia.com/cuda-downloads.  Be sure to install both the
    drivers and toolkit.  OpenCL is included with the CUDA drivers.
  * If you have an AMD GPU and are using Linux or Windows, download the latest
    version of the Catalyst driver from http://support.amd.com.  On OS X, OpenCL
    is included with the operating system and is supported on OS X 10.10.3 or
    later.

3. Open a command line terminal and type the following command
::

    conda install -c omnia -c conda-forge openmm

This installs a version of OpenMM that is compiled to work with CUDA 10.1.
Alternatively you can request a version that is compiled for a specific CUDA
version with the command
::

    conda install -c omnia/label/cuda92 -c conda-forge openmm

where :code:`cuda92` should be replaced with the particular CUDA version
installed on your computer.  Supported values are :code:`cuda75`, :code:`cuda80`,
:code:`cuda90`, :code:`cuda91`, :code:`cuda92`, :code:`cuda100`, and :code:`cuda101`.  Because
different CUDA releases are not binary compatible with each other, OpenMM can
only work with the particular CUDA version it was compiled with.

4. Verify your installation by typing the following command:
::

    python -m simtk.testInstallation

This command confirms that OpenMM is installed, checks whether GPU acceleration
is available (via the OpenCL and/or CUDA platforms), and verifies that all
platforms produce consistent results.


.. _running-simulations:

Running Simulations
###################

.. _a-first-example:

A First Example
***************

Let’s begin with our first example of an OpenMM script. It loads a PDB file
called :file:`input.pdb` that defines a biomolecular system, parameterizes it using the Amber14 force field and TIP3P-FB water
model, energy minimizes it, simulates it for 10,000 steps with a Langevin
integrator, and saves a snapshot frame to a PDB file called :file:`output.pdb` every 1000 time
steps.

.. samepage::
    ::

        from simtk.openmm.app import *
        from simtk.openmm import *
        from simtk.unit import *
        from sys import stdout

        pdb = PDBFile('input.pdb')
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                nonbondedCutoff=1*nanometer, constraints=HBonds)
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter('output.pdb', 1000))
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                potentialEnergy=True, temperature=True))
        simulation.step(10000)

    .. caption::

        :autonumber:`Example,PDB example`

You can find this script in the :file:`examples` folder of your OpenMM installation.
It is called :file:`simulatePdb.py`.  To execute it from a command line, go to your
terminal/console/command prompt window (see Section :ref:`installing-openmm`
on setting up the window to use OpenMM).  Navigate to the :file:`examples` folder by typing
::

    cd <examples_directory>

where the typical directory is :file:`/usr/local/openmm/examples` on Linux
and Mac machines and  :file:`C:\\Program Files\\OpenMM\\examples` on Windows
machines.

Then type
::

    python simulatePdb.py

You can name your own scripts whatever you want.  Let’s go through the script line
by line and see how it works.
::

    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
    from sys import stdout

These lines are just telling the Python interpreter about some libraries we will
be using.  Don’t worry about exactly what they mean.  Just include them at the
start of your scripts.
::

    pdb = PDBFile('input.pdb')

This line loads the PDB file from disk.  (The :file:`input.pdb` file in the :file:`examples`
directory contains the villin headpiece in explicit solvent.)  More precisely,
it creates a :class:`PDBFile` object, passes the file name :file:`input.pdb` to it as an
argument, and assigns the object to a variable called :code:`pdb`\ .  The
:class:`PDBFile` object contains the information that was read from the file: the
molecular topology and atom positions.  Your file need not be called
:file:`input.pdb`.  Feel free to change this line to specify any file you want,
though it must contain all of the atoms needed by the force field.
(More information on how to add missing atoms and residues using OpenMM tools can be found in Chapter :ref:`model-building-and-editing`.)
Make sure you include the single quotes around the file name.  OpenMM also can load
files in the newer PDBx/mmCIF format: just change :class:`PDBFile` to :class:`PDBxFile`.
::

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

This line specifies the force field to use for the simulation.  Force fields are
defined by XML files.  OpenMM includes XML files defining lots of standard force fields (see Section :ref:`force-fields`).
If you find you need to extend the repertoire of force fields available,
you can find more information on how to create these XML files in Chapter :ref:`creating-force-fields`.
In this case we load two of those files: :file:`amber14-all.xml`, which contains the
Amber14 force field, and :file:`amber14/tip3pfb.xml`, which contains the TIP3P-FB water model.  The
:class:`ForceField` object is assigned to a variable called :code:`forcefield`\ .
::

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer, constraints=HBonds)

This line combines the force field with the molecular topology loaded from the
PDB file to create a complete mathematical description of the system we want to
simulate.  (More precisely, we invoke the :class:`ForceField` object’s :meth:`.createSystem`
function.  It creates a :class:`System` object, which we assign to the variable
:code:`system`\ .)  It specifies some additional options about how to do that:
use particle mesh Ewald for the long range electrostatic interactions
(:code:`nonbondedMethod=PME`\ ), use a 1 nm cutoff for the direct space
interactions (\ :code:`nonbondedCutoff=1*nanometer`\ ), and constrain the length
of all bonds that involve a hydrogen atom (\ :code:`constraints=HBonds`\ ).
Note the way we specified the cutoff distance 1 nm using :code:`1*nanometer`:
This is an example of the powerful units tracking and automatic conversion facility
built into the OpenMM Python API that makes specifying unit-bearing quantities
convenient and less error-prone.  We could have equivalently specified
:code:`10*angstrom` instead of :code:`1*nanometer` and achieved the same result.
The units system will be described in more detail later, in Section :ref:`units-and-dimensional-analysis`.
::

    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

This line creates the integrator to use for advancing the equations of motion.
It specifies a :class:`LangevinIntegrator`, which performs Langevin dynamics,
and assigns it to a variable called :code:`integrator`\ .  It also specifies
the values of three parameters that are specific to Langevin dynamics: the
simulation temperature (300 K), the friction coefficient (1 ps\ :sup:`-1`\ ), and
the step size (0.002 ps).
::

    simulation = Simulation(pdb.topology, system, integrator)

This line combines the molecular topology, system, and integrator to begin a new
simulation.  It creates a :class:`Simulation` object and assigns it to a variable called
\ :code:`simulation`\ .  A :class:`Simulation` object manages all the processes
involved in running a simulation, such as advancing time and writing output.
::

    simulation.context.setPositions(pdb.positions)

This line specifies the initial atom positions for the simulation: in this case,
the positions that were loaded from the PDB file.
::

    simulation.minimizeEnergy()

This line tells OpenMM to perform a local energy minimization.  It is usually a
good idea to do this at the start of a simulation, since the coordinates in the
PDB file might produce very large forces.
::

    simulation.reporters.append(PDBReporter('output.pdb', 1000))

This line creates a “reporter” to generate output during the simulation, and
adds it to the :class:`Simulation` object’s list of reporters.  A :class:`PDBReporter` writes
structures to a PDB file.  We specify that the output file should be called
:file:`output.pdb`, and that a structure should be written every 1000 time steps.
::

    simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
            potentialEnergy=True, temperature=True))

It can be useful to get regular status reports as a simulation runs so you can
monitor its progress.  This line adds another reporter to print out some basic
information every 1000 time steps: the current step index, the potential energy
of the system, and the temperature.  We specify :code:`stdout` (not in
quotes) as the output file, which means to write the results to the console.  We
also could have given a file name (in quotes), just as we did for the
:class:`PDBReporter`, to write the information to a file.
::

    simulation.step(10000)

Finally, we run the simulation, integrating the equations of motion for 10,000
time steps.  Once it is finished, you can load the PDB file into any program you
want for analysis and visualization (VMD_, PyMol_, AmberTools_, etc.).

.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _PyMol: http://www.pymol.org
.. _AmberTools: http://ambermd.org

.. _using_amber_files:

Using AMBER Files
*****************

OpenMM can build a system in several different ways.  One option, as shown
above, is to start with a PDB file and then select a force field with which to
model it.  Alternatively, you can use AmberTools_ to model your system.  In that
case, you provide a :class:`prmtop` file and an :class:`inpcrd` file.  OpenMM loads the files and
creates a :class:`System` from them.  This is illustrated in the following script.  It can be
found in OpenMM’s :file:`examples` folder with the name :file:`simulateAmber.py`.

.. samepage::
    ::

        from simtk.openmm.app import *
        from simtk.openmm import *
        from simtk.unit import *
        from sys import stdout

        prmtop = AmberPrmtopFile('input.prmtop')
        inpcrd = AmberInpcrdFile('input.inpcrd')
        system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                constraints=HBonds)
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(prmtop.topology, system, integrator)
        simulation.context.setPositions(inpcrd.positions)
        if inpcrd.boxVectors is not None:
            simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter('output.pdb', 1000))
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                potentialEnergy=True, temperature=True))
        simulation.step(10000)

    .. caption::

        :autonumber:`Example,AMBER example`

This script is very similar to the previous one.  There are just a few
significant differences:
::

    prmtop = AmberPrmtopFile('input.prmtop')
    inpcrd = AmberInpcrdFile('input.inpcrd')

In these lines, we load the prmtop file and inpcrd file.  More precisely, we
create :class:`AmberPrmtopFile` and :class:`AmberInpcrdFile` objects and assign them to the
variables :code:`prmtop` and :code:`inpcrd`\ , respectively.  As before,
you can change these lines to specify any files you want.  Be sure to include
the single quotes around the file names.

.. note::

    The :class:`AmberPrmtopFile` reader provided by OpenMM only supports "new-style"
    :file:`prmtop` files introduced in AMBER 7. The AMBER distribution still contains a number of
    example files that are in the "old-style" :file:`prmtop` format. These "old-style" files will
    not run in OpenMM.

Next, the :class:`System` object is created in a different way:
::

    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
            constraints=HBonds)

In the previous section, we loaded the topology
from a PDB file and then had the force field create a system based on it.  In
this case, we don’t need a force field; the :file:`prmtop` file already contains the
force field parameters, so it can create the system
directly.
::

    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)

Notice that we now get the topology from the :file:`prmtop` file and the atom positions
from the :file:`inpcrd` file.  In the previous section, both of these came from a PDB
file, but AMBER puts the topology and positions in separate files.  We also add the
following lines:
::

    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

For periodic systems, the :file:`prmtop` file specifies the periodic box vectors, just
as a PDB file does.  When we call :meth:`createSystem`, it sets those as the default
periodic box vectors, to be used automatically for all simulations.  However, the
:file:`inpcrd` may *also* specify periodic box vectors,
and if so we want to use those ones instead.  For example, if the system has been
equilibrated with a barostat, the box vectors may have changed during equilibration.
We therefore check to see if the :file:`inpcrd` file contained box vectors.  If so,
we call :meth:`setPeriodicBoxVectors` to tell it to use those ones, overriding the
default ones provided by the :class:`System`.

.. _using_gromacs_files:

Using Gromacs Files
*******************

A third option for creating your system is to use the Gromacs setup tools.  They
produce a :file:`gro` file containing the coordinates and a :file:`top` file containing the
topology.  OpenMM can load these exactly as it did the AMBER files.  This is
shown in the following script.  It can be found in OpenMM’s :file:`examples` folder
with the name :file:`simulateGromacs.py`.

.. samepage::
    ::

        from simtk.openmm.app import *
        from simtk.openmm import *
        from simtk.unit import *
        from sys import stdout

        gro = GromacsGroFile('input.gro')
        top = GromacsTopFile('input.top', periodicBoxVectors=gro.getPeriodicBoxVectors(),
                includeDir='/usr/local/gromacs/share/gromacs/top')
        system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                constraints=HBonds)
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(top.topology, system, integrator)
        simulation.context.setPositions(gro.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter('output.pdb', 1000))
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                potentialEnergy=True, temperature=True))
        simulation.step(10000)

    .. caption::

        :autonumber:`Example,Gromacs example`

This script is nearly identical to the previous one, just replacing
:class:`AmberInpcrdFile` and :class:`AmberPrmtopFile` with :class:`GromacsGroFile` and :class:`GromacsTopFile`.
Note that when we create the :class:`GromacsTopFile`, we specify values for two extra
options.  First, we specify
:code:`periodicBoxVectors=gro.getPeriodicBoxVectors()`\ .  Unlike OpenMM and
AMBER, which can store periodic unit cell information with the topology, Gromacs
only stores it with the coordinates.  To let :class:`GromacsTopFile` create a :class:`Topology`
object, we therefore need to tell it the periodic box vectors that were loaded
from the :file:`gro` file.  You only need to do this if you are simulating a periodic
system.  For implicit solvent simulations, it usually can be omitted.

Second, we specify :code:`includeDir='/usr/local/gromacs/share/gromacs/top'`\ .  Unlike AMBER,
which stores all the force field parameters directly in a :file:`prmtop` file, Gromacs just stores
references to force field definition files that are installed with the Gromacs
application.  OpenMM needs to know where to find these files, so the
:code:`includeDir` parameter specifies the directory containing them.  If you
omit this parameter, OpenMM will assume the default location :file:`/usr/local/gromacs/share/gromacs/top`,
which is often where they are installed on
Unix-like operating systems.  So in :autonumref:`Example,Gromacs example` we actually could have omitted
this parameter, but if the Gromacs files were installed in any other location,
we would need to include it.

.. _using-charmm-files:

Using CHARMM Files
******************

Yet another option is to load files created by the CHARMM setup tools, or other compatible
tools such as VMD.  Those include a :file:`psf` file containing topology information, and an
ordinary PDB file for the atomic coordinates.  (Coordinates can also be loaded from CHARMM
coordinate or restart files using the :class:`CharmmCrdFile` and :class:`CharmmRstFile` classes).  In addition,
you must provide a set of files containing the force
field definition to use.  This can involve several different files with varying formats and
filename extensions such as :file:`par`, :file:`prm`, :file:`top`, :file:`rtf`, :file:`inp`,
and :file:`str`.  To do this, load all the definition files into a :class:`CharmmParameterSet`
object, then include that object as the first parameter when you call :meth:`createSystem`
on the :class:`CharmmPsfFile`.

.. samepage::
    ::

        from simtk.openmm.app import *
        from simtk.openmm import *
        from simtk.unit import *
        from sys import stdout, exit, stderr

        psf = CharmmPsfFile('input.psf')
        pdb = PDBFile('input.pdb')
        params = CharmmParameterSet('charmm22.rtf', 'charmm22.prm')
        system = psf.createSystem(params, nonbondedMethod=NoCutoff,
                nonbondedCutoff=1*nanometer, constraints=HBonds)
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(psf.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter('output.pdb', 1000))
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                potentialEnergy=True, temperature=True))
        simulation.step(10000)

    .. caption::

        :autonumber:`Example,CHARMM example`

Note that both the CHARMM and XPLOR versions of the :file:`psf` file format are supported.

.. _the-script-builder-application:

The OpenMM-Setup Application
****************************

One way to create your own scripts is to start with one of the examples given
above and customize it to suit your needs, but there's an even easier option.
OpenMM-Setup is a graphical application that walks you through the whole process
of loading your input files and setting options.  It then generates a complete
script, and can even run it for you.

.. figure:: ../images/OpenMMSetup.png
   :align: center
   :width: 100%

   :autonumber:`Figure,openmm setup`:  The OpenMM-Setup application

To install OpenMM-Setup, open a command line terminal and type the following command
::

    conda install -c omnia openmm-setup

You can then launch it by typing the command
::

    openmm-setup

It will automatically open a window in your web browser displaying the user interface.

OpenMM-Setup is far more than just a script generator.  It can fix problems in
your input files, add missing atoms, build membranes and water boxes, and much
more.  It is a very easy way to quickly do all necessary preparation and setup.
We highly recommend it to all users of OpenMM, from novices to experts.

.. _simulation-parameters:

Simulation Parameters
*********************

Now let’s consider lots of ways you might want to customize your script.

Platforms
=========

When creating a :class:`Simulation`, you can optionally tell it what :class:`Platform` to use.
OpenMM includes four platforms: :class:`Reference`, :class:`CPU`, :class:`CUDA`, and :class:`OpenCL`.  For a
description of the differences between them, see Section :ref:`platforms`.  There are three ways in which
the :class:`Platform` can be chosen:

1. By default, OpenMM will try to select the fastest available :class:`Platform`.  Usually its choice will
be reasonable, but sometimes you may want to change it.

2. Alternatively, you can set the :envvar:`OPENMM_DEFAULT_PLATFORM` environment variable to the name
of the :class:`Platform` to use.  This overrides the default logic.

3. Finally, you can explicitly specify a :class:`Platform` object in your script when you create the
:class:`Simulation`.  The following lines specify to use the :class:`CUDA` platform:
::

    platform = Platform.getPlatformByName('CUDA')
    simulation = Simulation(prmtop.topology, system, integrator, platform)

The platform name should be one of :code:`OpenCL`, :code:`CUDA`, :code:`CPU`, or
:code:`Reference`.

You also can specify platform-specific properties that customize how
calculations should be done.  See Chapter :ref:`platform-specific-properties` for details of the
properties that each Platform supports.  For example, the following lines specify to parallelize
work across two different GPUs (CUDA devices 0 and 1), doing all computations in
double precision:
::

    platform = Platform.getPlatformByName('CUDA')
    properties = {'DeviceIndex': '0,1', 'Precision': 'double'}
    simulation = Simulation(prmtop.topology, system, integrator, platform, properties)

.. _force-fields:

Force Fields
============

When you create a force field, you specify one or more XML files from which to
load the force field definition.  Most often, there will be one file to define
the main force field, and possibly a second file to define the water model
(either implicit or explicit).  For example:
::

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

In some cases, one XML file may load several others.  For example, :file:`amber14-all.xml`
is really just a shortcut for loading several different files that together make up
the AMBER14 force field.  If you need finer grained control over which parameters
are loaded, you can instead specify the component files individually.

Be aware that some force fields and water models include "extra particles", such
as lone pairs or Drude particles.  Examples include the CHARMM polarizable force
field and all of the 4 and 5 site water models.  To use these force fields, you
must first add the extra particles to the :class:`Topology`.  See section
:ref:`adding-or-removing-extra-particles` for details.

The force fields described below are the ones that are bundled with OpenMM.
Additional force fields are available online at https://github.com/choderalab/openmm-forcefields.

Amber14
-------

The Amber14\ :cite:`Maier2015` force field is made up of various files that define
parameters for proteins, DNA, RNA, lipids, water, and ions.

.. tabularcolumns:: |l|L|

===================================  ============================================
File                                 Parameters
===================================  ============================================
:file:`amber14/protein.ff14SB.xml`   Protein (recommended)
:file:`amber14/protein.ff15ipq.xml`  Protein (alternative)
:file:`amber14/DNA.OL15.xml`         DNA (recommended)
:file:`amber14/DNA.bsc1.xml`         DNA (alternative)
:file:`amber14/RNA.OL3.xml`          RNA
:file:`amber14/lipid17.xml`          Lipid
:file:`amber14/tip3p.xml`            TIP3P water model\ :cite:`Jorgensen1983` and ions
:file:`amber14/tip3pfb.xml`          TIP3P-FB water model\ :cite:`Wang2014` and ions
:file:`amber14/tip4pew.xml`          TIP4P-Ew water model\ :cite:`Horn2004` and ions
:file:`amber14/tip4pfb.xml`          TIP4P-FB water model\ :cite:`Wang2014` and ions
:file:`amber14/spce.xml`             SPC/E water model\ :cite:`Berendsen1987` and ions
===================================  ============================================

As a convenience, the file :file:`amber14-all.xml` can be used as a shortcut to include
:file:`amber14/protein.ff14SB.xml`, :file:`amber14/DNA.OL15.xml`, :file:`amber14/RNA.OL3.xml`,
and :file:`amber14/lipid17.xml`.  In most cases, you can simply include that file,
plus one of the water models, such as :file:`amber14/tip3pfb.xml` for the TIP3P-FB
water model and ions\ :cite:`Wang2014`:
::

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

.. tip:: The solvent model XML files included under the :file:`amber14/` directory
         include both water *and* ions compatible with that water model, so if you
         mistakenly specify :file:`tip3p.xml` instead of :file:`amber14/tip3p.xml`,
         you run the risk of having :class:`ForceField` throw an exception since
         :file:`tip3p.xml` will be missing parameters for ions in your system.

The converted parameter sets come from the `AmberTools 17 release <http://ambermd.org/AmberTools17-get.html>`_
and were converted using the `openmm-forcefields <https://github.com/choderalab/openmm-forcefields>`_ package and `ParmEd <https://github.com/parmed/parmed>`_.

CHARMM36
--------

The CHARMM36\ :cite:`Best2012` force field provides parameters for proteins, DNA,
RNA, lipids, carbohydrates, water, ions, and various small molecules (see `here <http://mackerell.umaryland.edu/charmm_ff.shtml#refs>`_
for full references).

.. tabularcolumns:: |l|L|

=================================  ============================================
File                               Parameters
=================================  ============================================
:file:`charmm36.xml`               Protein, DNA, RNA, lipids, carbohydrates, and small molecules
:file:`charmm36/water.xml`         Default CHARMM water model (a modified version of TIP3P\ :cite:`Jorgensen1983`) and ions
:file:`charmm36/spce.xml`          SPC/E water model\ :cite:`Berendsen1987` and ions
:file:`charmm36/tip3p-pme-b.xml`   TIP3P-PME-B water model\ :cite:`Price2004` and ions
:file:`charmm36/tip3p-pme-f.xml`   TIP3P-PME-F water model\ :cite:`Price2004` and ions
:file:`charmm36/tip4pew.xml`       TIP4P-Ew water model\ :cite:`Horn2004` and ions
:file:`charmm36/tip4p2005.xml`     TIP4P-2005 water model\ :cite:`Abascal2005` and ions
:file:`charmm36/tip5p.xml`         TIP5P water model\ :cite:`Mahoney2000` and ions
:file:`charmm36/tip5pew.xml`       TIP5P-Ew water model\ :cite:`Rick2004` and ions
=================================  ============================================

The file :file:`charmm36.xml` bundles everything but the water and ions into a single
file.  In most cases, you can simply include that file, plus one of the water models,
such as :file:`charmm36/water.xml`, which specifies the default CHARMM water model
(a modified version of TIP3P\ :cite:`Jorgensen1983`) and ions:
::

    forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')

.. warning:: Drude polarizable sites and lone pairs are not yet supported
             by `ParmEd <https://github.com/parmed/parmed>`_ and the CHARMM36 forcefields
             that depend on these features are not included in this port.
             To use the CHARMM 2013 polarizable force field\ :cite:`Lopes2013`,
             include the single file :file:`charmm_polar_2013.xml`.

.. tip:: The solvent model XML files included under the :file:`charmm36/` directory
         include both water *and* ions compatible with that water model, so if you
         mistakenly specify :file:`tip3p.xml` instead of :file:`charmm36/water.xml`,
         you run the risk of having :class:`ForceField` raise an exception due to
         missing parameters for ions in your system.

.. tip:: CHARMM makes extensive use of patches, which are automatically combined with
         residue templates to create an expanded library of patched residue templates
         by :class:`ForceField`. That means that patched residues, such as ``ACE`` and
         ``NME`` patched termini, must occur as a single residue in order for :class:`ForceField`
         to correctly match the residue template and apply parameters. Since these
         patched residues are not standard PDB residues, :class:`Modeller` does not know
         how to add hydrogens to these nonstandard residues, and your input topologies
         must already contain appropriate hydrogens. This can often cause problems when
         trying to read in PDB files from sources such as `CHARMM-GUI <http://charmm-gui.org/>`_
         that do not generate PDB files that comply with the `PDB standard <http://www.wwpdb.org/documentation/file-format>`_.
         If you're using files from `CHARMM-GUI <http://charmm-gui.org/>`_, it's easiest to load
         the PSF file directly, as discussed in Section :ref:`using-charmm-files`.

.. tip:: Trying to read in PDB files from sources such as `CHARMM-GUI <http://charmm-gui.org/>`_
         that do not generate PDB files that comply with the `PDB standard <http://www.wwpdb.org/documentation/file-format>`_
         and omit ``CONECT`` records specifying bonds between residues (such as cysteines)
         or include ``CONECT`` records specifying non-chemical ``H-H`` bonds in waters
         can cause issues with the detection and parameter assignment for disulfide bonds.
         Make sure the files you read in comply with the appropriate standards regarding
         additional bonds and nonstandard residue definitions. If you're using files from
         `CHARMM-GUI <http://charmm-gui.org/>`_, it's easiest to load
         the PSF file directly, as discussed in Section :ref:`using-charmm-files`.

The converted parameter sets come from the `CHARMM36 July 2017 update <http://mackerell.umaryland.edu/charmm_ff.shtml>`_
and were converted using the `openmm-forcefields <https://github.com/choderalab/openmm-forcefields>`_ package and `parmed <https://github.com/parmed/parmed>`_.

AMOEBA
------

The AMOEBA polarizable force field provides parameters for proteins, water, and ions.

.. tabularcolumns:: |l|L|

=============================  ================================================================================
File                           Parameters
=============================  ================================================================================
:file:`amoeba2013.xml`         AMOEBA 2013\ :cite:`Shi2013`
:file:`amoeba2013_gk.xml`      Generalized Kirkwood solvation model\ :cite:`Schnieders2007` for use with AMOEBA 2013 force field
:file:`amoeba2009.xml`         AMOEBA 2009\ :cite:`Ren2002`.  This force field is deprecated.  It is
                               recommended to use AMOEBA 2013 instead.
:file:`amoeba2009_gk.xml`      Generalized Kirkwood solvation model for use with AMOEBA 2009 force field
=============================  ================================================================================

For explicit solvent simulations, just include the single file :file:`amoeba2013.xml`.
AMOEBA also supports implicit solvent using a Generalized Kirkwood model.  To enable
it, also include :file:`amoeba2013_gk.xml`.

The older AMOEBA 2009 force field is provided only for backward compatibility, and is not
recommended for most simulations.

CHARMM Polarizable Force Field
------------------------------

To use the CHARMM 2013 polarizable force field\ :cite:`Lopes2013`, include the
single file :file:`charmm_polar_2013.xml`.  It includes parameters for proteins,
water, and ions.  When using this force field, remember to add extra particles to
the :class:`Topology` as described in section :ref:`adding-or-removing-extra-particles`.

Older Amber Force Fields
------------------------

OpenMM includes several older Amber force fields as well.  For most simulations
Amber14 is preferred over any of these, but they are still useful for reproducing
older results.

.. tabularcolumns:: |l|L|

=============================  ================================================================================
File                           Force Field
=============================  ================================================================================
:code:`amber96.xml`            Amber96\ :cite:`Kollman1997`
:code:`amber99sb.xml`          Amber99\ :cite:`Wang2000` with modified backbone torsions\ :cite:`Hornak2006`
:code:`amber99sbildn.xml`      Amber99SB plus improved side chain torsions\ :cite:`Lindorff-Larsen2010`
:code:`amber99sbnmr.xml`       Amber99SB with modifications to fit NMR data\ :cite:`Li2010`
:code:`amber03.xml`            Amber03\ :cite:`Duan2003`
:code:`amber10.xml`            Amber10 (documented in the AmberTools_ manual as `ff10`)
=============================  ================================================================================

Several of these force fields support implicit solvent.  To enable it, also
include the corresponding OBC file.

.. tabularcolumns:: |l|L|

=========================  =================================================================================================
File                       Implicit Solvation Model
=========================  =================================================================================================
:code:`amber96_obc.xml`    GBSA-OBC solvation model\ :cite:`Onufriev2004` for use with Amber96 force field
:code:`amber99_obc.xml`    GBSA-OBC solvation model for use with Amber99 force field and its variants
:code:`amber03_obc.xml`    GBSA-OBC solvation model for use with Amber03 force field
:code:`amber10_obc.xml`    GBSA-OBC solvation model for use with Amber10 force field
=========================  =================================================================================================

Note that the GBSA-OBC parameters in these files are those used in TINKER.\ :cite:`Tinker`
They are designed for use with Amber force fields, but they are different from
the parameters found in the AMBER application.

Water Models
------------

The following files define popular water models.  They can be used with force fields
that do not provide their own water models.  When using Amber14 or CHARMM36, use
the water files included with those force fields instead, since they also include
ion parameters.

.. tabularcolumns:: |l|L|

===================  ============================================
File                 Water Model
===================  ============================================
:code:`tip3p.xml`    TIP3P water model\ :cite:`Jorgensen1983`
:code:`tip3pfb.xml`  TIP3P-FB water model\ :cite:`Wang2014`
:code:`tip4pew.xml`  TIP4P-Ew water model\ :cite:`Horn2004`
:code:`tip4pfb.xml`  TIP4P-FB water model\ :cite:`Wang2014`
:code:`tip5p.xml`    TIP5P water model\ :cite:`Mahoney2000`
:code:`spce.xml`     SPC/E water model\ :cite:`Berendsen1987`
:code:`swm4ndp.xml`  SWM4-NDP water model\ :cite:`Lamoureux2006`
===================  ============================================


AMBER Implicit Solvent
======================


When creating a system from a prmtop file you do not specify force field files,
so you need a different way to tell it to use implicit solvent.  This is done
with the :code:`implicitSolvent` parameter:
::

    system = prmtop.createSystem(implicitSolvent=OBC2)

OpenMM supports all of the Generalized Born models used by AMBER.  Here are the
allowed values for :code:`implicitSolvent`\ :

.. tabularcolumns:: |l|L|

=============  ==================================================================================================================================
Value          Meaning
=============  ==================================================================================================================================
:code:`None`   No implicit solvent is used.
:code:`HCT`    Hawkins-Cramer-Truhlar GBSA model\ :cite:`Hawkins1995` (corresponds to igb=1 in AMBER)
:code:`OBC1`   Onufriev-Bashford-Case GBSA model\ :cite:`Onufriev2004` using the GB\ :sup:`OBC`\ I parameters (corresponds to igb=2 in AMBER).
:code:`OBC2`   Onufriev-Bashford-Case GBSA model\ :cite:`Onufriev2004` using the GB\ :sup:`OBC`\ II parameters (corresponds to igb=5 in AMBER).
               This is the same model used by the GBSA-OBC files described in Section :ref:`force-fields`.
:code:`GBn`    GBn solvation model\ :cite:`Mongan2007` (corresponds to igb=7 in AMBER).
:code:`GBn2`   GBn2 solvation model\ :cite:`Nguyen2013` (corresponds to igb=8 in AMBER).
=============  ==================================================================================================================================


You can further control the solvation model in a few ways.  First, you can
specify the dielectric constants to use for the solute and solvent:
::

    system = prmtop.createSystem(implicitSolvent=OBC2, soluteDielectric=1.0,
            solventDielectric=80.0)

If they are not specified, the solute and solvent dielectrics default to 1.0 and
78.5, respectively.  These values were chosen for consistency with AMBER, and
are slightly different from those used elsewhere in OpenMM: when building a
system from a force field, the solvent dielectric defaults to 78.3.

You also can model the effect of a non-zero salt concentration by specifying the
Debye-Huckel screening parameter\ :cite:`Srinivasan1999`:
::

    system = prmtop.createSystem(implicitSolvent=OBC2, implicitSolventKappa=1.0/nanometer)


Nonbonded Interactions
======================


When creating the system (either from a force field or a prmtop file), you can
specify options about how nonbonded interactions should be treated:
::

    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer)

The :code:`nonbondedMethod` parameter can have any of the following values:

.. tabularcolumns:: |l|L|

=========================  ===========================================================================================================================================================================================================================================
Value                      Meaning
=========================  ===========================================================================================================================================================================================================================================
:code:`NoCutoff`           No cutoff is applied.
:code:`CutoffNonPeriodic`  The reaction field method is used to eliminate all interactions beyond a cutoff distance.  Not valid for AMOEBA.
:code:`CutoffPeriodic`     The reaction field method is used to eliminate all interactions beyond a cutoff distance.  Periodic boundary conditions are applied, so each atom interacts only with the nearest periodic copy of every other atom.  Not valid for AMOEBA.
:code:`Ewald`              Periodic boundary conditions are applied.  Ewald summation is used to compute long range Coulomb interactions.  (This option is rarely used, since PME is much faster for all but the smallest systems.)  Not valid for AMOEBA.
:code:`PME`                Periodic boundary conditions are applied.  The Particle Mesh Ewald method is used to compute long range Coulomb interactions.
:code:`LJPME`              Periodic boundary conditions are applied.  The Particle Mesh Ewald method is used to compute long range interactions for both Coulomb and Lennard-Jones.
=========================  ===========================================================================================================================================================================================================================================


When using any method other than :code:`NoCutoff`\ , you should also specify a
cutoff distance.  Be sure to specify units, as shown in the examples above. For
example, :code:`nonbondedCutoff=1.5*nanometers` or
:code:`nonbondedCutoff=12*angstroms` are legal values.

When using :code:`Ewald`, :code:`PME`, or :code:`LJPME`\ , you can optionally specify an
error tolerance for the force computation.  For example:
::

    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
            ewaldErrorTolerance=0.00001)

The error tolerance is roughly equal to the fractional error in the forces due
to truncating the Ewald summation.  If you do not specify it, a default value of
0.0005 is used.

Another optional parameter when using a cutoff is :code:`switchDistance`.  This
causes Lennard-Jones interactions to smoothly go to zero over some finite range,
rather than being sharply truncated at the cutoff distance.  This can improve
energy conservation.  To use it, specify a distance at which the interactions
should start being reduced.  For example:
::

    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
            switchDistance=0.9*nanometer)


Nonbonded Forces for AMOEBA
---------------------------

For the AMOEBA force field, the valid values for the :code:`nonbondedMethod`
are :code:`NoCutoff` and :code:`PME`\ .  The other nonbonded methods,
:code:`CutoffNonPeriodic`\ , :code:`CutoffPeriodic`\ , and :code:`Ewald`
are unavailable for this force field.

For implicit solvent runs using AMOEBA, only the :code:`nonbondedMethod`
option :code:`NoCutoff` is available.

Lennard-Jones Interaction Cutoff Value
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition, for the AMOEBA force field a cutoff for the Lennard-Jones
interaction independent of the value used for the electrostatic interactions may
be specified using the keyword :code:`vdwCutoff`\ .
::

    system = forcefield.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
            ewaldErrorTolerance=0.00001, vdwCutoff=1.2*nanometer)

If :code:`vdwCutoff` is not specified, then the value of
:code:`nonbondedCutoff` is used for the Lennard-Jones interactions.

Specifying the Polarization Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using the AMOEBA force field, OpenMM allows the induced dipoles to be
calculated in any of three different ways.  The slowest but potentially most
accurate method is to iterate the calculation until the dipoles converge to a
specified tolerance.  To select this, specify :code:`polarization='mutual'`.
Use the :code:`mutualInducedTargetEpsilon` option to select the tolerance; for
most situations, a value of 0.00001 works well.  Alternatively you can specify
:code:`polarization='extrapolated'`.  This uses an analytic approximation
:cite:`Simmonett2015` to estimate what the fully converged dipoles will be without
actually continuing the calculation to convergence.  In many cases this can be
significantly faster with only a small loss in accuracy.  Finally, you can
specify :code:`polarization='direct'` to use the direct polarization
approximation, in which induced dipoles depend only on the fixed multipoles, not
on other induced dipoles.  This is even faster, but it produces very different
forces from mutual polarization, so it should only be used with force fields
that have been specifically parameterized for use with this approximation.

Here are examples of using each method:
::

    system = forcefield.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        vdwCutoff=1.2*nanometer, polarization='mutual', mutualInducedTargetEpsilon=0.00001)

    system = forcefield.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        vdwCutoff=1.2*nanometer, polarization='extrapolated')

    system = forcefield.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        vdwCutoff=1.2*nanometer, polarization='direct')


Implicit Solvent and Solute Dielectrics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For implicit solvent simulations using the AMOEBA force field, the
:file:`amoeba2013_gk.xml` file should be included in the initialization of the force
field:
::

    forcefield = ForceField('amoeba2009.xml', 'amoeba2009_gk.xml')

Only the :code:`nonbondedMethod` option :code:`NoCutoff` is available
for implicit solvent runs using AMOEBA.  In addition, the solvent and solute
dielectric values can be specified for implicit solvent simulations:
::

    system=forcefield.createSystem(nonbondedMethod=NoCutoff, soluteDielectric=2.0,
            solventDielectric=80.0)

The default values are 1.0 for the solute dielectric and 78.3 for the solvent
dielectric.

Constraints
===========


When creating the system (either from a force field or an AMBER :file:`prmtop` file), you can
optionally tell OpenMM to constrain certain bond lengths and angles.  For
example,
::

    system = prmtop.createSystem(nonbondedMethod=NoCutoff, constraints=HBonds)

The :code:`constraints` parameter can have any of the following values:

.. tabularcolumns:: |l|L|

================  =============================================================================================================================================
Value             Meaning
================  =============================================================================================================================================
:code:`None`      No constraints are applied.  This is the default value.
:code:`HBonds`    The lengths of all bonds that involve a hydrogen atom are constrained.
:code:`AllBonds`  The lengths of all bonds are constrained.
:code:`HAngles`   The lengths of all bonds are constrained.  In addition, all angles of the form H-X-H or H-O-X (where X is an arbitrary atom) are constrained.
================  =============================================================================================================================================


The main reason to use constraints is that it allows one to use a larger
integration time step.  With no constraints, one is typically limited to a time
step of about 1 fs for typical biomolecular force fields like AMBER or CHARMM.  With :code:`HBonds` constraints, this can be increased
to about 2 fs.  With :code:`HAngles`\ , it can be further increased to 3.5 or
4 fs.

Regardless of the value of this parameter, OpenMM makes water molecules
completely rigid, constraining both their bond lengths and angles.  You can
disable this behavior with the :code:`rigidWater` parameter:
::

    system = prmtop.createSystem(nonbondedMethod=NoCutoff, constraints=None, rigidWater=False)

Be aware that flexible water may require you to further reduce the integration
step size, typically to about 0.5 fs.

.. note::

   The AMOEBA forcefield is intended to be used without constraints.

Heavy Hydrogens
===============


When creating the system (either from a force field or an AMBER :file:`prmtop` file), you can
optionally tell OpenMM to increase the mass of hydrogen atoms.  For example,
::

    system = prmtop.createSystem(hydrogenMass=4*amu)

This applies only to hydrogens that are bonded to heavy atoms, and any mass
added to the hydrogen is subtracted from the heavy atom.  This keeps their total
mass constant while slowing down the fast motions of hydrogens.  When combined
with constraints (typically :code:`constraints=AllBonds`\ ), this allows a
further increase in integration step size.

Integrators
===========


OpenMM offers a choice of several different integration methods.  You select
which one to use by creating an integrator object of the appropriate type.

Langevin Integrator
-------------------

In the examples of the previous sections, we used Langevin integration:
::

    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

The three parameter values in this line are the simulation temperature (300 K),
the friction coefficient (1 ps\ :sup:`-1`\ ), and the step size (0.002 ps).  You
are free to change these to whatever values you want.  Be sure to specify units
on all values.  For example, the step size could be written either as
:code:`0.002*picoseconds` or :code:`2*femtoseconds`\ .  They are exactly
equivalent.

Leapfrog Verlet Integrator
--------------------------

A leapfrog Verlet integrator can be used for running constant energy dynamics.
The command for this is:
::

    integrator = VerletIntegrator(0.002*picoseconds)

The only option is the step size.

Brownian Integrator
-------------------

Brownian (diffusive) dynamics can be used by specifying the following:
::

    integrator = BrownianIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

The parameters are the same as for Langevin dynamics: temperature (300 K),
friction coefficient (1 ps\ :sup:`-1`\ ), and step size (0.002 ps).

Variable Time Step Langevin Integrator
--------------------------------------

A variable time step Langevin integrator continuously adjusts its step size to
keep the integration error below a specified tolerance.  In some cases, this can
allow you to use a larger average step size than would be possible with a fixed
step size integrator.  It also is very useful in cases where you do not know in
advance what step size will be stable, such as when first equilibrating a
system.  You create this integrator with the following command:
::

    integrator = VariableLangevinIntegrator(300*kelvin, 1/picosecond, 0.001)

In place of a step size, you specify an integration error tolerance (0.001 in
this example).  It is best not to think of this value as having any absolute
meaning.  Just think of it as an adjustable parameter that affects the step size
and integration accuracy.  Smaller values will produce a smaller average step
size.  You should try different values to find the largest one that produces a
trajectory sufficiently accurate for your purposes.

Variable Time Step Leapfrog Verlet Integrator
---------------------------------------------

A variable time step leapfrog Verlet integrator works similarly to the variable
time step Langevin integrator in that it continuously adjusts its step size to
keep the integration error below a specified tolerance.  The command for this
integrator is:
::

    integrator = VariableVerletIntegrator(0.001)

The parameter is the integration error tolerance (0.001), whose meaning is the
same as for the Langevin integrator.

Multiple Time Step Integrator
-----------------------------

The :class:`MTSIntegrator` class implements the rRESPA multiple time step
algorithm\ :cite:`Tuckerman1992`.  This allows some forces in the system to be evaluated more
frequently than others.  For details on how to use it, consult the API
documentation.

Compound Integrator
-------------------

The :class:`CompoundIntegrator` class is useful for cases where you want to use
multiple integration algorithms within a single simulation.  It allows you to
create multiple integrators, then switch back and forth between them.  For
details on how to use it, consult the API documentation.

Temperature Coupling
====================


If you want to run a simulation at constant temperature, using a Langevin
integrator (as shown in the examples above) is usually the best way to do it.
OpenMM does provide an alternative, however: you can use a Verlet integrator,
then add an Andersen thermostat to your system to provide temperature coupling.

To do this, we can add an :class:`AndersenThermostat` object to the :class:`System` as shown below.
::

    ...
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
            constraints=HBonds)
    system.addForce(AndersenThermostat(300*kelvin, 1/picosecond))
    integrator = VerletIntegrator(0.002*picoseconds)
    ...

The two parameters of the Andersen thermostat are the temperature (300 K) and
collision frequency (1 ps\ :sup:`-1`\ ).

Pressure Coupling
=================


All the examples so far have been constant volume simulations.  If you want to
run at constant pressure instead, add a Monte Carlo barostat to your system.
You do this exactly the same way you added the Andersen thermostat in the
previous section:
::

    ...
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
            constraints=HBonds)
    system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    ...

The parameters of the Monte Carlo barostat are the pressure (1 bar) and
temperature (300 K).  The barostat assumes the simulation is being run at
constant temperature, but it does not itself do anything to regulate the
temperature.

.. warning::

    It is therefore critical that you always use it along with a Langevin integrator or
    Andersen thermostat, and that you specify the same temperature for both the barostat
    and the integrator or thermostat.  Otherwise, you will get incorrect results.

There also is an anisotropic barostat that scales each axis of the periodic box
independently, allowing it to change shape.  When using the anisotropic
barostat, you can specify a different pressure for each axis.  The following
line applies a pressure of 1 bar along the X and Y axes, but a pressure of 2 bar
along the Z axis:
::

    system.addForce(MonteCarloAnisotropicBarostat((1, 1, 2)*bar, 300*kelvin))

Another feature of the anisotropic barostat is that it can be applied to only
certain axes of the periodic box, keeping the size of the other axes fixed.
This is done by passing three additional parameters that specify whether the
barostat should be applied to each axis.  The following line specifies that the
X and Z axes of the periodic box should not be scaled, so only the Y axis can
change size.
::

    system.addForce(MonteCarloAnisotropicBarostat((1, 1, 1)*bar, 300*kelvin,
            False, True, False))

There is a third barostat designed specifically for simulations of membranes.
It assumes the membrane lies in the XY plane, and treats the X and Y axes of the
box differently from the Z axis.  It also applies a uniform surface tension in
the plane of the membrane.  The following line adds a membrane barostat that
applies a pressure of 1 bar and a surface tension of 200 bar*nm.  It specifies
that the X and Y axes are treated isotropically while the Z axis is free to
change independently.
::

    system.addForce(MonteCarloMembraneBarostat(1*bar, 200*bar*nanometer,
        MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree, 300*kelvin))

See the API documentation for details about the allowed parameter values and
their meanings.


Energy Minimization
===================


As seen in the examples, performing a local energy minimization takes a single
line in the script:
::

    simulation.minimizeEnergy()

In most cases, that is all you need.  There are two optional parameters you can
specify if you want further control over the minimization.  First, you can
specify a tolerance for when the energy should be considered to have converged:
::

    simulation.minimizeEnergy(tolerance=5*kilojoule/mole)

If you do not specify this parameter, a default tolerance of 10 kJ/mole is used.

Second, you can specify a maximum number of iterations:
::

    simulation.minimizeEnergy(maxIterations=100)

The minimizer will exit once the specified number of iterations is reached, even
if the energy has not yet converged.  If you do not specify this parameter, the
minimizer will continue until convergence is reached, no matter how many
iterations it takes.

These options are independent.  You can specify both if you want:
::

    simulation.minimizeEnergy(tolerance=0.1*kilojoule/mole, maxIterations=500)

Removing Center of Mass Motion
==============================


By default, :class:`System` objects created with the OpenMM application tools add
a :class:`CMMotionRemover` that removes all center of mass motion at every time step so the
system as a whole does not drift with time.  This is almost always what you
want.  In rare situations, you may want to allow the system to drift with time.
You can do this by specifying the :code:`removeCMMotion` parameter when you
create the System:
::

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff,
            removeCMMotion=False)

Writing Trajectories
====================


OpenMM can save simulation trajectories to disk in three formats: PDB_,
`PDBx/mmCIF`_, and DCD_.  All of these are widely supported formats, so you
should be able to read them into most analysis and visualization programs.

.. _PDB: http://www.wwpdb.org/documentation/format33/v3.3.html
.. _PDBx/mmCIF: http://mmcif.wwpdb.org
.. _DCD: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html

To save a trajectory, just add a “reporter” to the simulation, as shown in the
example scripts above:
::

    simulation.reporters.append(PDBReporter('output.pdb', 1000))

The two parameters of the :class:`PDBReporter` are the output filename and how often (in
number of time steps) output structures should be written.  To use PDBx/mmCIF or
DCD format, just replace :class:`PDBReporter` with :class:`PDBxReporter` or
:class:`DCDReporter`.  The parameters represent the same values:
::

    simulation.reporters.append(DCDReporter('output.dcd', 1000))

Recording Other Data
====================


In addition to saving a trajectory, you may want to record other information
over the course of a simulation, such as the potential energy or temperature.
OpenMM provides a reporter for this purpose also.  Create a :class:`StateDataReporter`
and add it to the simulation:
::

    simulation.reporters.append(StateDataReporter('data.csv', 1000, time=True,
            kineticEnergy=True, potentialEnergy=True))

The first two parameters are the output filename and how often (in number of
time steps) values should be written.  The remaining arguments specify what
values should be written at each report.  The available options are
:code:`step` (the index of the current time step), :code:`time`\ ,
:code:`progress` (what percentage of the simulation has completed),
:code:`remainingTime` (an estimate of how long it will take the simulation to
complete), :code:`potentialEnergy`\ , :code:`kineticEnergy`\ ,
:code:`totalEnergy`\ , :code:`temperature`\ , :code:`volume` (the volume
of the periodic box), :code:`density` (the total system mass divided by the
volume of the periodic box), and :code:`speed` (an estimate of how quickly
the simulation is running).  If you include either the :code:`progress` or
:code:`remainingTime` option, you must also include the :code:`totalSteps`
parameter to specify the total number of time steps that will be included in the
simulation.  One line is written to the file for each report containing the
requested values.  By default the values are written in comma-separated-value
(CSV) format.  You can use the :code:`separator` parameter to choose a
different separator.  For example, the following line will cause values to be
separated by spaces instead of commas:
::

    simulation.reporters.append(StateDataReporter('data.txt', 1000, progress=True,
            temperature=True, totalSteps=10000, separator=' '))


Saving Simulation Progress and Results
==========================================

There are three built-in ways to save the results of your simulation in OpenMM
(additional methods can be written yourself or imported through other packages
like mdtraj or parmed). If you are simply interested in saving the structure,
you can write it out as a PDB file using :code:`PDBFile.writeFile()`.  You can
see an example of this in the modeller section :ref:`saving-the-results`.

If you are hoping to save more information than just positions, you can use
:code:`simulation.saveState()`. This will save the entire state of the
simulation, including positions, velocities, box dimensions and much more in an
XML file. This same file can be loaded back into OpenMM and used to continue
the simulation. Importantly, because this file is a text file, it can be
transfered between different platforms and different versions of OpenMM. A
potential downside to this approach is that state files are often quite large,
and may not fit all use cases. Here's an example of how to use it:
::

    simulation.saveState('output.xml')

To load the simulation back in:
::

    simulation.loadState('output.xml')

There is a third way to save your simulation, known as a checkpoint file, which
will save the entire simulation as a binary file. It will allow you to exactly
continue a simulation if the need arises (though whether the simulation is
deterministic depends on platform and methods, see
:ref:`platform-specific-properties-determinism`). There are important caveats
to this approach, however. This binary can only be used to restart simulations
on machines with the same hardware and the same OpenMM version as the one that
saved it. Therefore, it should only be used when it's clear that won't be an
issue.
::

    simulation.saveCheckpoint('state.chk')

And can be loaded back in like this:
::

    simulation.loadCheckpoint('state.chk')

Finally, OpenMM comes with a built-in reporter for saving checkpoints, the
:class:`CheckpointReporter`, which can be helpful in restarting simulations
that failed unexpectedly or due to outside reasons (e.g. server crash). To save
a checkpoint file every 5,000 steps, for example:
::

    simulation.reporters.append(CheckpointReporter('checkpnt.chk', 5000))

Note that the checkpoint reporter will overwrite the last checkpoint file.


Enhanced Sampling Methods
=========================

In many situations, the goal of a simulation is to sample the range of configurations
accessible to a system.  It does not matter whether the simulation represents a
single, physically realistic trajectory, only whether it produces a correct distribution
of states.  In this case, a variety of methods can be used to sample configuration
space much more quickly and efficiently than a single physical trajectory would.
These are known as enhanced sampling methods.  OpenMM offers several that you
can choose from.  They are briefly described here.  Consult the API documentation
for more detailed descriptions and example code.

Simulated Tempering
-------------------

Simulated tempering\ :cite:`Marinari1992` involves making frequent changes to the
temperature of a simulation.  At high temperatures, it can quickly cross energy barriers
and explore a wide range of configurations.  At lower temperatures, it more thoroughly
explores each local region of configuration space.  This is a powerful method to
speed up sampling when you do not know in advance what motions you want to sample.
Simply specify the range of temperatures to simulate and the algorithm handles
everything for you mostly automatically.

Metadynamics
------------

Metadynamics\ :cite:`Barducci2008` is used when you do know in advance what
motions you want to sample.  You specify one or more collective variables, and the
algorithm adds a biasing potential to make the simulation explore a wide range of
values for those variables.  It does this by periodically adding "bumps" to the biasing
potential at the current values of the collective variables.  This encourages the simulation
to move away from regions it has already explored and sample a wide range of values.
At the end of the simulation, the biasing potential can be used to calculate the
free energy of the system as a function of the collective variables.

Accelerated Molecular Dynamics (aMD)
------------------------------------

aMD\ :cite:`Hamelberg2007` is another method that can be used when you do not know in
advance what motions you want to accelerate.  It alters the potential energy surface
by adding a "boost" potential whenever the potential energy is below a threshold.
This makes local minima shallower and allows more frequent transitions between them.
The boost can be applied to the total potential energy, to just a subset of interactions
(typically the dihedral torsions), or both.  There are separate integrator classes
for each of these options: :class:`AMDIntegrator`, :class:`AMDForceGroupIntegrator`,
and :class:`DualAMDIntegrator`.


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

Allowed values for the :code:`model` option are ``'tip3p'``, ``'tip3pfb'``, ``'spce'``,
``'tip4pew'``, ``'tip4pfb'``, and ``'tip5p'``.  Be sure to include the single quotes
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

        from simtk.openmm.app import *
        from simtk.openmm import *
        from simtk.unit import *

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
        positions = simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(simulation.topology, positions, open('output.pdb', 'w'))
        print('Done')

    .. caption::

        :autonumber:`Example,Modeller complete`


Advanced Simulation Examples
############################

In the previous chapter, we looked at some basic scripts for running simulations
and saw lots of ways to customize them.  If that is all you want to do—run
straightforward molecular simulations—you already know everything you need to
know.  Just use the example scripts and customize them in the ways described in
Section :ref:`simulation-parameters`.

OpenMM can do far more than that.  Your script has the full OpenMM API at its
disposal, along with all the power of the Python language and libraries.  In
this chapter, we will consider some examples that illustrate more advanced
techniques.  Remember that these are still only examples; it would be impossible
to give an exhaustive list of everything OpenMM can do.  Hopefully they will
give you a sense of what is possible, and inspire you to experiment further on
your own.

Starting in this section, we will assume some knowledge of programming, as well
as familiarity with the OpenMM API.  Consult the OpenMM Users Guide and API
documentation if you are uncertain about how something works.   You can also use
the Python :code:`help` command.  For example,
::

    help(Simulation)

will print detailed documentation on the :class:`Simulation` class.

Simulated Annealing
*******************

Here is a very simple example of how to do simulated annealing.  The following
lines linearly reduce the temperature from 300 K to 0 K in 100 increments,
executing 1000 time steps at each temperature:

.. samepage::
    ::

        ...
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        for i in range(100):
            integrator.setTemperature(3*(100-i)*kelvin)
            simulation.step(1000)

    .. caption::

        :autonumber:`Example,simulated annealing`

This code needs very little explanation.  The loop is executed 100 times.  Each
time through, it adjusts the temperature of the :class:`LangevinIntegrator` and then
calls :code:`step(1000)` to take 1000 time steps.

Applying an External Force to Particles: a Spherical Container
**************************************************************

In this example, we will simulate a non-periodic system contained inside a
spherical container with radius 2 nm.  We implement the container by applying a
harmonic potential to every particle:

.. math::
    E(r) = \begin{cases}
           0          & r\le2\\
           100(r-2)^2 & r>2
           \end{cases}

where *r* is the distance of the particle from the origin, measured in nm.
We can easily do this using OpenMM’s :class:`CustomExternalForce` class.  This class
applies a force to some or all of the particles in the system, where the energy
is an arbitrary function of each particle’s (\ *x*\ , *y*\ , *z*\ )
coordinates.  Here is the code to do it:

.. samepage::
    ::

        ...
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic,
                nonbondedCutoff=1*nanometer, constraints=None)
        force = CustomExternalForce('100*max(0, r-2)^2; r=sqrt(x*x+y*y+z*z)')
        system.addForce(force)
        for i in range(system.getNumParticles()):
            force.addParticle(i, [])
        integrator = LangevinIntegrator(300*kelvin, 91/picosecond, 0.002*picoseconds)
        ...

    .. caption::

        :autonumber:`Example,spherical container`

The first thing it does is create a :class:`CustomExternalForce` object and add it to the
:class:`System`.  The argument to :class:`CustomExternalForce` is a mathematical expression
specifying the potential energy of each particle.  This can be any function of *x*\ ,
*y*\ , and *z* you want.  It also can depend on global or per-particle
parameters.  A wide variety of restraints, steering forces, shearing forces,
etc. can be implemented with this method.

Next it must specify which particles to apply the force to.  In this case, we
want it to affect every particle in the system, so we loop over them and call
:meth:`addParticle` once for each one.  The two arguments are the index of
the particle to affect, and the list of per-particle parameter values (an empty
list in this case).  If we had per-particle parameters, such as to make the
force stronger for some particles than for others, this is where we would
specify them.

Notice that we do all of this immediately after creating the :class:`System`.  That is
not an arbitrary choice.

.. warning::

    If you add new forces to a :class:`System`, you must do so before creating the :class:`Simulation`.
    Once you create a :class:`Simulation`, modifying the :class:`System` will have no effect on that :class:`Simulation`.

Extracting and Reporting Forces (and other data)
************************************************

OpenMM provides reporters for three output formats: PDB_, `PDBx/mmCIF`_ and DCD_.
All of those formats store only positions, not velocities, forces, or other data.  In this
section, we create a new reporter that outputs forces.  This illustrates two
important things: how to write a reporter, and how to query the simulation for
forces or other data.

Here is the definition of the :class:`ForceReporter` class:

.. samepage::
    ::

        class ForceReporter(object):
            def __init__(self, file, reportInterval):
                self._out = open(file, 'w')
                self._reportInterval = reportInterval

            def __del__(self):
                self._out.close()

            def describeNextReport(self, simulation):
                steps = self._reportInterval - simulation.currentStep%self._reportInterval
                return (steps, False, False, True, False, None)

            def report(self, simulation, state):
                forces = state.getForces().value_in_unit(kilojoules/mole/nanometer)
                for f in forces:
                    self._out.write('%g %g %g\n' % (f[0], f[1], f[2]))

    .. caption::

        :autonumber:`Example,ForceReporter`

The constructor and destructor are straightforward.  The arguments to the
constructor are the output filename and the interval (in time steps) at which it
should generate reports.  It opens the output file for writing and records the
reporting interval.  The destructor closes the file.

We then have two methods that every reporter must implement:
:meth:`describeNextReport()` and :meth:`report()`.  A Simulation object
periodically calls :meth:`describeNextReport()` on each of its reporters to
find out when that reporter will next generate a report, and what information
will be needed to generate it.  The return value should be a six element tuple,
whose elements are as follows:

* The number of time steps until the next report.  We calculate this as
  *(report interval)*\ -\ *(current step)*\ %\ *(report interval)*\ .  For
  example, if we want a report every 100 steps and the simulation is currently on
  step 530, we will return 100-(530%100) = 70.
* Whether the next report will need particle positions.
* Whether the next report will need particle velocities.
* Whether the next report will need forces.
* Whether the next report will need energies.
* Whether the positions should be wrapped to the periodic box.  If None, it will
  automatically decide whether to wrap positions based on whether the System uses
  periodic boundary conditions.


When the time comes for the next scheduled report, the :class:`Simulation` calls
:meth:`report()` to generate the report.  The arguments are the :class:`Simulation`
object, and a :class:`State` that is guaranteed to contain all the information that was
requested by :meth:`describeNextReport()`\ .  A State object contains a
snapshot of information about the simulation, such as forces or particle
positions.  We call :meth:`getForces()` to retrieve the forces and convert
them to the units we want to output (kJ/mole/nm).  Then we loop over each value
and write it to the file.  To keep the example simple, we just print the values
in text format, one line per particle.  In a real program, you might choose a
different output format.

Now that we have defined this class, we can use it exactly like any other
reporter.  For example,
::

    simulation.reporters.append(ForceReporter('forces.txt', 100))

will output forces to a file called “forces.txt” every 100 time steps.

Computing Energies
******************

This example illustrates a different sort of analysis.  Instead of running a
simulation, assume we have already identified a set of structures we are
interested in.  These structures are saved in a set of PDB files.  We want to
loop over all the files in a directory, load them in one at a time, and compute
the potential energy of each one.  Assume we have already created our :class:`System` and
:class:`Simulation`.  The following lines perform the analysis:

.. samepage::
    ::

        import os
        for file in os.listdir('structures'):
            pdb = PDBFile(os.path.join('structures', file))
            simulation.context.setPositions(pdb.positions)
            state = simulation.context.getState(getEnergy=True)
            print(file, state.getPotentialEnergy())

    .. caption::

        :autonumber:`Example,computing energies`

We use Python’s :code:`listdir()` function to list all the files in the
directory.  We create a :class:`PDBFile` object for each one and call
:meth:`setPositions()` on the Context to specify the particle positions loaded
from the PDB file.  We then compute the energy by calling :meth:`getState()`
with the option :code:`getEnergy=True`\ , and print it to the console along
with the name of the file.


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
OutOfPlaneSite, and LocalCoordinatesSite classes respectively.  The :code:`siteName`
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
   exact matching to replicate the beheavior of `SMIRNOFF <https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html>`_
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

CustomNonbondedForce also allows you to define tabulated functions.  See section
:ref:`tabulated-functions` for details.

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

CustomGBForce also allows you to define tabulated functions.  See section
:ref:`tabulated-functions` for details.

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

CustomHbondForce also allows you to define tabulated functions.  See section
:ref:`tabulated-functions` for details.

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

CustomManyParticleForce also allows you to define tabulated functions.  See section
:ref:`tabulated-functions` for details.

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
definitions are combined as follows.

* A file may refer to atom types and classes that it defines, as well as those
  defined in previous files.  It may not refer to ones defined in later files.
  This means that the order in which files are listed when calling the ForceField
  constructor is potentially significant.
* Forces that involve per-atom parameters (such as NonbondedForce or
  GBSAOBCForce) require parameter values to be defined for every atom type.  It
  does not matter which file those types are defined in.  For example, files that
  define explicit water models generally define a small number of atom types, as
  well as nonbonded parameters for those types.  In contrast, files that define
  implicit solvent models do not define any new atom types, but provide parameters
  for all the atom types that were defined in the main force field file.
* For other forces, the files are effectively independent.  For example, if two
  files each include a :code:`<HarmonicBondForce>` tag, bonds will be created
  based on the rules in the first file, and then more bonds will be created based
  on the rules in the second file.  This means you could potentially end up with
  multiple bonds between a single pair of atoms.


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
        forcefield : simtk.openmm.app.ForceField
            The ForceField object to which residue templates and/or parameters are to be added.
        residue : simtk.openmm.app.Topology.Residue
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

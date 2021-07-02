.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++

.. _the-opencl-platform:

The OpenCL Platform
###################

The OpenCL Platform is much more complicated than the reference Platform.  It
also provides many more tools to simplify your work, but those tools themselves
can be complicated to use correctly.  This chapter will attempt to explain how
to use some of the most important ones.  It will *not* teach you how to
program with OpenCL.  There are many tutorials on that subject available
elsewhere, and this guide assumes you already understand it.

Overview
********

When using the OpenCL Platform, the “platform-specific data” stored in
ContextImpl is of type OpenCLPlatform::PlatformData, which is declared in
OpenCLPlatform.h.  The most important field of this class is :code:`contexts`
, which is a vector of OpenCLContexts.  (There is one OpenCLContext for each
device you are using.  The most common case is that you are running everything
on a single device, in which case there will be only one OpenCLContext.
Parallelizing computations across multiple devices is not discussed here.)  The
OpenCLContext stores most of the important information about a simulation:
positions, velocities, forces, an OpenCL CommandQueue used for executing
kernels, workspace buffers of various sorts, etc.  It provides many useful
methods for compiling and executing kernels, clearing and reducing buffers, and
so on.  It also provides access to three other important objects: the
OpenCLIntegrationUtilities, OpenCLNonbondedUtilities, and OpenCLBondedUtilities.
These are discussed below.

Allocation of device memory is generally done through the OpenCLArray class.  It
takes care of much of the work of memory management, and provides a simple
interface for transferring data between host and device memory.

Every kernel is specific to a particular OpenCLContext, which in turn is
specific to a particular OpenMM::Context.  This means that kernel source code
can be customized for a particular simulation.  For example, values such as the
number of particles can be turned into compile-time constants, and specific
versions of kernels can be selected based on the device being used or on
particular aspects of the system being simulated.
:code:`OpenCLContext::createProgram()` makes it easy to specify a list of
preprocessor definitions to use when compiling a kernel.

The normal way to execute a kernel is by calling :code:`executeKernel()` on
the OpenCLContext.  It allows you to specify the total number of work-items to
execute, and optionally the size of each work-group.  (If you do not specify a
work-group size, it uses 64 as a default.)  The number of work-groups to launch
is selected automatically based on the work-group size, the total number of
work-items, and the number of compute units in the device it will execute on.

Numerical Precision
*******************

The OpenCL platform supports three precision modes:

#. **Single**\ : All values are stored in single precision, and nearly all
   calculations are done in single precision.  The arrays of positions, velocities,
   forces, and energies (returned by the OpenCLContext’s :code:`getPosq()`\ ,
   :code:`getVelm()`\ , :code:`getForce()`\ , :code:`getForceBuffers()`\ , and
   :code:`getEnergyBuffer()` methods) are all of type :code:`float4` (or
   :code:`float` in the case of :code:`getEnergyBuffer()`\ ).
#. **Mixed**\ : Forces are computed and stored in single precision, but
   integration is done in double precision.  The velocities have type
   :code:`double4`\ .  The positions are still stored in single precision to avoid
   adding overhead to the force calculations, but a second array of type
   :code:`float4` is created to store “corrections” to the positions (returned by
   the OpenCLContext’s getPosqCorrection() method).  Adding the position and the
   correction together gives the full double precision position.
#. **Double**\ : Positions, velocities, forces, and energies are all stored in
   double precision, and nearly all calculations are done in double precision.


You can call :code:`getUseMixedPrecision()` and
:code:`getUseDoublePrecision()` on the OpenCLContext to determine which mode
is being used.  In addition, when you compile a kernel by calling
:code:`createKernel()`\ , it automatically defines two types for you to make it
easier to write kernels that work in any mode:

#. :code:`real` is defined as :code:`float` in single or mixed precision
   mode, :code:`double` in double precision mode.
#. :code:`mixed` is defined as :code:`float` in single precision mode,
   :code:`double` in mixed or double precision mode.


It also defines vector versions of these types (\ :code:`real2`\ ,
:code:`real4`\ , etc.).

.. _computing-forces:

Computing Forces
****************

When forces are computed, they can be stored in either of two places.  There is
an array of :code:`long` values storing them as 64 bit fixed point values, and
a collection of buffers of :code:`real4` values storing them in floating point
format.  Most GPUs support atomic operations on 64 bit integers, which allows
many threads to simultaneously record forces without a danger of conflicts.
Some low end GPUs do not support this, however, especially the embedded GPUs
found in many laptops.  These devices write to the floating point buffers, with
careful coordination to make sure two threads will never write to the same
memory location at the same time.

At the start of a force calculation, all forces in all buffers are set to zero.
Each Force is then free to add its contributions to any or all of the buffers.
Finally, the buffers are summed to produce the total force on each particle.
The total is recorded in both the floating point and fixed point arrays.

The size of each floating point buffer is equal to the number of particles, rounded up to the
next multiple of 32.  Call :code:`getPaddedNumAtoms()` on the OpenCLContext
to get that number.  The actual force buffers are obtained by calling
:code:`getForceBuffers()`\ .  The first *n* entries (where *n* is the
padded number of atoms) represent the first force buffer, the next *n*
represent the second force buffer, and so on.  More generally, the *i*\ ’th
force buffer’s contribution to the force on particle *j* is stored in
element :code:`i*context.getPaddedNumAtoms()+j`\ .

The fixed point buffer is ordered differently.  For atom *i*\ , the x component
of its force is stored in element :code:`i`\ , the y component in element
:code:`i+context.getPaddedNumAtoms()`\ , and the z component in element
:code:`i+2*context.getPaddedNumAtoms()`\ .  To convert a value from floating
point to fixed point, multiply it by 0x100000000 (2\ :sup:`32`\ ),
then cast it to a :code:`long`\ .  Call :code:`getLongForceBuffer()` to get the
array of fixed point values.

The potential energy is also accumulated in a set of buffers, but this one is
simply a list of floating point values.  All of them are set to zero at the
start of a computation, and they are summed at the end of the computation to
yield the total energy.

The OpenCL implementation of each Force object should define a subclass of
ComputeForceInfo, and register an instance of it by calling :code:`addForce()` on
the OpenCLContext.  It implements methods for determining whether particular
particles or groups of particles are identical.  This is important when
reordering particles, and is discussed below.


Nonbonded Forces
****************

Computing nonbonded interactions efficiently is a complicated business in the
best of cases.  It is even more complicated on a GPU.  Furthermore, the
algorithms must vary based on the type of processor being used, whether there is
a distance cutoff, and whether periodic boundary conditions are being applied.

The OpenCLNonbondedUtilities class tries to simplify all of this.  To use it you
need provide only a piece of code to compute the interaction between two
particles.  It then takes responsibility for generating a neighbor list, looping
over interacting particles, loading particle parameters from global memory, and
writing the forces and energies to the appropriate buffers.  All of these things
are done using an algorithm appropriate to the processor you are running on and
high level aspects of the interaction, such as whether it uses a cutoff and
whether particular particle pairs need to be excluded.

Of course, this system relies on certain assumptions, the most important of
which is that the Force can be represented as a sum of independent pairwise
interactions.  If that is not the case, things become much more complicated.
You may still be able to use features of OpenCLNonbondedUtilities, but you
cannot use the simple mechanism outlined above.  That is beyond the scope of
this guide.

To define a nonbonded interaction, call :code:`addInteraction()` on the
OpenCLNonbondedUtilities, providing a block of OpenCL source code for computing
the interaction.  This block of source code will be inserted into the middle of
an appropriate kernel.  At the point where it is inserted, various variables
will have been defined describing the interaction to compute:

#. :code:`atom1` and :code:`atom2` are the indices of the two
   interacting particles.
#. :code:`r`\ , :code:`r2`\ , and :code:`invR` are the distance *r*
   between the two particles, *r*\ :sup:`2`\ , and 1/\ *r* respectively.
#. :code:`isExcluded` is a :code:`bool` specifying whether this pair of
   particles is marked as an excluded interaction.  (Excluded pairs are not skipped
   automatically, because in some cases they still need to be processed, just
   differently from other pairs.)
#. :code:`posq1` and :code:`posq2` are :code:`real4`\ s containing the
   positions (in the xyz fields) and charges (in the w fields) of the two
   particles.
#. Other per-particle parameters may be specified, as described below.


The following preprocessor macros will also have been defined:

#. :code:`NUM_ATOMS` is the total number of particles in the system.
#. :code:`PADDED_NUM_ATOMS` is the padded number of particles in the system.
#. :code:`USE_CUTOFF` is defined if and only if a cutoff is being used
#. :code:`USE_PERIODIC` is defined if and only if periodic boundary
   conditions are being used.
#. :code:`CUTOFF` and :code:`CUTOFF_SQUARED` are the cutoff distance and
   its square respectively (but only defined if a cutoff is being used).


Finally, two output variables will have been defined:

#. You should add the energy of the interaction to :code:`tempEnergy`\ .
#. You should add the derivative of the energy with respect to the inter-particle
   distance to :code:`dEdR`\ .


You can also define arbitrary per-particle parameters by calling
:code:`addParameter()` on the OpenCLNonbondedUtilities.  You provide an array
in device memory containing the set of values, and the values for the two
interacting particles will be loaded and stored into variables called
:code:`<name>1` and :code:`<name>2`\ , where <name> is the name you specify
for the parameter.  Note that nonbonded interactions are not computed until
after :code:`calcForcesAndEnergy()` has been called on every ForceImpl, so
it is possible to make the parameter values change with time by modifying them
inside :code:`calcForcesAndEnergy()`\ .  Also note that the length of the
array containing the parameter values must equal the *padded* number of
particles in the system.

Finally, you can specify arbitrary other memory objects that should be passed as
arguments to the interaction kernel by calling :code:`addArgument()`\ .  The
rest of the kernel ignores these arguments, but you can make use of them in your
interaction code.

Consider a simple example.  Suppose we want to implement a nonbonded interaction
of the form *E*\ =\ *k*\ :sub:`1`\ *k*\ :sub:`2`\ *r*\ :sup:`2`\ ,
where *k* is a per-particle parameter.  First we create a parameter as
follows
::

    nb.addParameter(ComputeParameterInfo(kparam, "kparam", "float", 1));

where :code:`nb` is the OpenCLNonbondedUtilities for the context.  Now we
call :code:`addInteraction()` to define an interaction with the following
source code:
::

    #ifdef USE_CUTOFF
    if (!isExcluded && r2 < CUTOFF_SQUARED) {
    #else
    if (!isExcluded) {
    #endif
        tempEnergy += kparam1*kparam2*r2;
        dEdR += 2*kparam1*kparam2*r;
    }

An important point is that this code is executed for every pair of particles in
the *padded* list of atoms.  This means that some interactions involve
padding atoms, and should not actually be included.  You might think, then, that
the above code is incorrect and we need another check to filter out the extra
interactions:
::

    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)

This is not necessary in our case, because the :code:`isExcluded` flag is
always set for interactions that involve a padding atom.  If our force did not
use excluded interactions (and so did not check :code:`isExcluded`\ ), then we
would need to add this extra check.  Self interactions are a similar case: we do
not check for :code:`(atom1 == atom2)` because the exclusion flag prevents
them from being processed, but for some forces that check is necessary.

Bonded Forces
*************

Just as OpenCLNonbondedUtilities simplifies the task of creating nonbonded
interactions, OpenCLBondedUtilities simplifies the process for many types of
bonded interactions.  A “bonded interaction” means one that is applied to small,
fixed groups of particles.  This includes bonds, angles, torsions, etc.  The
important point is that the list of particles forming a “bond” is known in
advance and does not change with time.

Using OpenCLBondedUtilities is very similar to the process described above.  You
provide a block of OpenCL code for evaluating a single interaction.  This block
of code will be inserted into the middle of a kernel that loops over all
interactions and evaluates each one.  At the point where it is inserted, the
following variables will have been defined describing the interaction to
compute:

#. :code:`index` is the index of the interaction being evaluated.
#. :code:`atom1`\ , :code:`atom2`\ , ... are the indices of the interacting
   particles.
#. :code:`pos1`\ , :code:`pos2`\ , ... are :code:`real4`\ s containing the
   positions (in the xyz fields) of the interacting particles.


A variable called :code:`energy` will have been defined for accumulating the
total energy of all interactions.  Your code should add the energy of the
interaction to it.  You also should define :code:`real4` variables called
:code:`force1`\ , :code:`force2`\ , ... and store the force on each atom into
them.

As a simple example, the following source code implements a pairwise interaction
of the form *E*\ =\ *r*\ :sup:`2`\ :
::

    real4 delta = pos2-pos1;
    energy += delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
    real4 force1 = 2.0f*delta;
    real4 force2 = -2.0f*delta;

To use it, call :code:`addInteraction()` on the Context’s
OpenCLBondedUtilities object.  You also provide a list of the particles involved
in every bonded interaction.

Exactly as with nonbonded interactions, you can call :code:`addArgument()`
to specify arbitrary memory objects that should be passed as arguments to the
interaction kernel.  These might contain per-bond parameters (use
:code:`index` to look up the appropriate element) or any other information you
want.

Reordering of Particles
***********************

Nonbonded calculations are done a bit differently in the OpenCL Platform than in
most CPU based codes.  In particular, interactions are computed on blocks of 32
particles at a time (which is why the number of particles needs to be padded to
bring it up to a multiple of 32), and the neighbor list actually lists pairs of
\ *blocks*\ , not pairs of individual particles, that are close enough to
interact with each other.

This only works well if sequential particles tend to be close together so that
blocks are spatially compact.  This is generally true of particles in a
macromolecule, but it is not true for solvent molecules.  Each water molecule,
for example, can move independently of other water molecules, so particles that
happen to be sequential in whatever order the molecules were defined in need not
be spatially close together.

The OpenCL Platform addresses this by periodically reordering particles so that
sequential particles are close together.  This means that what the OpenCL
Platform calls particle *i* need not be the same as what the System calls
particle *i*\ .

This reordering is done frequently, so it must be very fast.  If all the data
structures describing the structure of the System and the Forces acting on it
needed to be updated, that would make it prohibitively slow.  The OpenCL
Platform therefore only reorders particles in ways that do not alter any part of
the System definition.  In practice, this means exchanging entire molecules; as
long as two molecules are truly identical, their positions and velocities can be
exchanged without affecting the System in any way.

Every Force can contribute to defining the boundaries of molecules, and to
determining whether two molecules are identical.  This is done through the
ComputeForceInfo it adds to the OpenCLContext.  It can specify two types of
information:

#. Given a pair of particles, it can say whether those two particles are
   identical (as far as that Force is concerned).  For example, a Force object
   implementing a Coulomb force would check whether the two particles had equal
   charges.
#. It can define *particle groups*\ .  The OpenCL Platform will ensure that
   all the particles in a group are part of the same molecule.  It also can specify
   whether two groups are identical to each other.  For example, in a Force
   implementing harmonic bonds, each group would consist of the two particles
   connected by a bond, and two groups would be identical if they had the same
   spring constants and equilibrium lengths.


Integration Utilities
*********************

The OpenCLContext’s OpenCLIntegrationUtilities provides features that are used
by many integrators.  The two most important are random number generation and
constraint enforcement.

If you plan to use random numbers, you should call
:code:`initRandomNumberGenerator()` during initialization, specifying the
random number seed to use.  Be aware that there is only one random number
generator, even if multiple classes make use of it.  If two classes each call
:code:`initRandomNumberGenerator()` and request different seeds, an exception
will be thrown.  If they each request the same seed, the second call will simply
be ignored.

For efficiency, random numbers are generated in bulk and stored in an array in
device memory, which you can access by calling :code:`getRandom()`\ .  Each
time you need to use a block of random numbers, call
:code:`prepareRandomNumbers()`\ , specifying how many values you need.  It will
register that many values as having been used, and return the index in the array
at which you should start reading values.  If not enough unused values remain in
the array, it will generate a new batch of random values before returning.

To apply constraints, simply call :code:`applyConstraints()`\ .  For numerical
accuracy, the constraint algorithms do not work on particle positions directly,
but rather on the *displacements* taken by the most recent integration step.
These displacements must be stored in an array which you can get by calling
:code:`getPosDelta()`\ .  That is, the constraint algorithms assume the actual
(unconstrained) position of each particle equals the position stored in the
OpenCLContext plus the delta stored in the OpenCLIntegrationUtilities.  It then
modifies the deltas so that all distance constraints are satisfied.  The
integrator must then finish the time step by adding the deltas to the positions
and storing them into the main position array.

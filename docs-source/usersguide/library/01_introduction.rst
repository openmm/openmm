.. _the-openmm-library-introduction:

Introduction
############


What Is the OpenMM Library?
***************************

OpenMM consists of two parts.  First, there is a set of libraries for performing
many types of computations needed for molecular simulations: force evaluation,
numerical integration, energy minimization, etc.  These libraries provide an
interface targeted at developers of simulation software, allowing them to easily
add simulation features to their programs.

Second, there is an “application layer”, a set of Python libraries providing a
high level interface for running simulations.  This layer is targeted at
computational biologists or other people who want to run simulations, and who
may or may not be programmers.

The first part of this guide focused on the application layer and described how to run
simulations with it.  We now turn to the lower level libraries.  We will assume
you are a programmer, that you are writing your own applications, and that you
want to add simulation features to those applications.  The following chapters
describe how to do that with OpenMM.

How to get started
==================

We have provided a number of files that make it easy to get started with OpenMM.
Pre-compiled binaries are provided for quickly getting OpenMM onto your computer
(See Chapter :numref:`installing-openmm` for set-up instructions).  We recommend that you then
compile and run some of the tutorial examples, described in Chapter :numref:`openmm-tutorials`.
These highlight key functions within OpenMM and teach you the basic programming concepts for using
OpenMM.  Once you are ready to begin integrating OpenMM into a specific software package, read
through Chapter :numref:`examples-of-openmm-integration` to see how other software developers have
done this.

License
========

Two different licenses are used for different parts of OpenMM.  The public API,
the low level API, the reference platform, the CPU platform, and the application
layer are all distributed under the MIT
license.  This is a very permissive license which allows them to be used in
almost any way, requiring only that you retain the copyright notice and
disclaimer when distributing them.

The CUDA, HIP, and OpenCL platforms are distributed under the GNU Lesser General
Public License (LGPL).  This also allows you to use, modify, and distribute them
in any way you want, but it requires you to also distribute the source code for
your modifications.  This restriction applies only to modifications to OpenMM
itself; you need not distribute the source code to applications that use it.

OpenMM also uses several pieces of code that were written by other people and
are covered by other licenses.  All of these licenses are similar in their terms
to the MIT license, and do not significantly restrict how OpenMM can be used.

All of these licenses may be found in the “licenses” directory included with
OpenMM.


Design Principles
*****************

The design of the OpenMM API is guided by the following principles.

1. The API must support efficient implementations on a variety of architectures.

The most important consequence of this goal is that the API cannot provide
direct access to state information (particle positions, velocities, etc.) at all
times.  On some architectures, accessing this information is expensive.  With a
GPU, for example, it will be stored in video memory, and must be transferred to
main memory before outside code can access it.  On a distributed architecture,
it might not even be present on the local computer.  OpenMM therefore only
allows state information to be accessed in bulk, with the understanding that
doing so may be a slow operation.

2. The API should be easy to understand and easy to use.

This seems obvious, but it is worth stating as an explicit goal.  We are
creating OpenMM with the hope that many other people will use it.  To achieve
that goal, it should be possible for someone to learn it without an enormous
amount of effort.  An equally important aspect of being “easy to use” is being
easy to use *correctly*\ .  A well designed API should minimize the
opportunities for a programmer to make mistakes.  For both of these reasons,
clarity and simplicity are essential.

3. It should be modular and extensible.

We cannot hope to provide every feature any user will ever want.  For that
reason, it is important that OpenMM be easy to extend.  If a user wants to add a
new molecular force field, a new thermostat algorithm, or a new hardware
platform, the API should make that easy to do.

4. The API should be hardware independent.

Computer architectures are changing rapidly, and it is impossible to predict
what hardware platforms might be important to support in the future.  One of the
goals of OpenMM is to separate the API from the hardware.  The developers of a
simulation application should be able to write their code once, and have it
automatically take advantage of any architecture that OpenMM supports, even
architectures that do not yet exist when they write it.

Choice of Language
******************

Molecular modeling and simulation tools are written in a variety of languages:
C, C++, Fortran, Python, TCL, etc.  It is important that any of these tools be
able to use OpenMM.  There are two possible approaches to achieving this goal.

One option is to provide a separate version of the API for each language.  These
could be created by hand, or generated automatically with a wrapper generator
such as SWIG.  This would require the API to use only “lowest common
denominator” features that can be reasonably supported in all languages.  For
example, an object oriented API would not be an option, since it could not be
cleanly expressed in C or Fortran.

The other option is to provide a single version of the API written in a single
language.  This would permit a cleaner, simpler API, but also restrict the
languages it could be directly called from.  For example, a C++ API could not be
invoked directly from Fortran or Python.

We have chosen to use a hybrid of these two approaches.  OpenMM is based on an
object oriented C++ API.  This is the primary way to invoke OpenMM, and is the
only API that fully exposes all features of the library.  We believe this will
ultimately produce the best, easiest to use API and create the least work for
developers who use it.  It does require that any code which directly invokes
this API must itself be written in C++, but this should not be a significant
burden.  Regardless of what language we had chosen, developers would need to
write a thin layer for translating between their own application’s data model
and OpenMM.  That layer is the only part which needs to be written in C++.

In addition, we have created wrapper APIs that allow OpenMM to be invoked from
other languages.  The current release includes wrappers for C, Fortran, and
Python.  These wrappers support as many features as reasonably possible given
the constraints of the particular languages, but some features cannot be fully
supported.  In particular, writing plug-ins to extend the OpenMM API can only be
done in C++.

We are also aware that some features of C++ can easily lead to compatibility and
portability problems, and we have tried to avoid those features.  In particular,
we make minimal use of templates and avoid multiple inheritance altogether.  Our
goal is to support OpenMM on all major compilers and operating systems.

Architectural Overview
**********************

OpenMM is based on a layered architecture, as shown in the following diagram:


.. figure:: ../../images/ArchitectureLayers.jpg
   :align: center
   :width: 100%

   :autonumber:`Figure,OpenMM architecture`:  OpenMM architecture

At the highest level is the OpenMM public API.  This is the API developers
program against when using OpenMM within their own applications.  It is designed
to be simple, easy to understand, and completely platform independent.  This is
the only layer that many users will ever need to look at.

The public API is implemented by a layer of platform independent code.  It
serves as the interface to the lower level, platform specific code.  Most users
will never need to look at it.

The next level down is the OpenMM Low Level API (OLLA).  This acts as an
abstraction layer to hide the details of each hardware platform.  It consists of
a set of C++ interfaces that each platform must implement.  Users who want to
extend OpenMM will need to write classes at the OLLA level.  Note the different
roles played by the public API and the low level API: the public API defines an
interface for users to invoke in their own code, while OLLA defines an interface
that users must implement, and that is invoked by the OpenMM implementation
layer.

At the lowest level is hardware specific code that actually performs
computations.  This code may be written in any language and use any technologies
that are appropriate.  For example, code for GPUs will be written in stream
processing languages such as OpenCL or CUDA, code written to run on clusters
will use MPI or other distributed computing tools, code written for multicore
processors will use threading tools such as Pthreads or OpenMP, etc.  OpenMM
sets no restrictions on how these computational kernels are written.  As long as
they are wrapped in the appropriate OLLA interfaces, OpenMM can use them.

.. _the-openmm-public-api:

The OpenMM Public API
*********************

The public API is based on a small number of classes:

**System**\ : A System specifies generic properties of the system to be
simulated: the number of particles it contains, the mass of each one, the size
of the periodic box, etc.  The interactions between the particles are specified
through a set of Force objects (see below) that are added to the System.  Force
field specific parameters, such as particle charges, are not direct properties
of the System.  They are properties of the Force objects contained within the
System.

**Force**\ : The Force objects added to a System define the behavior of the
particles.  Force is an abstract class; subclasses implement specific behaviors.
The Force class is actually slightly more general than its name suggests.  A
Force can, indeed, apply forces to particles, but it can also directly modify
particle positions and velocities in arbitrary ways.  Some thermostats and
barostats, for example, can be implemented as Force classes.  Examples of Force
subclasses include HarmonicBondForce, NonbondedForce, and MonteCarloBarostat.

**Context**\ : This stores all of the state information for a simulation:
particle positions and velocities, as well as arbitrary parameters defined by
the Forces in the System.  It is possible to create multiple Contexts for a
single System, and thus have multiple simulations of that System in progress at
the same time.

**Integrator**\ : This implements an algorithm for advancing the simulation
through time.  It is an abstract class; subclasses implement specific
algorithms.  Examples of Integrator subclasses include LangevinIntegrator,
VerletIntegrator, and BrownianIntegrator.

**State**\ : A State stores a snapshot of the simulation at a particular point
in time.  It is created by calling a method on a Context.  As discussed earlier,
this is a potentially expensive operation.  This is the only way to query the
values of state variables, such as particle positions and velocities; Context
does not provide methods for accessing them directly.

Here is an example of what the source code to create a System and run a
simulation might look like:

.. code-block:: c

    System system;
    for (int i = 0; i < numParticles; ++i)
        system.addParticle(particle[i].mass);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    system.addForce(bonds);
    for (int i = 0; i < numBonds; ++i)
        bonds->addBond(bond[i].particle1, bond[i].particle2,
            bond[i].length, bond[i].k);
    HarmonicAngleForce* angles = new HarmonicAngleForce();
    system.addForce(angles);
    for (int i = 0; i < numAngles; ++i)
        angles->addAngle(angle[i].particle1, angle[i].particle2,
            angle[i].particle3, angle[i].angle, angle[i].k);
    // ...create and initialize other force field terms in the same way
    LangevinMiddleIntegrator integrator(temperature, friction, stepSize);
    Context context(system, integrator);
    context.setPositions(initialPositions);
    context.setVelocities(initialVelocities);
    integrator.step(10000);

We create a System, add various Forces to it, and set parameters on both the
System and the Forces.  We then create a LangevinMiddleIntegrator, initialize a
Context in which to run a simulation, and instruct the Integrator to advance the
simulation for 10,000 time steps.

The OpenMM Low Level API
************************

The OpenMM Low Level API (OLLA) defines a set of interfaces that users must
implement in their own code if they want to extend OpenMM, such as to create a
new Force subclass or support a new hardware platform.  It is based on the
concept of “kernels” that define particular computations to be performed.

More specifically, there is an abstract class called **KernelImpl**\ .
Instances of this class (or rather, of its subclasses) are created by
**KernelFactory** objects.  These classes provide the concrete implementations
of kernels for a particular platform.  For example, to perform calculations on a
GPU, one would create one or more KernelImpl subclasses that implemented the
computations with GPU kernels, and one or more KernelFactory subclasses to
instantiate the KernelImpl objects.

All of these objects are encapsulated in a single object that extends
**Platform**\ . KernelFactory objects are registered with the Platform to be
used for creating specific named kernels.  The choice of what implementation to
use (a GPU implementation, a multithreaded CPU implementation, an MPI-based
distributed implementation, etc.) consists entirely of choosing what Platform to
use.

As discussed so far, the low level API is not in any way specific to molecular
simulation; it is a fairly generic computational API.  In addition to defining
the generic classes, OpenMM also defines abstract subclasses of KernelImpl
corresponding to specific calculations.  For example, there is a class called
CalcHarmonicBondForceKernel to implement HarmonicBondForce and a class called
IntegrateLangevinMiddleStepKernel to implement LangevinMiddleIntegrator.  It is
these classes for which each Platform must provide a concrete subclass.

This architecture is designed to allow easy extensibility.  To support a new
hardware platform, for example, you create concrete subclasses of all the
abstract kernel classes, then create appropriate factories and a Platform
subclass to bind everything together.  Any program that uses OpenMM can then use
your implementation simply by specifying your Platform subclass as the platform
to use.

Alternatively, you might want to create a new Force subclass to implement a new
type of interaction.  To do this, define an abstract KernelImpl subclass
corresponding to the new force, then write the Force class to use it.  Any
Platform can support the new Force by providing a concrete implementation of
your KernelImpl subclass.  Furthermore, you can easily provide that
implementation yourself, even for existing Platforms created by other people.
Simply create a new KernelFactory subclass for your kernel and register it with
the Platform object.  The goal is to have a completely modular system.  Each
module, which might be distributed as an independent library, can either add new
features to existing platforms or support existing features on new platforms.

In fact, there is nothing “special” about the kernel classes defined by OpenMM.
They are simply KernelImpl subclasses that happen to be used by Forces and
Integrators that happen to be bundled with OpenMM.  They are treated exactly
like any other KernelImpl, including the ones you define yourself.

It is important to understand that OLLA defines an interface, not an
implementation.  It would be easy to assume a one-to-one correspondence between
KernelImpl objects and the pieces of code that actually perform calculations,
but that need not be the case.  For a GPU implementation, for example, a single
KernelImpl might invoke several GPU kernels.  Alternatively, a single GPU kernel
might perform the calculations of several KernelImpl subclasses.

.. _platforms:

Platforms
*********

This release of OpenMM contains the following Platform subclasses:

**ReferencePlatform**\ : This is designed to serve as reference code for
writing other platforms.  It is written with simplicity and clarity in mind, not
performance.

**CpuPlatform**\ : This platform provides high performance when running on
conventional CPUs.

**CudaPlatform**\ : This platform is implemented using the CUDA language, and
performs calculations on Nvidia GPUs.

**HipPlatform**\ : This platform is implemented using the HIP language, and
performs calculations on ROCm-compatible AMD GPUs.

**OpenCLPlatform**\ : This platform is implemented using the OpenCL language,
and performs calculations on a variety of types of GPUs and CPUs.

The choice of which platform to use for a simulation depends on various factors:

#. The Reference platform is much slower than the others, and therefore is
   rarely used for production simulations.
#. The CPU platform is usually the fastest choice when a fast GPU is not
   available.  However, it requires the CPU to support SSE 4.1.  That includes most
   CPUs made in the last several years, but this platform may not be available on
   some older computers.  Also, for simulations that use certain features
   (primarily the various “custom” force classes), it may be faster to use the
   OpenCL platform running on the CPU.
#. The CUDA platform can be used with NVIDIA GPUs.  For using an AMD GPU,
   use the HIP platform (or the OpenCL platform which is usually slower), for
   using an Intel GPU, use the OpenCL platform.
#. The AMOEBA force field works with all platforms, but the performance
   of the Reference and CPU platforms is usually too slow to be useful.

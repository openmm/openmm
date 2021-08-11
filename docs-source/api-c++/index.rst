==============
OpenMM C++ API
==============

The C++ API provides information about the classes and methods available in OpenMM for C++ developers. OpenMM uses an object-oriented API that makes all its functionality available through a small number of classes.

Core classes
============

.. toctree::
    :maxdepth: 1
    :hidden:

    generated/System
    generated/Context
    generated/State
    generated/Platform

:cpp:class:`OpenMM::System`
---------------------------

A ``System`` specifies generic properties of the molecular system to be
simulated: the number of particles it contains, the mass of each one, the size
of the periodic box, and so on.  The interactions between the particles are
specified through a set of :ref:`Force <forces>` objects that are added to the
``System``.  Force field specific parameters, such as particle charges, are
stored in these ``Force`` objects, not as direct properties of the ``System``.

:cpp:class:`OpenMM::Context`
----------------------------

A ``Context`` stores all of the state information for a simulation: particle
positions and velocities, as well as arbitrary parameters defined by the
``Forces`` in the System.  It is possible to create multiple ``Contexts`` for a
single ``System``, and thus have multiple simulations of that ``System`` in
progress at the same time. ``Context`` does not provide methods for accessing
state variables directly; they must be read via a ``State`` object.


:cpp:class:`OpenMM::State`
--------------------------

A ``State`` object must be constructed before data can be read from a
simulation. State variables are not accessible directly via a ``Context`` in
order to make explicit the precise time that a variable reflects. A ``State``
is created by calling a method on a ``Context`` and stores only the information
requested at invocation.

:cpp:class:`OpenMM::Platform`
-----------------------------

A ``Platform`` is a single implementation of OpenMM at a low level. This allows
the same high level API documented here to be used on all sorts of compute
hardware, from GPUs to supercomputers. A ``Platform`` implements some set of
kernels, which define which operations it supports. Writing a new ``Platform``
allows OpenMM to be ported to new hardware or to be implemented in a new way
without rewriting the entire application.

Forces
======

``Force`` objects define the behavior of the particles in a ``System``. The
``Force`` class is actually slightly more general than its name suggests.  A
``Force`` can, indeed, apply forces to particles, but it can also directly
modify particle positions and velocities in arbitrary ways.  Some thermostats
and barostats, for example, can be implemented as ``Force`` classes.  Examples
of Force subclasses include :cpp:class:`HarmonicBondForce
<OpenMM::HarmonicBondForce>`, :cpp:class:`NonbondedForce
<OpenMM::NonbondedForce>`, and :cpp:class:`MonteCarloBarostat
<OpenMM::MonteCarloBarostat>`.

.. toctree::
    :maxdepth: 2

    forces

Integrators
===========

An ``Integrator`` implements an algorithm for advancing the simulation through
time.  They provide a ``Context`` a means of stepping the simulation forward,
and must be coupled to a ``Context`` to function. Examples of Integrator
subclasses include :cpp:class:`LangevinIntegrator <OpenMM::LangevinIntegrator>`,
:cpp:class:`VerletIntegrator <OpenMM::VerletIntegrator>`, and :cpp:class:`BrownianIntegrator <OpenMM::BrownianIntegrator>`.

.. toctree::
    :maxdepth: 2

    integrators

Extras
======

OpenMM's public API includes a few more classes that support the above.

.. toctree::
    :maxdepth: 2

    extras

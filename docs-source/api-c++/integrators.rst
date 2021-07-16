===========
Integrators
===========

The ``Integrator`` abstract class
=================================

.. toctree::
    :maxdepth: 2

    generated/OpenMM::Integrator

An ``Integrator`` implements an algorithm for advancing the simulation
through time.  It is an abstract class; subclasses implement specific
algorithms.


Meta ``Integrator`` classes
===========================

These classes facilitate customisation of the integrator. ``CustomIntegrator``
allows a wide variety of integration algorithms to be implemented efficiently
without writing any low-level code. The integrator is built up as a series of
steps, each defined as an algebraic expression. ``CompoundIntegrator`` allows
different integrators to be combined by making it possible to switch the active
integrator in the middle of a simulation.

.. toctree::
    :maxdepth: 2

    generated/OpenMM::CustomIntegrator
    generated/OpenMM::CompoundIntegrator

General purpose integrators
===========================

.. toctree::
    :maxdepth: 2

    generated/OpenMM::BrownianIntegrator
    generated/OpenMM::LangevinIntegrator
    generated/OpenMM::LangevinMiddleIntegrator
    generated/OpenMM::NoseHooverIntegrator
    generated/OpenMM::VariableLangevinIntegrator
    generated/OpenMM::VariableVerletIntegrator
    generated/OpenMM::VerletIntegrator

Drude integrators
=================

These integrators permit modelling polarization with a Drude particle.

.. toctree::
    :maxdepth: 2

    generated/OpenMM::DrudeIntegrator
    generated/OpenMM::DrudeLangevinIntegrator
    generated/OpenMM::DrudeNoseHooverIntegrator
    generated/OpenMM::DrudeSCFIntegrator

Ring Polymer Molecular Dynamics integrators
===========================================

.. toctree::
    :maxdepth: 2

    generated/OpenMM::RPMDIntegrator



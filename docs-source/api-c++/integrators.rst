===========
Integrators
===========

The ``Integrator`` abstract class
=================================

An ``Integrator`` implements an algorithm for advancing the simulation through
time.  ``Integrator`` is an abstract class; subclasses implement specific
algorithms.

.. toctree::
    :maxdepth: 2

    generated/Integrator


General purpose integrators
===========================

These are integrators appropriate for traditional MD and BD simulations.

.. toctree::
    :maxdepth: 2

    generated/BrownianIntegrator
    generated/LangevinIntegrator
    generated/LangevinMiddleIntegrator
    generated/NoseHooverIntegrator
    generated/VariableLangevinIntegrator
    generated/VariableVerletIntegrator
    generated/VerletIntegrator


Drude integrators
=================

These integrators permit modelling polarization with a Drude particle.

.. toctree::
    :maxdepth: 2

    generated/DrudeIntegrator
    generated/DrudeLangevinIntegrator
    generated/DrudeNoseHooverIntegrator
    generated/DrudeSCFIntegrator


Ring Polymer Molecular Dynamics integrators
===========================================

The RPMD integrator implements Ring Polymer MD.

.. toctree::
    :maxdepth: 2

    generated/RPMDIntegrator


Dissipative Particle Dynamics integrators
=========================================

The DPD integrator implements Dissipative Particle Dynamics.

.. toctree::
    :maxdepth: 2

    generated/DPDIntegrator


Custom integrators
==================

These classes facilitate customisation of the integrator. ``CustomIntegrator``
allows a wide variety of integration algorithms to be implemented efficiently
without writing any low-level code. The integrator is built up as a series of
steps, each defined as an algebraic expression. ``CompoundIntegrator`` allows
different integrators to be combined by making it possible to switch the active
integrator in the middle of a simulation.

.. toctree::
    :maxdepth: 2

    generated/CustomIntegrator
    generated/CompoundIntegrator

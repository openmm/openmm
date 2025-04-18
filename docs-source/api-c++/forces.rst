.. _forces:

======
Forces
======

The ``Force`` abstract class
============================

The ``Force`` objects added to a ``System`` define the behavior of the
particles. ``Force`` is an abstract class; subclasses implement specific behaviors. Classes that extend ``Force`` may implement actual physical forces, or any number of processes that either actually apply forces to particles or directly modify their positions or momenta.

.. toctree::
    :maxdepth: 2

    generated/Force


Common bonded and non-bonded forces
===================================

These classes implement forces that are widely used in biomolecular simulation.

.. toctree::
    :maxdepth: 2

    generated/CMAPTorsionForce
    generated/DrudeForce
    generated/GBSAOBCForce
    generated/GayBerneForce
    generated/HarmonicAngleForce
    generated/HarmonicBondForce
    generated/NonbondedForce
    generated/PeriodicTorsionForce
    generated/RBTorsionForce


AMOEBA forces
=============

These forces are used to implement the polarizable AMOEBA force fields.

.. toctree::
    :maxdepth: 2

    generated/AmoebaGeneralizedKirkwoodForce
    generated/AmoebaMultipoleForce
    generated/AmoebaTorsionTorsionForce
    generated/AmoebaVdwForce
    generated/AmoebaWcaDispersionForce
    generated/HippoNonbondedForce


Pseudo-forces
=============

These inherit from ``Force``, but do not describe physical forces. They are used
to implement thermostats or barostats, or otherwise modify the simulation from
step to step. They are conceptually closer to modifications to the integrator,
but providing them as a ``Force`` simplifies implementation and allows them to
be combined in arbitrary ways.

.. toctree::
    :maxdepth: 2

    generated/AndersenThermostat
    generated/ATMForce
    generated/CMMotionRemover
    generated/MonteCarloAnisotropicBarostat
    generated/MonteCarloBarostat
    generated/MonteCarloFlexibleBarostat
    generated/MonteCarloMembraneBarostat
    generated/RMSDForce
    generated/RPMDMonteCarloBarostat


.. _custom-forces:

Customizing ``Force``
=====================

OpenMM provides a number of classes that make it easier to implement custom
forces for common scenarios. These classes implement constructors that take an
algebraic expression as a string. The class is instantiated (not extended) to
provide a ``Force`` object that efficiently implements the provided
expression.

.. toctree::
    :maxdepth: 2

    generated/CustomAngleForce
    generated/CustomBondForce
    generated/CustomCVForce
    generated/CustomCentroidBondForce
    generated/CustomCompoundBondForce
    generated/CustomExternalForce
    generated/CustomGBForce
    generated/CustomHbondForce
    generated/CustomManyParticleForce
    generated/CustomNonbondedForce
    generated/CustomTorsionForce
    generated/CustomVolumeForce

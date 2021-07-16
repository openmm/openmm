======
Forces
======

The ``Force`` abstract class
============================

.. toctree::
    :maxdepth: 2

    generated/OpenMM::Force

The ``Force`` objects added to a ``System`` define the behavior of the
particles. ``Force`` is an abstract class; subclasses implement specific behaviors. Classes that extend ``Force`` may implement actual physical forces, or any number of processes that either actually apply forces to particles or directly modify their positions or momenta.

.. _custom-forces:

The ``Custom`` force classes
============================

.. toctree::
    :maxdepth: 2
    :glob:

    generated/OpenMM::CustomAngleForce
    generated/OpenMM::CustomBondForce
    generated/OpenMM::CustomCVForce
    generated/OpenMM::CustomCentroidBondForce
    generated/OpenMM::CustomCompoundBondForce
    generated/OpenMM::CustomExternalForce
    generated/OpenMM::CustomGBForce
    generated/OpenMM::CustomHbondForce
    generated/OpenMM::CustomManyParticleForce
    generated/OpenMM::CustomNonbondedForce
    generated/OpenMM::CustomTorsionForce

OpenMM provides a number of classes that make it easier to implement custom
forces for common scenarios. These classes implement constructors that take an
algebraic expression as a string. The class is instantiated (not extended) to
provide a ``Force`` object that efficiently implements the provided
expression.

Common bonded and non-bonded forces
===================================

.. toctree::
    :maxdepth: 2

    generated/OpenMM::CMAPTorsionForce
    generated/OpenMM::DrudeForce
    generated/OpenMM::GBSAOBCForce
    generated/OpenMM::GayBerneForce
    generated/OpenMM::HarmonicAngleForce
    generated/OpenMM::HarmonicBondForce
    generated/OpenMM::HippoNonbondedForce
    generated/OpenMM::NonbondedForce
    generated/OpenMM::PeriodicTorsionForce
    generated/OpenMM::RBTorsionForce

AMOEBA forces
=============

These forces are used to implement the polarizable AMOEBA force fields.

.. toctree::
    :maxdepth: 2

    generated/OpenMM::AmoebaGeneralizedKirkwoodForce
    generated/OpenMM::AmoebaMultipoleForce
    generated/OpenMM::AmoebaTorsionTorsionForce
    generated/OpenMM::AmoebaVdwForce
    generated/OpenMM::AmoebaWcaDispersionForce

Pseudo-forces
=============

These inherit from ``Force``, but do not describe physical forces. They are used
to implement thermostats or barostats, or otherwise modify the simulation from
step to step. They are conceptually closer to modifications to the integrator,
but providing them as a ``Force`` simplifies implementation and allows them to
be combined in arbitrary ways.

.. toctree::
    :maxdepth: 2

    generated/OpenMM::AndersenThermostat
    generated/OpenMM::CMMotionRemover
    generated/OpenMM::MonteCarloAnisotropicBarostat
    generated/OpenMM::MonteCarloBarostat
    generated/OpenMM::MonteCarloMembraneBarostat
    generated/OpenMM::RMSDForce
    generated/OpenMM::RPMDMonteCarloBarostat

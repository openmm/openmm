.. _the-theory-behind-openmm-introduction:

Introduction
############

Overview
********

This guide describes the mathematical theory behind OpenMM.  For each
computational class, it describes what computations the class performs and how
it should be used.  This serves two purposes.  If you are using OpenMM within an
application, this guide teaches you how to use it correctly.  If you are
implementing the OpenMM API for a new Platform, it teaches you how to correctly
implement the required kernels.

On the other hand, many details are intentionally left unspecified.  Any
behavior that is not specified either in this guide or in the API documentation
is left up to the Platform, and may be implemented in different ways by
different Platforms.  For example, an Integrator is required to produce a
trajectory that satisfies constraints to within the user-specified tolerance,
but the algorithm used to enforce those constraints is left up to the Platform.
Similarly, this guide provides the functional form of each Force, but does not
specify what level of numerical precision it must be calculated to.

This is an essential feature of the design of OpenMM, because it allows the API
to be implemented efficiently on a wide variety of hardware and software
platforms, using whatever methods are most appropriate for each platform.  On
the other hand, it means that a single program may produce meaningfully
different results depending on which Platform it uses.  For example, different
constraint algorithms may have different regions of convergence, and thus a time
step that is stable on one platform may be unstable on a different one.  It is
essential that you validate your simulation methodology on each Platform you
intend to use, and do not assume that good results on one Platform will
guarantee good results on another Platform when using identical parameters.


.. _units:

Units
*****

There are several different sets of units widely used in molecular simulations.
For example, energies may be measured in kcal/mol or kJ/mol, distances may be in
Angstroms or nm, and angles may be in degrees or radians.  OpenMM uses the
following units everywhere.

===========  =================
Quantity     Units
===========  =================
distance     nm
time         ps
mass         atomic mass units
charge       proton charge
temperature  Kelvin
angle        radians
energy       kJ/mol
===========  =================

These units have the important feature that they form an internally consistent
set.  For example, a force always has the same units (kJ/mol/nm) whether it is
calculated as the gradient of an energy or as the product of a mass and an
acceleration.  This is not true in some other widely used unit systems, such as
those that express energy in kcal/mol.

The header file Units.h contains predefined constants for converting between the
OpenMM units and some other common units.  For example, if your application
expresses distances in Angstroms, you should multiply them by
OpenMM::NmPerAngstrom before passing them to OpenMM, and positions calculated by
OpenMM should be multiplied by OpenMM::AngstromsPerNm before passing them back
to your application.

.. _physical-constants:

Physical Constants
******************

OpenMM uses the CODATA 2018 values for all physical constants.  Here are the
specific values it uses for the constants that frequently come up in molecular
simulations.

========================================  =================================
Quantity                                  Value
========================================  =================================
Elementary Charge (:math:`e`)             1.602176634·10\ :sup:`-19` C
Boltzmann's Constant (:math:`k_B`)        1.380649·10\ :sup:`-23` J/K
Avogadro's Number (:math:`N_A`)           6.02214076·10\ :sup:`23`
Vacuum Permittivity (:math:`\epsilon_0`)  8.8541878128·10\ :sup:`-12` F/m
========================================  =================================


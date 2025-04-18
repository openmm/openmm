.. _custom-forces:

Custom Forces
#############

In addition to the standard forces described in the previous chapter, OpenMM
provides a number of “custom” force classes.   These classes provide detailed
control over the mathematical form of the force by allowing the user to specify
one or more arbitrary algebraic expressions.  The details of how to write these
custom expressions are described in section :numref:`writing-custom-expressions`\ .

CustomBondForce
***************

CustomBondForce is similar to HarmonicBondForce in that it represents an
interaction between certain pairs of particles as a function of the distance
between them, but it allows the precise form of the interaction to be specified
by the user.  That is, the interaction energy of each bond is given by


.. math::
   E=f\left(r\right)


where *f*\ (\ *r*\ ) is a user defined mathematical expression.

In addition to depending on the inter-particle distance *r*\ , the energy may
also depend on an arbitrary set of user defined parameters.  Parameters may be
specified in two ways:

* Global parameters have a single, fixed value.
* Per-bond parameters are defined by specifying a value for each bond.


CustomAngleForce
****************

CustomAngleForce is similar to HarmonicAngleForce in that it represents an
interaction between sets of three particles as a function of the angle between
them, but it allows the precise form of the interaction to be specified by the
user.  That is, the interaction energy of each angle is given by


.. math::
   E=f\left(\theta\right)


where :math:`f(\theta)` is a user defined mathematical expression.

In addition to depending on the angle :math:`\theta`\ , the energy may also depend on an
arbitrary set of user defined parameters.  Parameters may be specified in two
ways:

* Global parameters have a single, fixed value.
* Per-angle parameters are defined by specifying a value for each angle.


CustomTorsionForce
******************

CustomTorsionForce is similar to PeriodicTorsionForce in that it represents an
interaction between sets of four particles as a function of the dihedral angle
between them, but it allows the precise form of the interaction to be specified
by the user.  That is, the interaction energy of each angle is given by


.. math::
   E=f(\theta)


where :math:`f(\theta)` is a user defined mathematical expression.  The angle
:math:`\theta` is guaranteed to be in the range :math:`[-\pi, +\pi]`\ .  Like PeriodicTorsionForce, it
is defined to be zero when the first and last particles are on the same side of
the bond formed by the middle two particles (the *cis* configuration).

In addition to depending on the angle :math:`\theta`\ , the energy may also depend on an
arbitrary set of user defined parameters.  Parameters may be specified in two
ways:

* Global parameters have a single, fixed value.
* Per-torsion parameters are defined by specifying a value for each torsion.


.. _customnonbondedforce:

CustomNonbondedForce
********************

CustomNonbondedForce is similar to NonbondedForce in that it represents a
pairwise interaction between all particles in the System, but it allows the
precise form of the interaction to be specified by the user.  That is, the
interaction energy between each pair of particles is given by


.. math::
   E=f(r)


where *f*\ (\ *r*\ ) is a user defined mathematical expression.

In addition to depending on the inter-particle distance *r*\ , the energy may
also depend on an arbitrary set of user defined parameters.  Parameters may be
specified in two ways:

* Global parameters have a single, fixed value.
* Per-particle parameters are defined by specifying a value for each particle.


A CustomNonbondedForce can optionally be restricted to only a subset of particle
pairs in the System.  This is done by defining “interaction groups”.  See the
API documentation for details.

When using a cutoff, a switching function can optionally be applied to make the
energy go smoothly to 0 at the cutoff distance.  When
:math:`r_\mathit{switch} < r < r_\mathit{cutoff}`\ , the energy is multiplied by



.. math::
   S=1-{6x}^{5}+15{x}^{4}-10{x}^{3}


where :math:`x=(r-r_\mathit{switch})/(r_\mathit{cutoff}-r_\mathit{switch})`\ .
This function decreases smoothly from 1 at :math:`r=r_\mathit{switch}`
to 0 at :math:`r=r_\mathit{cutoff}`\ , and has continuous first and
second derivatives at both ends.

When using periodic boundary conditions, CustomNonbondedForce can optionally add
a term (known as a *long range truncation correction*\ ) to the energy that
approximately represents the contribution from all interactions beyond the
cutoff distance:\ :cite:`Shirts2007`


.. math::
   {E}_{cor}=\frac{2\pi N^2}{V}\left\langle\underset{{r}_\mathit{cutoff}}{\overset{\infty}{\int}}E(r)r^{2}dr\right\rangle


where *N* is the number of particles in the system, *V* is the volume of
the periodic box, and :math:`\langle \text{...} \rangle` represents an average over all pairs of particles in
the system.  When a switching function is in use, there is an additional
contribution to the correction given by


.. math::
   E_{cor}^\prime=\frac{2\pi N^2}{V}\left\langle\underset{{r}_\mathit{switch}}{\overset{{r}_\mathit{cutoff}}{\int }}E(r)(1-S(r))r^{2}dr\right\rangle


The long range dispersion correction is primarily useful when running
simulations at constant pressure, since it produces a more accurate variation in
system energy with respect to volume.

CustomExternalForce
*******************

CustomExternalForce represents a force that is applied independently to each
particle as a function of its position.   That is, the energy of each particle
is given by


.. math::
   E=f(x,y,z)


where *f*\ (\ *x*\ , *y*\ , *z*\ ) is a user defined mathematical
expression.

In addition to depending on the particle’s (\ *x*\ , *y*\ , *z*\ )
coordinates, the energy may also depend on an arbitrary set of user defined
parameters.  Parameters may be specified in two ways:

* Global parameters have a single, fixed value.
* Per-particle parameters are defined by specifying a value for each particle.


CustomCompoundBondForce
***********************

CustomCompoundBondForce supports a wide variety of bonded interactions.  It
defines a “bond” as a single energy term that depends on the positions of a
fixed set of particles.  The number of particles involved in a bond, and how the
energy depends on their positions, is configurable.  It may depend on the
positions of individual particles, the distances between pairs of particles, the
angles formed by sets of three particles, and the dihedral angles formed by sets
of four particles.  That is, the interaction energy of each bond is given by


.. math::
   E=f(\{x_i\},\{r_i\},\{\theta_i\},\{\phi_i\})


where *f*\ (\ *...*\ ) is a user defined mathematical expression.  It may
depend on an arbitrary set of positions {\ :math:`x_i`\ }, distances {\ :math:`r_i`\ },
angles {\ :math:`\theta_i`\ }, and dihedral angles {\ :math:`\phi_i`\ }
guaranteed to be in the range :math:`[-\pi, +\pi]`\ .

Each distance, angle, or dihedral is defined by specifying a sequence of
particles chosen from among the particles that make up the bond.  A distance
variable is defined by two particles, and equals the distance between them.  An
angle variable is defined by three particles, and equals the angle between them.
A dihedral variable is defined by four particles, and equals the angle between
the first and last particles about the axis formed by the middle two particles.
It is equal to zero when the first and last particles are on the same side of
the axis.

In addition to depending on positions, distances, angles, and dihedrals, the
energy may also depend on an arbitrary set of user defined parameters.
Parameters may be specified in two ways:

* Global parameters have a single, fixed value.
* Per-bond parameters are defined by specifying a value for each bond.


CustomCentroidBondForce
***********************

CustomCentroidBondForce is very similar to CustomCompoundBondForce, but instead
of creating bonds between individual particles, the bonds are between the
centers of groups of particles.  This is useful for purposes such as restraining
the distance between two molecules or pinning the center of mass of a single
molecule.

The first step in computing this force is to calculate the center position of
each defined group of particles.  This is calculated as a weighted average of
the positions of all the particles in the group, with the weights being user
defined.  The computation then proceeds exactly as with CustomCompoundBondForce,
but the energy of each "bond" is now calculated based on the centers of a set
of groups, rather than on the positions of individual particles.

This class supports all the same function types and features as
CustomCompoundBondForce.  In fact, any interaction that could be implemented
with CustomCompoundBondForce can also be implemented with this class, simply by
defining each group to contain only a single atom.


CustomManyParticleForce
***********************

CustomManyParticleForce is similar to CustomNonbondedForce in that it represents
a custom nonbonded interaction between particles, but it allows the interaction
to depend on more than two particles.  This allows it to represent a wide range
of non-pairwise interactions.  It is defined by specifying the number of
particles :math:`N` involved in the interaction and how the energy depends on
their positions.  More specifically, it takes a user specified energy function

.. math::
   E=f(\{x_i\},\{r_i\},\{\theta_i\},\{\phi_i\})

that may depend on an arbitrary set of positions {\ :math:`x_i`\ }, distances
{\ :math:`r_i`\ }, angles {\ :math:`\theta_i`\ }, and dihedral angles
{\ :math:`\phi_i`\ } from a particular set of :math:`N` particles.

Each distance, angle, or dihedral is defined by specifying a sequence of
particles chosen from among the particles in the set.  A distance
variable is defined by two particles, and equals the distance between them.  An
angle variable is defined by three particles, and equals the angle between them.
A dihedral variable is defined by four particles, and equals the angle between
the first and last particles about the axis formed by the middle two particles.
It is equal to zero when the first and last particles are on the same side of
the axis.

In addition to depending on positions, distances, angles, and dihedrals, the
energy may also depend on an arbitrary set of user defined parameters.
Parameters may be specified in two ways:

* Global parameters have a single, fixed value.
* Per-particle parameters are defined by specifying a value for each particle.

The energy function is evaluated one or more times for every unique set of
:math:`N` particles in the system.  The exact number of times depends on the
*permutation mode*\ .  A set of :math:`N` particles has :math:`N!` possible
permutations.  In :code:`SinglePermutation` mode, the function is evaluated
for a single arbitrarily chosen one of those permutations.  In
:code:`UniqueCentralParticle` mode, the function is evaluated for :math:`N` of
those permutations, once with each particle as the "central particle".

The number of times the energy function is evaluated can be further restricted
by specifying *type filters*\ .  Each particle may have a "type" assigned to it,
and then each of the :math:`N` particles involved in an interaction may be
restricted to only a specified set of types.  This provides a great deal of
flexibility in controlling which particles interact with each other.


CustomGBForce
*************

CustomGBForce implements complex, multiple stage nonbonded interactions between
particles.  It is designed primarily for implementing Generalized Born implicit
solvation models, although it is not strictly limited to that purpose.

The interaction is specified as a series of computations, each defined by an
arbitrary algebraic expression.  These computations consist of some number of
per-particle *computed values*\ , followed by one or more *energy terms*\ .
A computed value is a scalar value that is computed for each particle in the
system.  It may depend on an arbitrary set of global and per-particle
parameters, and well as on other computed values that have been calculated
before it.  Once all computed values have been calculated, the energy terms and
their derivatives are evaluated to determine the system energy and particle
forces.  The energy terms may depend on global parameters, per-particle
parameters, and per-particle computed values.

Computed values can be calculated in two different ways:

* *Single particle* values are calculated by evaluating a user defined
  expression for each particle:

.. math::
  {value}_{i}=f\left(\text{.}\text{.}\text{.}\right)
..

  where *f*\ (...) may depend only on properties of particle *i* (its
  coordinates and parameters, as well as other computed values that have already
  been calculated).

* *Particle pair* values are calculated as a sum over pairs of particles:

.. math::
  {value}_{i}=\sum _{j\ne i}f\left(r,\text{...}\right)
..

  where the sum is over all other particles in the System, and *f*\ (\ *r*\ ,
  ...) is a function of the distance *r* between particles *i* and *j*\,
  as well as their parameters and computed values.

Energy terms may similarly be calculated per-particle or per-particle-pair.

* *Single particle* energy terms are calculated by evaluating a user
  defined expression for each particle:

.. math::
  E=f\left(\text{.}\text{.}\text{.}\right)
..

  where *f*\ (...) may depend only on properties of that particle (its
  coordinates, parameters, and computed values).

* *Particle pair* energy terms are calculated by evaluating a user defined
  expression once for every pair of particles in the System:

.. math::
  E=\sum _{i,j}f\left(r,\text{.}\text{.}\text{.}\right)
..

  where the sum is over all particle pairs *i* *< j*\ , and *f*\ (\ *r*\ ,
  ...) is a function of the distance *r* between particles *i* and *j*\,
  as well as their parameters and computed values.

Note that energy terms are assumed to be symmetric with respect to the two
interacting particles, and therefore are evaluated only once per pair.  In
contrast, expressions for computed values need not be symmetric and therefore
are calculated twice for each pair: once when calculating the value for the
first particle, and again when calculating the value for the second particle.

Be aware that, although this class is extremely general in the computations it
can define, particular Platforms may only support more restricted types of
computations.  In particular, all currently existing Platforms require that the
first computed value *must* be a particle pair computation, and all computed
values after the first *must* be single particle computations. This is
sufficient for most Generalized Born models, but might not permit some other
types of calculations to be implemented.

CustomHbondForce
****************

CustomHbondForce supports a wide variety of energy functions used to represent
hydrogen bonding.  It computes interactions between "donor" particle groups and
"acceptor" particle groups, where each group may include up to three particles.
Typically a donor group consists of a hydrogen atom and the atoms it is bonded
to, and an acceptor group consists of a negatively charged atom and the atoms it
is bonded to.  The interaction energy between each donor group and each acceptor
group is given by


.. math::
   E=f(\{r_i\},\{\theta_i\},\{\phi_i\})


where *f*\ (\ *...*\ ) is a user defined mathematical expression.  It may
depend on an arbitrary set of distances {\ :math:`r_i`\ }, angles {\ :math:`\theta_i`\ },
and dihedral angles {\ :math:`\phi_i`\ }.

Each distance, angle, or dihedral is defined by specifying a sequence of
particles chosen from the interacting donor and acceptor groups (up to six atoms
to choose from, since each group may contain up to three atoms).  A distance
variable is defined by two particles, and equals the distance between them.  An
angle variable is defined by three particles, and equals the angle between them.
A dihedral variable is defined by four particles, and equals the angle between
the first and last particles about the axis formed by the middle two particles.
It is equal to zero when the first and last particles are on the same side of
the axis.

In addition to depending on distances, angles, and dihedrals, the energy may
also depend on an arbitrary set of user defined parameters.  Parameters may be
specified in three ways:

* Global parameters have a single, fixed value.
* Per-donor parameters are defined by specifying a value for each donor group.
* Per-acceptor parameters are defined by specifying a value for each acceptor group.

.. _customcvforce:

CustomCVForce
*************

CustomCVForce computes an energy as a function of "collective variables".  A
collective variable may be any scalar valued function of the particle positions
and other parameters.  Each one is defined by a :code:`Force` object, so any
function that can be defined via any force class (either standard or custom) can
be used as a collective variable.  The energy is then computed as

.. math::
   E=f(...)

where *f*\ (...) is a user supplied mathematical expression of the collective
variables.  It also may depend on user defined global parameters.

CustomVolumeForce
*****************

CustomVolumeForce computes an energy that depends only on the box vectors
defining the periodic box (see Section :numref:`periodic-boundary-conditions`).

.. math::
   E = f(\mathbf{a}, \mathbf{b}, \mathbf{c})

Because the energy does not depend on particle positions, it does not apply any
forces to particles.  It is primarily useful for constant pressure simulations,
where the volume-dependent energy can influence the behavior of the barostat.
Energy terms of this sort are often used for pressure matching in coarse grained
force fields.

ATMForce
********

ATMForce implements the Alchemical Transfer Method for free energy calculations.\ :cite:`Azimi2022`
It contains one or more :code:`Force` objects whose energy is evaluated twice,
before and after displacing some particles to new positions.  The final energy
is determined by a user supplied mathematical function of the two energies.  See
the API documentation and the publication for more details.

.. _writing-custom-expressions:

Writing Custom Expressions
**************************

The custom forces described in this chapter involve user defined algebraic
expressions.  These expressions are specified as character strings, and may
involve a variety of standard operators and mathematical functions.

The following operators are supported: + (add), - (subtract), * (multiply), /
(divide), and ^ (power).  Parentheses “(“ and “)” may be used for grouping.

The following standard functions are supported: sqrt, exp, log, sin, cos, sec,
csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs,
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

Setting Parameters
******************

Most custom forces have two types of parameters you can define.  The simplest type
are global parameters, which represent a single number.  The value is stored in
the :class:`Context`, and can be changed at any time by calling :meth:`setParameter`
on it.  Global parameters are designed to be very inexpensive to change.  Even if
you set a new value for a global parameter on every time step, the overhead will
usually be quite small.  There can be exceptions to this rule, however.  For
example, if a :class:`CustomNonbondedForce` uses a long range correction, changing
a global parameter may require the correction coefficient to be recalculated,
which is expensive.

It is possible for multiple forces to depend on the same global parameter.  To do this,
simply have each force specify a parameter with the same name.  This can be useful
in certain cases.  For example, in an alchemical simulation, you might have a
parameter that interpolates between two endpoints corresponding to different molecules.
Changing the one parameter would simultaneously modify multiple bonded and nonbonded
forces.

The other type of parameter is ones that record many values, one for each element
of the force, such as per-particle or per-bond parameters.  These values are stored
directly in the force object itself, and hence are part of the system definition.
When a :class:`Context` is created, the values are copied over to it, and thereafter
the two are disconnected.  Modifying the force will have no effect on any
:class:`Context` that already exists.

Some forces do provide a way to modify these parameters via an :meth:`updateParametersInContext`
method.  These methods tend to be somewhat expensive, so it is best not to call
them too often.  On the other hand, they are still much less expensive than calling
:meth:`reinitialize` on the :class:`Context`, which is the other way of updating
the system definition for a running simulation.

Parameter Derivatives
*********************

Many custom forces have the ability to compute derivatives of the potential energy
with respect to global parameters.  To use this feature, first define a global
parameter that the energy depends on.  Then instruct the custom force to compute
the derivative with respect to that parameter by calling :meth:`addEnergyParameterDerivative()`
on it.  Whenever forces and energies are computed, the specified derivative will
then also be computed at the same time.  You can query it by calling :meth:`getState()`
on a :class:`Context`, just as you would query forces or energies.

An important application of this feature is to use it in combination with a
:class:`CustomIntegrator` (described in section :numref:`custom-integrator`\ ).  The
derivative can appear directly in expressions that define the integration
algorithm.  This can be used to implement algorithms such as lambda-dynamics,
where a global parameter is integrated as a dynamic variable.


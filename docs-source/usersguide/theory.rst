.. include:: header.rst

.. _the-theory-behind-openmm-introduction:

The Theory Behind OpenMM: Introduction
######################################


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



Standard Forces
###############

The following classes implement standard force field terms that are widely used
in molecular simulations.

HarmonicBondForce
*****************

Each harmonic bond is represented by an energy term of the form



.. math::
   E=\frac{1}{2}k{\left(x-{x}_{0}\right)}^{2}


where *x* is the distance between the two particles, *x*\ :sub:`0` is
the equilibrium distance, and *k* is the force constant.  This produces a
force of magnitude *k*\ (\ *x*\ -\ *x*\ :sub:`0`\ ).

Be aware that some force fields define their harmonic bond parameters in a
slightly different way: *E* = *k*\ ´(\ *x*\ -\ *x*\ :sub:`0`\ )\
:sup:`2`\ , leading to a force of magnitude 2\ *k*\ ´(\ *x*\ -\ *x*\ :sub:`0`\ ).
Comparing these two forms, you can see that *k* = 2\ *k*\ ´.  Be sure to
check which form a particular force field uses, and if necessary multiply the
force constant by 2.

HarmonicAngleForce
******************

Each harmonic angle is represented by an energy term of the form


.. math::
   E=\frac{1}{2}k{\left(\theta-\theta_0\right)}^{2}


where :math:`\theta` is the angle formed by the three particles, :math:`\theta_0` is
the equilibrium angle, and *k* is the force constant.

As with HarmonicBondForce, be aware that some force fields define their harmonic
angle parameters as *E* = *k*\ ´(\ :math:`\theta`\ -\ :math:`\theta`\ :sub:`0`\ )\ :sup:`2`\ .
Be sure to check which form a particular force field uses, and if necessary
multiply the force constant by 2.

PeriodicTorsionForce
********************

Each torsion is represented by an energy term of the form


.. math::
   E=k\left(1+\text{cos}\left(n\theta-\theta_0\right)\right)


where :math:`\theta` is the dihedral angle formed by the four particles, :math:`\theta_0`
is the phase offset, *n* is the periodicity, and *k* is
the force constant.

RBTorsionForce
**************

Each torsion is represented by an energy term of the form


.. math::
   E=\sum _{i=0}^{5}{C}_{i}{\left(\text{cos}\phi\right)}^{i}


where :math:`\phi` is the dihedral angle formed by the four particles and
*C*\ :sub:`0` through *C*\ :sub:`5` are constant coefficients.

For reason of convention, PeriodicTorsionForce and RBTorsonForce define the
torsion angle differently. :math:`\theta` is zero when the first and last particles are
on the *same* side of the bond formed by the middle two particles (the
*cis* configuration), whereas :math:`\phi` is zero when they are on *opposite*
sides (the *trans* configuration).  This means that :math:`\theta` = :math:`\phi` - :math:`\pi`.

CMAPTorsionForce
****************

Each torsion pair is represented by an energy term of the form


.. math::
   E=f\left(\theta_1,\theta_2\right)


where :math:`\theta_1` and :math:`\theta_2` are the two dihedral angles
coupled by the term, and *f*\ (\ *x*\ ,\ *y*\ ) is defined by a user-supplied
grid of tabulated values.  A natural cubic spline surface is fit through the
tabulated values, then evaluated to determine the energy for arbitrary (\ :math:`\theta_1`\ ,
:math:`\theta_2`\ ) pairs.

NonbondedForce
**************

.. _lennard-jones-interaction:

Lennard-Jones Interaction
=========================

The Lennard-Jones interaction between each pair of particles is represented by
an energy term of the form


.. math::
   E=4\epsilon\left({\left(\frac{\sigma}{r}\right)}^{12}-{\left(\frac{\sigma}{r}\right)}^{6}\right)


where *r* is the distance between the two particles, :math:`\sigma` is the distance
at which the energy equals zero, and :math:`\epsilon` sets the strength of the
interaction.  If the NonbondedMethod in use is anything other than NoCutoff and
\ *r* is greater than the cutoff distance, the energy and force are both set
to zero.  Because the interaction decreases very quickly with distance, the
cutoff usually has little effect on the accuracy of simulations.

Optionally you can use a switching function to make the energy go smoothly to 0
at the cutoff distance.  When :math:`r_\mathit{switch} < r < r_\mathit{cutoff}`\ , the energy is multiplied by

.. math::
   S=1-{6x}^{5}+15{x}^{4}-10{x}^{3}


where :math:`x = (r-r_\mathit{switch})/(r_\mathit{cutoff}-r_\mathit{switch})`. This function decreases smoothly from 1 at
:math:`r = r_\mathit{switch}` to 0 at :math:`r = r_\mathit{cutoff}`, and has continuous first and
second derivatives at both ends.

When an exception has been added for a pair of particles, :math:`\sigma` and :math:`\epsilon`
are the parameters specified by the exception.  Otherwise they are determined
from the parameters of the individual particles using the Lorentz-Berthelot
combining rule:

.. math::
   \sigma=\frac{\sigma_1+\sigma_2}{2}

.. math::
   \epsilon=\sqrt{\epsilon_1 \epsilon_2}

When using periodic boundary conditions, NonbondedForce can optionally add a
term (known as a *long range dispersion correction*\ ) to the energy that
approximately represents the contribution from all interactions beyond the
cutoff distance:\ :cite:`Shirts2007`\

.. math::
   {E}_{\text{cor}}=\frac{{8\pi N}^{2}}{V}\left(\frac{\langle \epsilon_{ij}\sigma_{ij}^{12}\rangle}{{9r_c}^9}-\frac{\langle \epsilon_{ij}\sigma_{ij}^{6}\rangle}{{3r_c}^3}\right)

where *N* is the number of particles in the system, *V* is the volume of
the periodic box, :math:`r_c` is the cutoff distance, :math:`\sigma_{ij}` and
:math:`\epsilon_{ij}` are the interaction parameters between particle *i* and
particle *j*\ , and :math:`\langle \text{...} \rangle` represents an average over all pairs of particles in
the system.  When a switching function is in use, there is also a contribution
to the correction that depends on the integral of *E*\ ·(1-\ *S*\ ) over the
switching interval.  The long range dispersion correction is primarily useful
when running simulations at constant pressure, since it produces a more accurate
variation in system energy with respect to volume.

The Lennard-Jones interaction is often parameterized in two other equivalent
ways.  One is


.. math::
   E=\epsilon\left({\left(\frac{{r}_{\mathit{min}}}{r}\right)}^{\text{12}}-2{\left(\frac{{r}_{\mathit{min}}}{r}\right)}^{6}\right)


where :math:`r_\mathit{min}` (sometimes known as :math:`d_\mathit{min}`; this is not a
radius) is the center-to-center distance at which the energy is minimum.  It is
related to :math:`\sigma` by


.. math::
   \sigma=\frac{{r}_{\mathit{min}}}{{2}^{1/6}}


In turn, :math:`r_\mathit{min}` is related to the van der Waals radius by :math:`r_\mathit{min} = 2r_\mathit{vdw}`\ .

Another common form is



.. math::
   E=\frac{A}{{r}^{\text{12}}}-\frac{B}{{r}^{6}}


The coefficients A and B are related to :math:`\sigma` and :math:`\epsilon` by



.. math::
   \sigma={\left(\frac{A}{B}\right)}^{1/6}



.. math::
   \epsilon=\frac{{B}^{2}}{4A}


Coulomb Interaction Without Cutoff
==================================

The form of the Coulomb interaction between each pair of particles depends on
the NonbondedMethod in use.  For NoCutoff, it is given by


.. math::
   E=\frac{1}{4{\pi}{\epsilon}_{0}}\frac{{q}_{1}{q}_{2}}{r}


where *q*\ :sub:`1` and *q*\ :sub:`2` are the charges of the two
particles, and *r* is the distance between them.

Coulomb Interaction With Cutoff
===============================

For CutoffNonPeriodic or CutoffPeriodic, it is modified using the reaction field
approximation.  This is derived by assuming everything beyond the cutoff
distance is a solvent with a uniform dielectric constant.\ :cite:`Tironi1995`


.. math::
   E=\frac{{q}_{1}{q}_{2}}{4\pi\epsilon_0}\left(\frac{1}{r}+{k}_{\mathit{rf}}{r}^{2}-{c}_{\mathit{rf}}\right)


.. math::
   {k}_{\mathit{rf}}=\left(\frac{1}{{r_\mathit{cutoff}}^3}\right)\left(\frac{{\epsilon}_{\mathit{solvent}}-1}{2{\epsilon}_{\mathit{solvent}}+1}\right)


.. math::
   {c}_{\mathit{rf}}=\left(\frac{1}{{r}_{\mathit{cutoff}}}\right)\left(\frac{3{\epsilon}_{\mathit{solvent}}}{2{\epsilon}_{\mathit{solvent}}+1}\right)


where :math:`r_\mathit{cutoff}` is the cutoff distance and :math:`\epsilon_\mathit{solvent}` is
the dielectric constant of the solvent.  In the limit :math:`\epsilon_\mathit{solvent}` >> 1,
this causes the force to go to zero at the cutoff.

Coulomb Interaction With Ewald Summation
========================================

For Ewald, the total Coulomb energy is the sum of three terms: the *direct
space sum*\ , the *reciprocal space sum*\ , and the *self-energy term*\ .\
:cite:`Toukmaji1996`


.. math::
   E=E_{\mathit{dir}}+{E}_{\mathit{rec}}+{E}_{\mathit{self}}


.. math::
   E_{\mathit{dir}}=\frac{1}{2}\sum _{i,j}\sum_\mathbf{n}{q}_{i}{q}_{j}\frac{\text{erfc}\left({\mathit{\alpha r}}_{ij,\mathbf{n}}\right)}{r_{ij,\mathbf{n}}}


.. math::
   E_{\mathit{rec}}=\frac{1}{2{\pi}V}\sum _{i,j}q_i q_j\sum _{\mathbf{k}{\neq}0}\frac{\text{exp}(-(\pi \mathbf{k}/\alpha)^2+2\pi i \mathbf{k} \cdot (\mathbf{r}_{i}-\mathbf{r}_{j}))}{\mathbf{m}^2}


.. math::
   E_{\mathit{self}}=-\frac{\alpha}{\sqrt{\pi}}\sum _{i}{q}_{{i}^{2}}


In the above expressions, the indices *i* and *j* run over all
particles, **n** = (n\ :sub:`1`\ , n\ :sub:`2`\ , n\ :sub:`3`\ ) runs over
all copies of the periodic cell, and **k** = (k\ :sub:`1`\ , k\ :sub:`2`\ ,
k\ :sub:`3`\ ) runs over all integer wave vectors from (-k\ :sub:`max`\ ,
-k\ :sub:`max`\ , -k\ :sub:`max`\ ) to (k\ :sub:`max`\ , k\ :sub:`max`\ ,
k\ :sub:`max`\ ) excluding (0, 0, 0).  :math:`\mathbf{r}_i` is the position of
particle i , while :math:`r_{ij}` is the distance between particles *i* and *j*\ .
*V* is the volume of the periodic cell, and :math:`\alpha` is an internal parameter.

In the direct space sum, all pairs that are further apart than the cutoff
distance are ignored.  Because the cutoff is required to be less than half the
width of the periodic cell, the number of terms in this sum is never greater
than the square of the number of particles.

The error made by applying the direct space cutoff depends on the magnitude of
:math:`\text{erfc}({\alpha}r_\mathit{cutoff})`\ .  Similarly, the error made in the reciprocal space
sum by ignoring wave numbers beyond k\ :sub:`max` depends on the magnitude
of :math:`\text{exp}(-({\pi}k_{max}/{\alpha})^2`\ ).  By changing :math:`\alpha`, one can decrease the
error in either term while increasing the error in the other one.

Instead of having the user specify :math:`\alpha` and -k\ :sub:`max`\ , NonbondedForce
instead asks the user to choose an error tolerance :math:`\delta`.  It then calculates :math:`\alpha` as


.. math::
   \alpha =\sqrt{-\text{log}\left(2{\delta}\right)}/{r}_{\mathit{cutoff}}


Finally, it estimates the error in the reciprocal space sum as


.. math::
   \mathit{error}=\frac{k_{\mathit{max}}\sqrt{d\alpha}}{20}\text{exp}(-(\pi k_\mathit{max}/d\alpha)^2)


where *d* is the width of the periodic box, and selects the smallest value
for k\ :sub:`max` which gives *error* < :math:`\delta`\ .  (If the box is not square,
k\ :sub:`max` will have a different value along each axis.)

This means that the accuracy of the calculation is determined by :math:`\delta`\ .
:math:`r_\mathit{cutoff}` does not affect the accuracy of the result, but does affect the speed
of the calculation by changing the relative costs of the direct space and
reciprocal space sums.  You therefore should test different cutoffs to find the
value that gives best performance; this will in general vary both with the size
of the system and with the Platform being used for the calculation.  When the
optimal cutoff is used for every simulation, the overall cost of evaluating the
nonbonded forces scales as O(N\ :sup:`3/2`\ ) in the number of particles.

Be aware that the error tolerance :math:`\delta` is not a rigorous upper bound on the errors.
The formulas given above are empirically found to produce average relative
errors in the forces that are less than or similar to :math:`\delta` across a variety of
systems and parameter values, but no guarantees are made.  It is important to
validate your own simulations, and identify parameter values that produce
acceptable accuracy for each system.

Coulomb Interaction With Particle Mesh Ewald
============================================

The Particle Mesh Ewald (PME) algorithm\ :cite:`Essmann1995` is similar to
Ewald summation, but instead of calculating the reciprocal space sum directly,
it first distributes the particle charges onto nodes of a rectangular mesh using
5th order B-splines.  By using a Fast Fourier Transform, the sum can then be
computed very quickly, giving performance that scales as O(N log N) in the
number of particles (assuming the volume of the periodic box is proportional to
the number of particles).

As with Ewald summation, the user specifies the direct space cutoff :math:`r_\mathit{cutoff}`
and error tolerance :math:`\delta`\ .  NonbondedForce then selects :math:`\alpha` as


.. math::
   \alpha =\sqrt{-\text{log}\left(2\delta\right)}/{r}_\mathit{cutoff}


and the number of nodes in the mesh along each dimension as


.. math::
   n_\mathit{mesh}=\frac{2\alpha d}{{3\delta}^{1/5}}


where *d* is the width of the periodic box along that dimension.  Alternatively,
the user may choose to explicitly set values for these parameters.  (Note that
some Platforms may choose to use a larger value of :math:`n_\mathit{mesh}` than that
given by this equation.  For example, some FFT implementations require the mesh
size to be a multiple of certain small prime numbers, so a Platform might round
it up to the nearest permitted value.  It is guaranteed that :math:`n_\mathit{mesh}`
will never be smaller than the value given above.)

The comments in the previous section regarding the interpretation of :math:`\delta` for Ewald
summation also apply to PME, but even more so.  The behavior of the error for
PME is more complicated than for simple Ewald summation, and while the above
formulas will usually produce an average relative error in the forces less than
or similar to :math:`\delta`\ , this is not a rigorous guarantee.  PME is also more sensitive
to numerical round-off error than Ewald summation.  For Platforms that do
calculations in single precision, making :math:`\delta` too small (typically below about
5·10\ :sup:`-5`\ ) can actually cause the error to increase.

Lennard-Jones Interaction With Particle Mesh Ewald
==================================================

The PME algorithm can also be used for Lennard-Jones interactions.  Usually this
is not necessary, since Lennard-Jones forces are short ranged, but there are
situations (such as membrane simulations) where neglecting interactions beyond
the cutoff can measurably affect results.

For computational efficiency, certain approximations are made\ :cite:`Wennberg2015`.
Interactions beyond the cutoff distance include only the attractive :math:`1/r^6`
term, not the repulsive :math:`1/r^{12}` term.  Since the latter is much smaller
than the former at long distances, this usually has negligible effect.  Also,
the interaction between particles farther apart than the cutoff distance is
computed using geometric combination rules:

.. math::
   \sigma=\sqrt{\sigma_1 \sigma_2}

The effect of this approximation is also quite small, and it is still far more
accurate than ignoring the interactions altogether (which is what would happen
with PME).

The formula used to compute the number of nodes along each dimension of the mesh
is slightly different from the one used for Coulomb interactions:

.. math::
   n_\mathit{mesh}=\frac{\alpha d}{{3\delta}^{1/5}}

As before, this is an empirical formula.  It will usually produce an average
relative error in the forces less than or similar to :math:`\delta`\ , but that
is not guaranteed.

.. _gbsaobcforce:

GBSAOBCForce
************


Generalized Born Term
=====================

GBSAOBCForce consists of two energy terms: a Generalized Born Approximation term
to represent the electrostatic interaction between the solute and solvent, and a
surface area term to represent the free energy cost of solvating a neutral
molecule.  The Generalized Born energy is given by\ :cite:`Onufriev2004`


.. math::
   E\text{=-}\frac{1}{2}\left(\frac{1}{\epsilon_{\mathit{solute}}}-\frac{1}{\epsilon_{\mathit{solvent}}}\right)\sum _{i,j}\frac{{q}_{i}{q}_{j}}{{f}_{\text{GB}}\left({d}_{ij},{R}_{i},{R}_{j}\right)}


where the indices *i* and *j* run over all particles, :math:`\epsilon_\mathit{solute}`
and :math:`\epsilon_\mathit{solvent}` are the dielectric constants of the solute and solvent
respectively, :math:`q_i` is the charge of particle *i*\ , and :math:`d_{ij}` is the distance
between particles *i* and *j*\ .  :math:`f_\text{GB}(d_{ij}, R_i, R_j)` is defined as


.. math::
   {f}_{\text{GB}}\left({d}_{ij},{R}_{i},{R}_{j}\right)={\left[{d}_{ij}^2+{R}_{i}{R}_{j}\text{exp}\left(\frac{-{d}_{ij}^2}{{4R}_{i}{R}_{j}}\right)\right]}^{1/2}


:math:`R_i` is the Born radius of particle *i*\ , which calculated as


.. math::
   R_i=\frac{1}{\rho_i^{-1}-r_i^{-1}\text{tanh}\left(\alpha \Psi_{i}-{\beta \Psi}_i^2+{\gamma \Psi}_i^3\right)}


where :math:`\alpha`, :math:`\beta`, and :math:`\gamma` are the GB\ :sup:`OBC`\ II parameters :math:`\alpha` = 1, :math:`\beta` = 0.8, and :math:`\gamma` =
4.85.  :math:`\rho_i` is the adjusted atomic radius of particle *i*\ , which
is calculated from the atomic radius :math:`r_i` as :math:`\rho_i = r_i-0.009` nm.
:math:`\Psi_i` is calculated as an integral over the van der Waals
spheres of all particles outside particle *i*\ :


.. math::
   \Psi_i=\frac{\rho_i}{4\pi}\int_{\text{VDW}}\theta\left(|\mathbf{r}|-{\rho }_{i}\right)\frac{1}{{|\mathbf{r}|}^{4}}{d}^{3}\mathbf{r}


where :math:`\theta`\ (\ *r*\ ) is a step function that excludes the interior of particle
\ *i* from the integral.

Surface Area Term
=================

The surface area term is given by\ :cite:`Schaefer1998`\ :cite:`Ponder`


.. math::
   E=E_{SA} \cdot 4\pi \sum _{i}{\left({r}_{i}+{r}_{\mathit{solvent}}\right)}^{2}{\left(\frac{{r}_{i}}{{R}_{i}}\right)}^{6}


where :math:`r_i` is the atomic radius of particle *i*\ , :math:`r_i` is
its atomic radius, and :math:`r_\mathit{solvent}` is the solvent radius, which is taken
to be 0.14 nm.  The default value for the energy scale :math:`E_{SA}` is 2.25936 kJ/mol/nm\ :sup:`2`\ .


GayBerneForce
*************

This is similar to the Lennard-Jones interaction described in section :ref:`lennard-jones-interaction`,
but instead of being based on the distance between two point particles, it is based
on the distance of closest approach between two ellipsoids.\ :cite:`Everaers2003`
Let :math:`\mathbf{A}_1` and :math:`\mathbf{A}_2` be rotation matrices that transform
from the lab frame to the body frames of two interacting ellipsoids.  These rotations
are determined from the positions of other particles, as described in the API documentation.
Let :math:`\mathbf{r}_{12}` be the vector pointing from particle 1 to particle 2, and
:math:`\hat{\mathbf{r}}_{12}=\mathbf{r}_{12}/|\mathbf{r}_{12}|`.  Let :math:`\mathbf{S}_1`
and :math:`\mathbf{S}_2` be diagonal matrices containing the three radii of each particle:

.. math::
   \mathbf{S}_i=\begin{bmatrix}
   a_i & 0 & 0 \\
   0 & b_i & 0 \\
   0 & 0 & c_i
   \end{bmatrix}

The energy is computed as a product of three terms:

.. math::
   E=U_r(\mathbf{A}_1, \mathbf{A}_2, \mathbf{r}_{12}) \cdot \eta_{12}(\mathbf{A}_1, \mathbf{A}_2) \cdot \chi_{12}(\mathbf{A}_1, \mathbf{A}_2, \hat{\mathbf{r}}_{12})

The first term describes the distance dependence, and is very similar in form to
the Lennard-Jones interaction:

.. math::
   U_r=4\epsilon\left({\left(\frac{\sigma}{h_{12}+\sigma}\right)}^{12}-{\left(\frac{\sigma}{h_{12}+\sigma}\right)}^{6}\right)

where :math:`h_{12}` is an approximation to the distance of closest approach between
the two ellipsoids:

.. math::
   h_{12}=|\mathbf{r}_{12}|-\sigma_{12}(\mathbf{A}_1, \mathbf{A}_2, \hat{\mathbf{r}}_{12})

.. math::
   \sigma_{12}(\mathbf{A}_1, \mathbf{A}_2, \hat{\mathbf{r}}_{12})=\left[ \frac{1}{2} \hat{\mathbf{r}}_{12}^T \mathbf{G}_{12}^{-1} \hat{\mathbf{r}}_{12} \right]^{-1/2}

.. math::
   \mathbf{G}_{12}=\mathbf{A}_1^T \mathbf{S}_1^2 \mathbf{A}_1 + \mathbf{A}_2^T \mathbf{S}_2^2 \mathbf{A}_2

The second term adjusts the energy based on the relative orientations of the two ellipsoids:

.. math::
   \eta_{12}(\mathbf{A}_1, \mathbf{A}_2)=\left[ \frac{2 s_1 s_2}{\text{det}(\mathbf{G}_{12})} \right]^{1/2}

.. math::
   s_i=(a_i b_i + c_i^2)\sqrt{a_i b_i}

The third term applies the user-defined scale factors :math:`e_a`, :math:`e_b`,
and :math:`e_c` that adjust the strength of the interaction along each axis:

.. math::
   \chi_{12}(\mathbf{A}_1, \mathbf{A}_2, \hat{\mathbf{r}}_{12})=(2 \hat{\mathbf{r}}_{12}^T \mathbf{B}_{12}^{-1} \hat{\mathbf{r}}_{12})^2

.. math::
   \mathbf{B}_{12}=\mathbf{A}_1^T \mathbf{E}_1 \mathbf{A}_1 + \mathbf{A}_2^T \mathbf{E}_2 \mathbf{A}_2

.. math::
   \mathbf{E}_i=\begin{bmatrix}
   e_{ai}^{-1/2} & 0 & 0 \\
   0 & e_{bi}^{-1/2} & 0 \\
   0 & 0 & e_{ci}^{-1/2}
   \end{bmatrix}

When using a cutoff, you can optionally use a switching function to make the energy go smoothly to 0
at the cutoff distance.  When :math:`r_\mathit{switch} < r < r_\mathit{cutoff}`\ , the energy is multiplied by

.. math::
   S=1-{6x}^{5}+15{x}^{4}-10{x}^{3}

where :math:`x = (r-r_\mathit{switch})/(r_\mathit{cutoff}-r_\mathit{switch})`. This function decreases smoothly from 1 at
:math:`r = r_\mathit{switch}` to 0 at :math:`r = r_\mathit{cutoff}`, and has continuous first and
second derivatives at both ends.


AndersenThermostat
******************

AndersenThermostat couples the system to a heat bath by randomly selecting a
subset of particles at the start of each time step, then setting their
velocities to new values chosen from a Boltzmann distribution.  This represents
the effect of random collisions between particles in the system and particles in
the heat bath.\ :cite:`Andersen1980`

The probability that a given particle will experience a collision in a given
time step is


.. math::
   P=1-{e}^{-f\Delta t}


where *f* is the collision frequency and :math:`\Delta t` is the step size.
Each component of its velocity is then set to


.. math::
   {v}_{i}=\sqrt{\frac{{k}_{B}T}{m}}R


where *T* is the thermostat temperature, *m* is the particle mass, and
*R* is a random number chosen from a normal distribution with mean of zero and
variance of one.

MonteCarloBarostat
******************

MonteCarloBarostat models the effect of constant pressure by allowing the size
of the periodic box to vary with time.\ :cite:`Chow1995`\ :cite:`Aqvist2004`
At regular intervals, it attempts a Monte Carlo step by scaling the box vectors
and the coordinates of each molecule’s center by a factor *s*\ .  The scale
factor *s* is chosen to change the volume of the periodic box from *V*
to *V*\ +\ :math:`\Delta`\ *V*\ :


.. math::
   s={\left(\frac{V+\Delta V}{V}\right)}^{1/3}


The change in volume is chosen randomly as


.. math::
   \Delta V=A\cdot r


where *A* is a scale factor and *r* is a random number uniformly
distributed between -1 and 1.  The step is accepted or rejected based on the
weight function


.. math::
   \Delta W=\Delta E+P\Delta V-Nk_{B}T \text{ln}\left(\frac{V+\Delta V}{V}\right)


where :math:`\Delta E` is the change in potential energy resulting from the step,
\ *P* is the pressure being applied to the system, *N* is the number of molecules in the
system, :math:`k_B` is Boltzmann’s constant, and *T* is the system
temperature.  In particular, if :math:`\Delta W\le 0` the step is always accepted.
If :math:`\Delta W > 0`\ , the step is accepted with probability
:math:`\text{exp}(-\Delta W/k_B T)`\ .

This algorithm tends to be more efficient than deterministic barostats such as
the Berendsen or Parrinello-Rahman algorithms, since it does not require an
expensive virial calculation at every time step.  Each Monte Carlo step involves
two energy evaluations, but this can be done much less often than every time
step.  It also does not require you to specify the compressibility of the
system, which usually is not known in advance.

The scale factor *A* that determines the size of the steps is chosen
automatically to produce an acceptance rate of approximately 50%.  It is
initially set to 1% of the periodic box volume.  The acceptance rate is then
monitored, and if it varies too much from 50% then *A* is modified
accordingly.

Each Monte Carlo step modifies particle positions by scaling the centroid of
each molecule, then applying the resulting displacement to each particle in the
molecule.  This ensures that each molecule is translated as a unit, so bond
lengths and constrained distances are unaffected.

MonteCarloBarostat assumes the simulation is being run at constant temperature
as well as pressure, and the simulation temperature affects the step acceptance
probability.  It does not itself perform temperature regulation, however.  You
must use another mechanism along with it to maintain the temperature, such as
LangevinIntegrator or AndersenThermostat.

MonteCarloAnisotropicBarostat
*****************************

MonteCarloAnisotropicBarostat is very similar to MonteCarloBarostat, but instead
of scaling the entire periodic box uniformly, each Monte Carlo step scales only
one axis of the box.  This allows the box to change shape, and is useful for
simulating anisotropic systems whose compressibility is different along
different directions.  It also allows a different pressure to be specified for
each axis.

You can specify that the barostat should only be applied to certain axes of the
box, keeping the other axes fixed.  This is useful, for example, when doing
constant surface area simulations of membranes.

MonteCarloMembraneBarostat
**************************

MonteCarloMembraneBarostat is very similar to MonteCarloBarostat, but it is
specialized for simulations of membranes.  It assumes the membrane lies in the
XY plane.  In addition to applying a uniform pressure to regulate the volume of
the periodic box, it also applies a uniform surface tension to regulate the
cross sectional area of the periodic box in the XY plane.  The weight function
for deciding whether to accept a step is

.. math::
   \Delta W=\Delta E+P\Delta V-S\Delta A-Nk_{B}T \text{ln}\left(\frac{V+\Delta V}{V}\right)

where *S* is the surface tension and :math:`\Delta`\ *A* is the change in cross
sectional area.  Notice that pressure and surface tension are defined with
opposite senses: a larger pressure tends to make the box smaller, but a larger
surface tension tends to make the box larger.

MonteCarloMembraneBarostat offers some additional options to customize the
behavior of the periodic box:

* The X and Y axes can be either

  * isotropic (they are always scaled by the same amount, so their ratio remains fixed)
  * anisotropic (they can change size independently)

* The Z axis can be either

  * free (its size changes independently of the X and Y axes)
  * fixed (its size does not change)
  * inversely varying with the X and Y axes (so the total box volume does not
    change)

CMMotionRemover
***************

CMMotionRemover prevents the system from drifting in space by periodically
removing all center of mass motion.  At the start of every *n*\ ’th time step
(where *n* is set by the user), it calculates the total center of mass
velocity of the system:


.. math::
   \mathbf{v}_\text{CM}=\frac{\sum _{i}{m}_{i}\mathbf{v}_{i}}{\sum _{i}{m}_{i}}


where :math:`m_i` and :math:`\mathbf{v}_i` are the mass and velocity of particle
\ *i*\ .  It then subtracts :math:`\mathbf{v}_\text{CM}` from the velocity of every
particle.

RMSDForce
*********

RMSDForce computes the root-mean-squared deviation (RMSD) between the current
particle positions :math:`\mathbf{x}_i` and a set of reference positions
:math:`\mathbf{x}_i^\text{ref}`:

.. math::
   \text{RMSD} = \sqrt{\frac{\sum_{i} \| \mathbf{x}_i - \mathbf{x}_i^\text{ref} \|^2}{N}}

Before computing this, the reference positions are first translated and rotated
so as to minimize the RMSD.  The computed value is therefore :math:`argmin(\text{RMSD})`,
where the :math:`argmin` is taken over all possible translations and rotations.

This force is normally used with a CustomCVForce (see Section :ref:`customcvforce`).
One rarely wants a force whose energy exactly equals the RMSD, but there are many
situations where it is useful to have a restraining or biasing force that depends
on the RMSD in some way.


Custom Forces
#############

In addition to the standard forces described in the previous chapter, OpenMM
provides a number of “custom” force classes.   These classes provide detailed
control over the mathematical form of the force by allowing the user to specify
one or more arbitrary algebraic expressions.  The details of how to write these
custom expressions are described in section :ref:`writing-custom-expressions`\ .

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
:class:`CustomIntegrator` (described in section :ref:`custom-integrator`\ ).  The
derivative can appear directly in expressions that define the integration
algorithm.  This can be used to implement algorithms such as lambda-dynamics,
where a global parameter is integrated as a dynamic variable.


Integrators
###########


VerletIntegrator
****************

VerletIntegrator implements the leap-frog Verlet integration method.  The
positions and velocities stored in the context are offset from each other by
half a time step.  In each step, they are updated as follows:


.. math::
   \mathbf{v}_{i}(t+\Delta t/2)=\mathbf{v}_{i}(t-\Delta t/2)+\mathbf{f}_{i}(t)\Delta t/{m}_{i}


.. math::
   \mathbf{r}_{i}(t+\Delta t)=\mathbf{r}_{i}(t)+\mathbf{v}_{i}(t+\Delta t/2)\Delta t


where :math:`\mathbf{v}_i` is the velocity of particle *i*\ , :math:`\mathbf{r}_i` is
its position, :math:`\mathbf{f}_i` is the force acting on it, :math:`m_i` is its
mass, and :math:`\Delta t` is the time step.

Because the positions are always half a time step later than the velocities,
care must be used when calculating the energy of the system.  In particular, the
potential energy and kinetic energy in a State correspond to different times,
and you cannot simply add them to get the total energy of the system.  Instead,
it is better to retrieve States after two successive time steps, calculate the
on-step velocities as


.. math::
   \mathbf{v}_{i}(t)=\frac{\mathbf{v}_{i}\left(t-\Delta t/2\right)+\mathbf{v}_{i}\left(t+\Delta t/2\right)}{2}


then use those velocities to calculate the kinetic energy at time *t*\ .

LangevinIntegator
*****************

LangevinIntegator simulates a system in contact with a heat bath by integrating
the Langevin equation of motion:


.. math::
   m_i\frac{d\mathbf{v}_i}{dt}=\mathbf{f}_i-\gamma m_i \mathbf{v}_i+\mathbf{R}_i


where :math:`\mathbf{v}_i` is the velocity of particle *i*\ , :math:`\mathbf{f}_i` is
the force acting on it, :math:`m_i` is its mass, :math:`\gamma` is the friction
coefficient, and :math:`\mathbf{R}_i` is an uncorrelated random force whose
components are chosen from a normal distribution with mean zero and variance
:math:`2m_i \gamma k_B T`\ , where *T* is the temperature of
the heat bath.

The integration is done using a leap-frog method similar to VerletIntegrator.
:cite:`Izaguirre2010` The same comments about the offset between positions and
velocities apply to this integrator as to that one.

BrownianIntegrator
******************

BrownianIntegrator simulates a system in contact with a heat bath by integrating
the Brownian equation of motion:


.. math::
   \frac{d\mathbf{r}_i}{dt}=\frac{1}{\gamma m_i}\mathbf{f}_i+\mathbf{R}_i


where :math:`\mathbf{r}_i` is the position of particle *i*\ , :math:`\mathbf{f}_i` is
the force acting on it, :math:`\gamma` is the friction coefficient, and :math:`\mathbf{R}_i`
is an uncorrelated random force whose components are chosen from a normal
distribution with mean zero and variance :math:`2 k_B T/m_i  \gamma`,
where *T* is the temperature of the heat bath.

The Brownian equation of motion is derived from the Langevin equation of motion
in the limit of large :math:`\gamma`\ .  In that case, the velocity of a particle is
determined entirely by the instantaneous force acting on it, and kinetic energy
ceases to have much meaning, since it disappears as soon as the applied force is
removed.


VariableVerletIntegrator
************************

This is very similar to VerletIntegrator, but instead of using the same step
size for every time step, it continuously adjusts the step size to keep the
integration error below a user-specified tolerance.  It compares the positions
generated by Verlet integration with those that would be generated by an
explicit Euler integrator, and takes the difference between them as an estimate
of the integration error:


.. math::
   error={\left(\Delta t\right)}^{2}\sum _{i}\frac{|\mathbf{f}_{i}|}{m_i}


where :math:`\mathbf{f}_i` is the force acting on particle *i* and :math:`m_i`
is its mass.  (In practice, the error made by the Euler integrator is usually
larger than that made by the Verlet integrator, so this tends to overestimate
the true error.  Even so, it can provide a useful mechanism for step size
control.)

It then selects the value of :math:`\Delta t` that makes the error exactly equal the
specified error tolerance:


.. math::
   \Delta t=\sqrt{\frac{\delta}{\sum _{i}\frac{|\mathbf{f}_i|}{m_i}}}


where :math:`\delta` is the error tolerance.  This is the largest step that may be
taken consistent with the user-specified accuracy requirement.

(Note that the integrator may sometimes choose to use a smaller value for :math:`\Delta t`
than given above.  For example, it might restrict how much the step size
can grow from one step to the next, or keep the step size constant rather than
increasing it by a very small amount.  This behavior is not specified and may
vary between Platforms.  It is required, however, that :math:`\Delta t` never be larger
than the value given above.)

A variable time step integrator is generally superior to a fixed time step one
in both stability and efficiency.  It can take larger steps on average, but will
automatically reduce the step size to preserve accuracy and avoid instability
when unusually large forces occur.  Conversely, when each uses the same step
size on average, the variable time step one will usually be more accurate since
the time steps are concentrated in the most difficult areas of the trajectory.

Unlike a fixed step size Verlet integrator, variable step size Verlet is not
symplectic.  This means that for a given average step size, it will not conserve
energy as precisely over long time periods, even though each local region of the
trajectory is more accurate.  For this reason, it is most appropriate when
precise energy conservation is not important, such as when simulating a system
at constant temperature.  For constant energy simulations that must maintain the
energy accurately over long time periods, the fixed step size Verlet may be more
appropriate.

VariableLangevinIntegrator
**************************

This is similar to LangevinIntegrator, but it continuously adjusts the step size
using the same method as VariableVerletIntegrator.  It is usually preferred over
the fixed step size Langevin integrator for the reasons given above.
Furthermore, because Langevin dynamics involves a random force, it can never be
symplectic and therefore the fixed step size Verlet integrator’s advantages do
not apply to the Langevin integrator.

.. _custom-integrator:

CustomIntegrator
****************

CustomIntegrator is a very flexible class that can be used to implement a wide
range of integration methods.  This includes both deterministic and stochastic
integrators; Metropolized integrators; multiple time step integrators; and
algorithms that must integrate additional quantities along with the particle
positions and momenta.

The algorithm is specified as a series of computations that are executed in
order to perform a single time step.  Each computation computes the value (or
values) of a *variable*\ .  There are two types of variables: *global
variables* have a single value, while *per-DOF variables* have a separate
value for every degree of freedom (that is, every *x*\ , *y*\ , or *z*
component of a particle).  CustomIntegrator defines lots of variables you can
compute and/or use in computing other variables.  Some examples include the step
size (global), the particle positions (per-DOF), and the force acting on each
particle (per-DOF).  In addition, you can define as many variables as you want
for your own use.

The actual computations are defined by mathematical expressions as described in
section :ref:`writing-custom-expressions`\ .  Several types of computations are supported:

* *Global*\ : the expression is evaluated once, and the result is stored into
  a global variable.
* *Per-DOF*\ : the expression is evaluated once for every degree of freedom,
  and the results are stored into a per-DOF variable.
* *Sum*\ : the expression is evaluated once for every degree of freedom.  The
  results for all degrees of freedom are added together, and the sum is stored
  into a global variable.


There also are other, more specialized types of computations that do not involve
mathematical expressions.  For example, there are computations that apply
distance constraints, modifying the particle positions or velocities
accordingly.

CustomIntegrator is a very powerful tool, and this description only gives a
vague idea of the scope of its capabilities.  For full details and examples,
consult the API documentation.


.. _other-features:

Other Features
##############


Periodic Boundary Conditions
****************************

Many Force objects support periodic boundary conditions.  They act as if space
were tiled with infinitely repeating copies of the system, then compute the
forces acting on a single copy based on the infinite periodic copies.  In most
(but not all) cases, they apply a cutoff so that each particle only interacts
with a single copy of each other particle.

OpenMM supports triclinic periodic boxes.  This means the periodicity is defined
by three vectors, :math:`\mathbf{a}`\ , :math:`\mathbf{b}`\ , and
:math:`\mathbf{c}`\ .  Given a particle position, the infinite periodic copies
of that particle are generated by adding vectors of the form
:math:`i \mathbf{a}+j \mathbf{b}+k \mathbf{c}`\ , where :math:`i`\ ,
:math:`j`\ , and :math:`k` are arbitrary integers.

The periodic box vectors must be chosen to satisfy certain requirements.
Roughly speaking, :math:`\mathbf{a}`\ , :math:`\mathbf{b}`\ , and
:math:`\mathbf{c}` need to "mostly" correspond to the x, y, and z axes.  They
must have the form

.. math::
   \mathbf{a} = (a_x, 0, 0)

   \mathbf{b} = (b_x, b_y, 0)

   \mathbf{c} = (c_x, c_y, c_z)

It is always possible to put the box vectors into this form by rotating the
system until :math:`\mathbf{a}` is parallel to x and :math:`\mathbf{b}` lies in
the xy plane.

Furthermore, they must obey the following constraints:

.. math::
   a_x > 0, b_y > 0, c_z > 0

   a_x \ge 2 |b_x|

   a_x \ge 2 |c_x|

   b_y \ge 2 |c_y|

This effectively requires the box vectors to be specified in a particular
reduced form.  By forming combinations of box vectors (a process known as
"lattice reduction"), it is always possible to put them in this form without
changing the periodic system they represent.

These requirements have an important consequence: the periodic unit cell can
always be treated as an axis-aligned rectangular box of size
:math:`(a_x, b_y, c_z)`\ .  The remaining non-zero elements of the box vectors
cause the repeating copies of the system to be staggered relative to each other,
but they do not affect the shape or size of each copy.  The volume of the unit
cell is simply given by :math:`a_x b_y c_z`\ .

LocalEnergyMinimizer
********************

This provides an implementation of the L-BFGS optimization algorithm.
:cite:`Liu1989`  Given a Context specifying initial particle positions, it
searches for a nearby set of positions that represent a local minimum of the
potential energy.  Distance constraints are enforced during minimization by
adding a harmonic restraining force to the potential function.  The strength of
the restraining force is steadily increased until the minimum energy
configuration satisfies all constraints to within the tolerance specified by the
Context's Integrator.

XMLSerializer
*************

This provides the ability to “serialize” a System, Force, Integrator, or State
object to a portable XML format, then reconstruct it again later.  When
serializing a System, the XML data contains a complete copy of the entire system
definition, including all Forces that have been added to it.

Here are some examples of uses for this class:

#. A model building utility could generate a System in memory, then serialize it
   to a file on disk.  Other programs that perform simulation or analysis could
   then reconstruct the model by simply loading the XML file.
#. When running simulations on a cluster, all model construction could be done
   on a single node.  The Systems and Integrators could then be encoded as XML,
   allowing them to be easily transmitted to other nodes.


XMLSerializer is a templatized class that, in principle, can be used to
serialize any type of object.  At present, however, only System, Force,
Integrator, and State are supported.

Force Groups
************

It is possible to split the Force objects in a System into groups.  Those groups
can then be evaluated independently of each other.  Some Force classes also
provide finer grained control over grouping.  For example, NonbondedForce allows
direct space computations to be in one group and reciprocal space computations
in a different group.

The most important use of force groups is for implementing multiple time step
algorithms with CustomIntegrator.  For example, you might evaluate the slowly
changing nonbonded interactions less frequently than the quickly changing bonded
ones.  It also is useful if you want the ability to query a subset of the forces
acting on the system.

Virtual Sites
*************

A virtual site is a particle whose position is computed directly from the
positions of other particles, not by integrating the equations of motion.  An
important example is the “extra sites” present in 4 and 5 site water models.
These particles are massless, and therefore cannot be integrated.  Instead,
their positions are computed from the positions of the massive particles in the
water molecule.

Virtual sites are specified by creating a VirtualSite object, then telling the
System to use it for a particular particle.  The VirtualSite defines the rules
for computing its position.  It is an abstract class with subclasses for
specific types of rules.  They are:

* TwoParticleAverageSite: The virtual site location is computed as a weighted
  average of the positions of two particles:

.. math::
   \mathbf{r}={w}_{1}\mathbf{r}_{1}+{w}_{2}\mathbf{r}_{2}

* ThreeParticleAverageSite: The virtual site location is computed as a weighted
  average of the positions of three particles:

.. math::
   \mathbf{r}={w}_{1}\mathbf{r}_{1}+{w}_{2}\mathbf{r}_{2}+{w}_{3}\mathbf{r}_{3}

* OutOfPlaneSite: The virtual site location is computed as a weighted average
  of the positions of three particles and the cross product of their relative
  displacements:

.. math::
   \mathbf{r}=\mathbf{r}_{1}+{w}_{12}\mathbf{r}_{12}+{w}_{13}\mathbf{r}_{13}+{w}_\mathit{cross}\left(\mathbf{r}_{12}\times \mathbf{r}_{13}\right)
..

  where :math:`\mathbf{r}_{12} = \mathbf{r}_{2}-\mathbf{r}_{1}` and
  :math:`\mathbf{r}_{13} = \mathbf{r}_{3}-\mathbf{r}_{1}`\ .  This allows
  the virtual site to be located outside the plane of the three particles.

* LocalCoordinatesSite: The locations of several other particles are used to compute a local
  coordinate system, and the virtual site is placed at a fixed location in that coordinate
  system.  The number of particles used to define the coordinate system is user defined.
  The origin of the coordinate system and the directions of its x and y axes
  are each specified as a weighted sum of the locations of the other particles:

.. math::
   \mathbf{o}={w}^{o}_{1}\mathbf{r}_{1} + {w}^{o}_{2}\mathbf{r}_{2} + ...

   \mathbf{dx}={w}^{x}_{1}\mathbf{r}_{1} + {w}^{x}_{2}\mathbf{r}_{2} + ...

   \mathbf{dy}={w}^{y}_{1}\mathbf{r}_{1} + {w}^{y}_{2}\mathbf{r}_{2} + ...

   \mathbf{dz}=\mathbf{dx}\times \mathbf{dy}
..

   These vectors are then used to construct a set of orthonormal coordinate axes as follows:

.. math::
   \mathbf{\hat{x}}=\mathbf{dx}/|\mathbf{dx}|

   \mathbf{\hat{z}}=\mathbf{dz}/|\mathbf{dz}|

   \mathbf{\hat{y}}=\mathbf{\hat{z}}\times \mathbf{\hat{x}}
..

   Finally, the position of the virtual site is set to

.. math::
   \mathbf{r}=\mathbf{o}+p_1\mathbf{\hat{x}}+p_2\mathbf{\hat{y}}+p_3\mathbf{\hat{z}}
..

Random Numbers with Stochastic Integrators and Forces
*****************************************************

OpenMM includes many stochastic integrators and forces that make extensive use
of random numbers. It is impossible to generate truly random numbers on a
computer like you would with a dice roll or coin flip in real life---instead
programs rely on pseudo-random number generators (PRNGs) that take some sort of
initial "seed" value and steps through a sequence of seemingly random numbers.

The exact implementation of the PRNGs is not important (in fact, each platform
may have its own PRNG whose performance is optimized for that hardware).  What
*is* important, however, is that the PRNG must generate a uniform distribution
of random numbers between 0 and 1. Random numbers drawn from this distribution
can be manipulated to yield random integers in a desired range or even a random
number from a different type of probability distribution function (e.g., a
normal distribution).

What this means is that the random numbers used by integrators and forces within
OpenMM cannot have any discernible pattern to them.  Patterns can be induced in
PRNGs in two principal ways:

1. The PRNG uses a bad algorithm with a short period.
2. Two PRNGs are started using the same seed

All PRNG algorithms in common use are periodic---that is their sequence of
random numbers repeats after a given *period*, defined by the number of "unique"
random numbers in the repeating pattern.  As long as this period is longer than
the total number of random numbers your application requires (preferably by
orders of magnitude), the first problem described above is avoided. All PRNGs
employed by OpenMM have periods far longer than any current simulation can cycle
through.

Point two is far more common in biomolecular simulation, and can result in very
strange artifacts that may be difficult to detect. For example, with Langevin
dynamics, two simulations that use the same sequence of random numbers appear to
synchronize in their global movements.\ :cite:`Uberuaga2004`\
:cite:`Sindhikara2009` It is therefore very important that the stochastic forces
and integrators in OpenMM generate unique sequences of pseudo-random numbers not
only within a single simulation, but between two different simulations of the
same system as well (including any restarts of previous simulations).

Every stochastic force and integrator that does (or could) make use of random
numbers has two instance methods attached to it: :meth:`getRandomNumberSeed()`
and :meth:`setRandomNumberSeed(int seed)`. If you set a unique random seed for
two different simulations (or different forces/integrators if applicable),
OpenMM guarantees that the generated sequences of random numbers will be
different (by contrast, no guarantee is made that the same seed will result in
identical random number sequences).

Since breaking simulations up into pieces and/or running multiple replicates of
a system to obtain more complete statistics is common practice, a new strategy
has been employed for OpenMM versions 6.3 and later with the aim of trying to
ensure that each simulation will be started with a unique random seed. A random
seed value of 0 (the default) will cause a unique random seed to be generated
when a new :class:`Context` is instantiated.

Prior to the introduction of this feature, deserializing a serialized
:class:`System` XML file would result in each stochastic force or integrator
being assigned the same random seed as the original instance that was
serialized. If you use a :class:`System` XML file generated by a version of
OpenMM older than 6.3 to start a new simulation, you should manually set the
random number seed of each stochastic force or integrator to 0 (or another
unique value).

.. _standard-forces:

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

The reaction field approximation is not applied to nonbonded exceptions.  They
are always evaluated at full strength, regardless of the cutoff distance.  That
is because exceptions are primarily used to model 1-4 interactions, which are
really a type of bonded interaction.  They are typically parametrized without a
cutoff together with the other bonded interactions.

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

.. _coulomb-interaction-pme:

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

ConstantPotentialForce
**********************

This force is an implementation of the periodic finite field constant potential
method :cite:`Dufils2019`.  The constant potential method implements Coulomb
interactions between particles using the particle mesh Ewald method as described
in section :numref:`coulomb-interaction-pme`.  Unlike Coulomb interactions
implemented by a NonbondedForce that use fixed charges, the constant potential
method allows certain particles in a system to be designated electrode particles
whose charges can fluctuate.  Additionally, instead of point charges, these
electrode particles have Gaussian charge densities, are held at given electric
potentials, and can optionally employ a simple Thomas-Fermi metallicity model
:cite:`Scalfi2020`.  Finally, an electric field can be applied to all charged
particles in the system.  For a given configuration of all particles, the
charges on electrode particles are set to minimize the total electrostatic
energy, given by

.. math::
   E=E_{PME}+E_{gauss}+E_{self}+E_{field}+E_{potential}+E_{TF}

Here :math:`E_{PME}` is the interaction of point particles computed by the
particle mesh Ewald method.  Now

.. math::
   E_{gauss}=-\frac{1}{2}\sum_{i,j}q_iq_i\frac{\text{erfc}\left(\eta_{ij}r_{ij}\right)}{r_{ij}}

with :math:`\eta_{ij}=1/\sqrt{w_i^2+w_j^2}`, where :math:`w_i=0` for a
non-electrode (point) particle and :math:`w_i>0` for an electrode particle.  No
term is included in the sum when :math:`i=j` or :math:`w_i=w_j=0`.  Note that in
OpenMM's implementation, specification of an electrode requires the Gaussian
width :math:`w_i`, while the reciprocal width :math:`\eta_i=1/w_i` is often used
in the literature.  In addition to this correction to interactions between
particles, the presence of Gaussian charge densities on electrode particles
contributes a self-interaction energy:

.. math::
   E_{self}=\frac{1}{\sqrt{2\pi}}\sum_{i\in\text{elec}}\frac{q_i^2}{w_i}

Next, for the electric field :math:`\mathbf{E}` (on all particles) and applied
potential (:math:`\Psi_i` on electrode particles :math:`i`),

.. math::
   E_{field}=-\sum_{i}q_i\mathbf{r}_i\cdot\mathbf{E}

.. math::
   E_{potential}=-\sum_{i\in\text{elec}}q_i\Psi_i

Finally, for the Thomas-Fermi contribution on the electrode particles,

.. math::
   E_{TF}=2\pi\sum_{i\in\text{elec}}q_i^2\left(\frac{\ell_{TF,i}^2}{v_{TF,i}}\right)

where :math:`\ell_{TF,i}` is the Thomas-Fermi length and :math:`v_{TF,i}` is a
characteristic volume scale (often referred to as a "Voronoi volume" in other
implementations).  In OpenMM's implementation, the parameter
:math:`\ell_{TF,i}^2/v_{TF,i}` is specified as a single Thomas-Fermi parameter,
with units of reciprocal length, for all particles in an electrode.

Besides these additional terms, and the fluctuating nature of the electrode
particle charges, there are a few other important differences between the
implementation of nonbonded interactions in ConstantPotentialForce and that in
NonbondedForce using PME:

1. ConstantPotentialForce does not include any capability for computing
   Lennard-Jones interactions.  If Lennard-Jones interactions (as NonbondedForce
   computes them) are desired, a separate NonbondedForce should be added to the
   system with the appropriate sigma and epsilon parameters set, and with all
   charges set to zero.

2. For the solved charges on electrode particles to be valid, all particles in
   the system must use a single ConstantPotentialForce.  Setting charges in more
   than one ConstantPotentialForce, or setting any non-zero charges in a
   NonbondedForce, will produce invalid results, as solving for the electrode
   charges requires a global minimization of the system's electrostatic energy.

3. The Gaussian charge correction to pairwise interactions given by
   :math:`E_{gauss}` above is calculated using the minimum image convention and
   truncated with a fixed cutoff :math:`r_{cut}`.  If decreasing :math:`r_{cut}`
   to tune PME performance, ensure it stays large enough with respect to the
   largest Gaussian width :math:`w_i` in the system.  Using the formula for the
   direct space error in section :numref:`coulomb-interaction-pme` with
   :math:`1/w_i` in place of :math:`\alpha` (owing to the similarities in the
   functional forms of these terms), :math:`r_{cut}\ge2.7w_{max}` should be
   sufficient for the default error tolerance :math:`\delta=5\times10^{-4}`.

ConstantPotentialForce provides two algorithms to solve for electrode charges.
The matrix solver precomputes a capacitance matrix for the system before a
simulation and solves directly for charges, while the conjugate gradient (CG)
solver finds charges iteratively without requiring a matrix inversion.  The
matrix solver may be faster if the number of electrode particles is small
enough, although it cannot be used unless the electrode particles are fixed in
place during a simulation by setting their masses to zero.  Unless the positions
of electrode particles are allowed to fluctuate (in which case only the CG
solver can be used), both solvers can be benchmarked on a given problem to find
the one with the best performance.

By default, the CG solver uses a preconditioner that automatically activates if
any electrodes use Gaussian widths or Thomas-Fermi parameters different from
each other.  This will usually allow the solver to converge to a given tolerance
in fewer iterations.  However, if these electrode parameters differ but are very
close to each other, benefit from the preconditioner may be offset by additional
numerical error introduced due to floating-point roundoff.  This should not
affect the results, but may cause the CG solver to require more iterations to
converge.  Disabling the preconditioner may improve performance in such a case.

Because ConstantPotentialForce needs to compute all of the electrostatic
interactions in a system, it should not be used with AmoebaMultipoleForce, and
is thus incompatible with induced dipole-based polarizable force fields.
However, ConstantPotentialForce can be used with Drude polarizable force fields.
None of the Drude sites or their parent particles should belong to an electrode.
Also, DrudeSCFIntegrator should not be used with ConstantPotentialForce to
integrate a system with Drude particles, because the current implementations
do not permit simultaneously solving for the charges on electrode particles and
the positions of the Drude particles.  Constant potential simulations with Drude
polarizable force fields for the electrolyte should use DrudeLangevinIntegrator
or DrudeNoseHooverIntegrator instead.

GayBerneForce
*************

This is similar to the Lennard-Jones interaction described in section :numref:`lennard-jones-interaction`,
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

Another feature of MonteCarloBarostat is that it can compute the "instantaneous
pressure", which is defined based on the molecular virial:

.. math::
   P_{inst} = \frac{1}{3V} (\sum_i m_i |\mathbf{v}_i|^2 ) - \frac{dE}{dV}

where the sum is taken over molecules, :math:`m_i` is the mass of the i'th molecule,
and :math:`\mathbf{v}_i` is the velocity of its center of mass.  The derivative
of potential energy with respect to volume is approximated with a finite difference.

In most cases, the time average of the instantaneous pressure should equal the
pressure applied by the barostat.  Fluctuations around the average can be extremely
large, however, especially when simulating incompressible materials like water.
A very long simulation may be required to accurately compute the average.  There
also are situations where the average instantaneous pressure differs from the
applied pressure.  For example, if the system contains immobile massless particles,
they will reduce the kinetic energy below what would be expected based on the
temperature, and hence reduce the calculated pressure.

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

Like MonteCarloBarostat, the anisotropic barostat can compute an instantaneous
pressure, but in this case the pressure is computed separately along each axis
according to the formula

.. math::
   P_{inst} = \frac{1}{V} \sum_i m_i v_i^2 - \frac{dE}{dV}

where :math:`v_i` is the component of the i'th molecule's velocity along the
specified axis.  The derivative :math:`dE/dV` is again computed with a finite
difference, but in this case :math:`dE` refers to the change in potential energy
when scaling the box size along only a single axis.

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

Like other barostats, MonteCarloMembraneBarostat can compute the instantaneous
pressure.  It is computed separately for each axis using the same formula as
MonteCarloAnisotropicBarostat.

MonteCarloFlexibleBarostat
**************************

MonteCarloFlexibleBarostat is very similar to MonteCarloBarostat, but it allows
the periodic box to be fully flexible.\ :cite:`Vandenhaute2021`  Monte Carlo
steps can change not just the lengths of the box sides, but also the angles.  It
is especially useful for simulations of bulk materials where the shape of a
crystal's unit cell may not be known in advance, or could even change with time
as it transitions between phases.

Computing the instantaneous pressure with MonteCarloFlexibleBarostat is more
complicated than other barostats.  Because the box angles can change, it needs
to compute the full pressure tensor.  Let :math:`\mathbf{a}`, :math:`\mathbf{b}`,
and :math:`\mathbf{c}` be the vectors defining the periodic box.  They can be
assembled into a 3 by 3 box matrix

.. math::

    h = \left[ {\begin{array}{ccc}
        a_x & a_y & a_z \\
        b_x & b_y & b_z \\
        c_x & c_y & c_z \\
      \end{array} } \right]

The pressure tensor is similarly represented by a 3 by 3 matrix.  Its elements
are given by

.. math::

   P_{jk} = \frac{1}{V} \left( \sum_i m_i v_{ij} v_{ik} - \frac{\partial E}{\partial h_{jk}} h_{jk} \right)

Because of the particular reduced form OpenMM uses for the box vectors (see Section
:numref:`periodic-boundary-conditions`), the elements above the diagonal are
always zero.  In addition, the elements below the diagonal do not involve any
change to the box volume, and therefore are not affected by the applied pressure.
The six non-zero elements are as follows.

* The diagonal elements reflect the pressure acting to change the box size along
  each axis.  In equilibrium they should average to the applied pressure, although
  fluctuations around the average can be very large.
* The elements below the diagonal reflect the pressure acting to change the box
  angles.  In equilibrium they should average to zero.

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

This force is normally used with a CustomCVForce (see Section :numref:`customcvforce`).
One rarely wants a force whose energy exactly equals the RMSD, but there are many
situations where it is useful to have a restraining or biasing force that depends
on the RMSD in some way.

OrientationRestraintForce
*************************

OrientationRestraintForce is used to keep a group of particles in a fixed
orientation.  The calculation is closely related to the one in RMSD force.  Both
forces begin by finding the translation and rotation that optimally superimpose
the current particle positions on a set of reference positions.  Whereas
RMSDForce computes the deviation after removing the translation and rotation,
OrientationRestraintForce applies a force based on the rotation itself:

.. math::
   E = 2 k \mathrm{sin}^2(\theta/2)

where :math:`k` is the force constant and :math:`\theta` is the rotation angle.

For small rotations, :math:`E \approx \frac{k}{2}\theta^2`, giving a harmonic
restraint.  For large rotations, the restraint is weaker than harmonic, and the
force vanishes at :math:`\theta = \frac{\pi}{2}`.  This ensures that the force
is continuous everywhere; a simple harmonic restraint would have a discontinuous
force at :math:`\frac{\pi}{2}`.

RGForce
*********

RGForce computes the radius of gyration (Rg) of a set of particle positions:

.. math::
   \text{Rg} = \sqrt{\frac{\sum_{i} \| \mathbf{x}_i - \mathbf{x}_c \|^2}{N}}

where :math:`\mathbf{x}_c` is the center position,

.. math::
   \mathbf{x}_c = \frac{\sum_{i} \mathbf{x}_i}{N}

This force is normally used with a CustomCVForce (see Section :numref:`customcvforce`).
One rarely wants a force whose energy exactly equals Rg, but there are many
situations where it is useful to have a restraining or biasing force that depends
on Rg in some way.

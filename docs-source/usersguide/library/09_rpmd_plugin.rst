.. _ring-polymer-molecular-dynamics-plugin:

Ring Polymer Molecular Dynamics (RPMD) Plugin
#############################################

Ring Polymer Molecular Dynamics (RPMD) provides an efficient approach to include
nuclear quantum effects in molecular simulations.\ :cite:`Craig2004`  When
used to calculate static equilibrium properties, RPMD reduces to path integral
molecular dynamics and gives an exact description of the effect of quantum
fluctuations for a given potential energy model.\ :cite:`Parrinello1984`  For
dynamical properties RPMD is no longer exact but has shown to be a good
approximation in many cases.

For a system with a classical potential energy *E*\ (\ *q*\ ), the RPMD
Hamiltonian is given by


.. math::
   H=\sum _{k=1}^{n}\left(\frac{{p}_{{k}^{2}}}{2m}+E({q}_{k})+\frac{m({k}_{B}Tn)^{2}}{2\hbar^{2}}({q}_{k}-{q}_{k-1})^{2}\right)


This Hamiltonian resembles that of a system of classical ring polymers where
different copies of the system are connected by harmonic springs.  Hence each
copy of the classical system is commonly referred to as a “bead”.  The spread of
the ring polymer representing each particle is directly related to its De
Broglie thermal wavelength (uncertainty in its position).

RPMD calculations must be converged with respect to the number *n* of beads
used.  Each bead is evolved at the effective temperature *nT*\ , where *T*
is the temperature for which properties are required.  The number of beads
needed to converge a calculation can be estimated using\ :cite:`Markland2008`


.. math::
   n>\frac{\hbar\omega_{max}}{{k}_{B}T}


where :math:`\omega_{max}` is the highest frequency in the problem.  For example, for
flexible liquid water the highest frequency is the OH stretch at around 3000
cm\ :sup:`-1`\ , so around 24 to 32 beads are needed depending on the accuracy
required.  For rigid water where the highest frequency is only around 1000
cm\ :sup:`-1`\ , only 6 beads are typically needed.  Due to the replication needed
of the classical system, the extra cost of the calculation compared to a
classical simulation increases linearly with the number of beads used.

This cost can be reduced by “contracting” the ring polymer to a smaller number
of beads.\ :cite:`Markland2008`  The rapidly changing forces are then computed
for the full number of beads, while slower changing forces are computed on a
smaller set.  In the case of flexible water, for example, a common arrangement
would be to compute the high frequency bonded forces on all 32 beads, the direct
space nonbonded forces on only 6 beads, and the reciprocal space nonbonded
forces on only a single bead.

Due to the stiff spring terms between the beads, NVE RPMD trajectories can
suffer from ergodicity problems and hence thermostatting is highly recommended,
especially when dynamical properties are not required.\ :cite:`Hall1984`  The
thermostat implemented here is the path integral Langevin equation (PILE)
approach.\ :cite:`Ceriotti2010`  This method couples an optimal white noise
Langevin thermostat to the normal modes of each polymer, leaving only one
parameter to be chosen by the user which controls the friction applied to the
center of mass of each ring polymer.  A good choice for this is to use a value
similar to that used in a classical calculation of the same system.


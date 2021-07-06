.. _drude-plugin:

Drude Plugin
############

Drude oscillators are a method for incorporating electronic polarizability into
a model.\ :cite:`Lamoureux2003`  For each polarizable particle, a second
particle (the “Drude particle”) is attached to it by an anisotropic harmonic
spring.  When both particles are at the same location, they are equivalent to an
ordinary point particle.  Applying an electric field causes the Drude particle
to move a short distance away from its parent particle, creating an induced
dipole moment.  The polarizability :math:`\alpha` is related to the charge *q* on
the Drude particle and the spring constant *k* by



.. math::
   \alpha =\frac{{q}^{2}}{k}


A damped interaction\ :cite:`Thole1981` is used between dipoles that are
bonded to each other.

The equations of motion can be integrated with three different methods:

#. In the Self Consistent Field (SCF) method, the ordinary particles are first
   updated as usual.  A local energy minimization is then performed to select new
   positions for the Drude particles.  This ensures that the induced dipole moments
   respond instantly to changes in their environments.  This method is accurate but
   computationally expensive.
#. In the extended Lagrangian method, the positions of the Drude particles are
   treated as dynamical variables, just like any other particles.  A small amount
   of mass is transferred from the parent particles to the Drude particles,
   allowing them to be integrated normally.  A dual Langevin or Nose-Hoover integrator is used to
   maintain the center of mass of each Drude particle pair at the system
   temperature, while using a much lower temperature for their relative internal
   motion.  In practice, this produces dipole moments very close to those from the
   SCF solution while being much faster to compute.
#. The Nosé-Hoover dual thermostat method.  In this approach the motion of
   non-Drude sites and center of mass motion of Drude pairs are thermostated to
   the target temperature with one thermostat.  Another thermostat is used to keep
   relative motion of Drude pairs to a different, typically much lower,
   temperature to maintain separation of nuclear and electronic degrees of
   freedom.  The minimal specification is as follows::

      DrudeNoseHooverIntegrator integrator(temperature, frequency,
                                           temperatureDrude, frequencyDrude,
                                           1*femtoseconds)

   Where the first and third arguments specify the center-of-mass temperature and
   relative temperature for each Drude pair, respecitvely.  The second and fourth
   arguments describe the frequency of interaction with the center-of-mass and
   relative heat baths, respectively, and should be specified with inverse time
   units.  The fifth argument is the timestep.  The multi-timestep and Nosé-Hoover
   chain length may also be specified, but sensible defaults are provided.

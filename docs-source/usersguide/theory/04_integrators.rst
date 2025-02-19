.. _integrators-theory:

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
:math:`2 \gamma k_B T`. *T* is the temperature of the heat bath.

The integration is done using the LFMiddle discretization. :cite:`Zhang2019` In each step,
the positions and velocities are updated as follows:


.. math::
   \mathbf{v}_{i}(t+\Delta t/2) = \mathbf{v}_{i}(t-\Delta t/2) + \mathbf{f}_{i}(t)\Delta t/{m}_{i}


.. math::
   \mathbf{r}_{i}(t+\Delta t/2) = \mathbf{r}_{i}(t) + \mathbf{v}_{i}(t+\Delta t/2)\Delta t/2


.. math::
   \mathbf{v'}_{i}(t+\Delta t/2) = \mathbf{v}_{i}(t+\Delta t/2)\alpha + \sqrt{kT(1-\alpha^2)/m}R


.. math::
   \mathbf{r}_{i}(t+\Delta t) = \mathbf{r}_{i}(t+\Delta t/2) + \mathbf{v'}_{i}(t+\Delta t/2)\Delta t/2


where :math:`k` is Boltzmann's constant, :math:`T` is the temperature,
:math:`\gamma` is the friction coefficient, :math:`R` is a normally distributed
random number, and :math:`\alpha=\exp(-\gamma\Delta t)`.

The same comments about the offset between positions and velocities apply to
this integrator as to VerletIntegrator.

LangevinMiddleIntegrator
************************

This integrator is identical to LangevinIntegerator.  Separate classes exist only for historical reasons.

.. _nosehoover-integrators-theory:

NoseHooverIntegrator
********************

Like LangevinMiddleIntegerator, this uses the LFMiddle discretization.
:cite:`Zhang2019` In each step, the positions and velocities are updated as
follows:


.. math::
   \mathbf{v}_{i}(t+\Delta t/2) = \mathbf{v}_{i}(t-\Delta t/2) + \mathbf{f}_{i}(t)\Delta t/{m}_{i}


.. math::
   \mathbf{r}_{i}(t+\Delta t/2) = \mathbf{r}_{i}(t) + \mathbf{v}_{i}(t+\Delta t/2)\Delta t/2


.. math::
   \mathbf{v'}_{i}(t+\Delta t/2) = \mathrm{scale}\times\mathbf{v}_{i}(t+\Delta t/2)


.. math::
   \mathbf{r}_{i}(t+\Delta t) = \mathbf{r}_{i}(t+\Delta t/2) + \mathbf{v'}_{i}(t+\Delta t/2)\Delta t/2


The universal scale factor used in the third step is determined by propagating
auxilliary degrees of freedom alongside the regular particles.  The original
Nosé-Hoover formulation used a single harmonic oscillator for the heat bath,
but this is problematic in small or stiff systems, which are non-ergodic, so
the chain formulation extends this by replacing the single oscillator
thermostat with a chain of connected oscillators.  :cite:`Martyna1992`  For
large systems a single oscillator (*i.e.* a chain length of one) will suffice,
but longer chains are necessary to properly thermostat non-ergodic systems.
The OpenMM default is to use a chain length of three to cover the latter case,
but this can be safely reduced to increase efficiency in large systems.

The heat bath propagation is performed using a multi-timestep algorithm.  Each
propagation step is discretized into substeps using a factorization from
Yoshida and Suzuki; the default discretization uses a :math:`\mathcal{O}(\Delta
t^6)` approach that uses 7 points, but 1, 3 or 5 points may also be used to
increase performace, at the expense of accuracy.  Each step is further
subdivided into multi-timesteps with a default of 3 multi time steps per
propagation; as with the number of Yoshida-Suziki points this value may be
increase to increase accuracy but with additional computational expense.

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
section :numref:`writing-custom-expressions`\ .  Several types of computations are supported:

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

DPDIntegrator
*************

This implements dissipative particle dynamics (DPD). :cite:`Espanol1995`  It is
similar to a Langevin integrator, but instead of applying friction and noise to
the Cartesian coordinates of particles, it applies them to inter-particle
displacements.  This allows it to conserve momentum and efficiently model the
hydrodynamics of coarse grained models.

The system evolves according to the equation of motion

.. math::
   m_i\frac{d\mathbf{v}_i}{dt}=\mathbf{f}_i + \sum_{j \neq i } \mathbf{e}_{ij} [ -\gamma \omega^2(r_{ij})(\mathbf{e}_{ij} \cdot \mathbf{v}_{ij}) + \sqrt{2 \gamma k_BT} \omega(r_{ij}) R_{ij} ]


where :math:`\mathbf{v}_i` is the velocity of particle *i*\ , :math:`\mathbf{f}_i` is
the force acting on it, :math:`m_i` is its mass, :math:`r_{ij}` is the distance
between particles *i* and *j*\ , :math:`\mathbf{e}_{ij}` is a unit vector pointing from
particle *i* to particle *j*\ , :math:`\gamma` is the friction coefficient, and
:math:`\mathbf{R}_{ij} = \mathbf{R}_{ji}` is an uncorrelated random force whose
components are chosen from a normal distribution with mean zero and unit variance.
*T* is the temperature of the heat bath.  The friction and noise are limited to
pairs closer than a cutoff distance :math:`r_c` using the function
:math:`\omega(r) = 1-r/r_c`.  The friction coefficient :math:`\gamma` and cutoff
distance :math:`r_c` may be constants, or alternatively they can vary depending
on the types of the particles.

The integration is done using the same LFMiddle discretization used for
LangevinIntegrator. :cite:`Zhang2019`

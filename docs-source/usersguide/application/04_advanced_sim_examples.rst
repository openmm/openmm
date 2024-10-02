.. default-domain:: py

.. _advanced-simulation-examples:

Advanced Simulation Examples
############################

In the previous chapter, we looked at some basic scripts for running simulations
and saw lots of ways to customize them.  If that is all you want to do—run
straightforward molecular simulations—you already know everything you need to
know.  Just use the example scripts and customize them in the ways described in
Section :numref:`simulation-parameters`.

OpenMM can do far more than that.  Your script has the full OpenMM API at its
disposal, along with all the power of the Python language and libraries.  In
this chapter, we will consider some examples that illustrate more advanced
techniques.  Remember that these are still only examples; it would be impossible
to give an exhaustive list of everything OpenMM can do.  Hopefully they will
give you a sense of what is possible, and inspire you to experiment further on
your own.

Starting in this section, we will assume some knowledge of programming, as well
as familiarity with the OpenMM API.  Consult this User's Guide and the OpenMM
API documentation if you are uncertain about how something works. You can also
usethe Python :code:`help` command.  For example,
::

    help(Simulation)

will print detailed documentation on the :class:`Simulation` class.

Simulated Annealing
*******************

Here is a very simple example of how to do simulated annealing.  The following
lines linearly reduce the temperature from 300 K to 0 K in 100 increments,
executing 1000 time steps at each temperature:

.. samepage::
    ::

        ...
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        for i in range(100):
            integrator.setTemperature(3*(100-i)*kelvin)
            simulation.step(1000)

    .. caption::

        :autonumber:`Example,simulated annealing`

This code needs very little explanation.  The loop is executed 100 times.  Each
time through, it adjusts the temperature of the :class:`LangevinMiddleIntegrator` and then
calls :code:`step(1000)` to take 1000 time steps.

Applying an External Force to Particles: a Spherical Container
**************************************************************

In this example, we will simulate a non-periodic system contained inside a
spherical container with radius 2 nm.  We implement the container by applying a
harmonic potential to every particle:

.. math::
    E(r) = \begin{cases}
           0          & r\le2\\
           100(r-2)^2 & r>2
           \end{cases}

where *r* is the distance of the particle from the origin, measured in nm.
We can easily do this using OpenMM’s :class:`CustomExternalForce` class.  This class
applies a force to some or all of the particles in the system, where the energy
is an arbitrary function of each particle’s (\ *x*\ , *y*\ , *z*\ )
coordinates.  Here is the code to do it:

.. samepage::
    ::

        ...
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic,
                nonbondedCutoff=1*nanometer, constraints=None)
        force = CustomExternalForce('100*max(0, r-2)^2; r=sqrt(x*x+y*y+z*z)')
        system.addForce(force)
        for i in range(system.getNumParticles()):
            force.addParticle(i, [])
        integrator = LangevinMiddleIntegrator(300*kelvin, 91/picosecond, 0.004*picoseconds)
        ...

    .. caption::

        :autonumber:`Example,spherical container`

The first thing it does is create a :class:`CustomExternalForce` object and add it to the
:class:`System`.  The argument to :class:`CustomExternalForce` is a mathematical expression
specifying the potential energy of each particle.  This can be any function of *x*\ ,
*y*\ , and *z* you want.  It also can depend on global or per-particle
parameters.  A wide variety of restraints, steering forces, shearing forces,
etc. can be implemented with this method.

Next it must specify which particles to apply the force to.  In this case, we
want it to affect every particle in the system, so we loop over them and call
:meth:`addParticle` once for each one.  The two arguments are the index of
the particle to affect, and the list of per-particle parameter values (an empty
list in this case).  If we had per-particle parameters, such as to make the
force stronger for some particles than for others, this is where we would
specify them.

Notice that we do all of this immediately after creating the :class:`System`.  That is
not an arbitrary choice.

.. warning::

    If you add new forces to a :class:`System`, you must do so before creating the :class:`Simulation`.
    Once you create a :class:`Simulation`, modifying the :class:`System` will have no effect on that :class:`Simulation`.

Extracting and Reporting Forces (and other data)
************************************************

OpenMM provides reporters for three output formats: `PDB <https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html>`_, `PDBx/mmCIF <https://mmcif.wwpdb.org/pdbx-mmcif-home-page.html>`_ and `DCD <https://www.ks.uiuc.edu/Research/namd/2.6/ug/node13.html>`_.
All of those formats store only positions, not velocities, forces, or other data.  In this
section, we create a new reporter that outputs forces.  This illustrates two
important things: how to write a reporter, and how to query the simulation for
forces or other data.

Here is the definition of the :class:`ForceReporter` class:

.. samepage::
    ::

        class ForceReporter(object):
            def __init__(self, file, reportInterval):
                self._out = open(file, 'w')
                self._reportInterval = reportInterval

            def __del__(self):
                self._out.close()

            def describeNextReport(self, simulation):
                steps = self._reportInterval - simulation.currentStep%self._reportInterval
                return = {'steps': steps, 'periodic': None, 'include':['forces']}

            def report(self, simulation, state):
                forces = state.getForces().value_in_unit(kilojoules/mole/nanometer)
                for f in forces:
                    self._out.write('%g %g %g\n' % (f[0], f[1], f[2]))

    .. caption::

        :autonumber:`Example,ForceReporter`

The constructor and destructor are straightforward.  The arguments to the
constructor are the output filename and the interval (in time steps) at which it
should generate reports.  It opens the output file for writing and records the
reporting interval.  The destructor closes the file.

We then have two methods that every reporter must implement:
:meth:`describeNextReport()` and :meth:`report()`.  A Simulation object
periodically calls :meth:`describeNextReport()` on each of its reporters to
find out when that reporter will next generate a report, and what information
will be needed to generate it.  The return value should be a dictionary with the following content:

* :code:`steps` (int): The number of time steps until the next report.  We calculate this as
  *(report interval)*\ -\ *(current step)*\ %\ *(report interval)*\ .  For
  example, if we want a report every 100 steps and the simulation is currently on
  step 530, we will return 100-(530%100) = 70.
* :code:`include` (list of strings): The types of information that need to be included in the next report. The values in the list correspond to arguments to :meth:`Context.getState()`. Allowed values include :code:`'positions'`, :code:`'velocities'`, :code:`'forces'`, :code:`'energy'`, :code:`'parameters'`, :code:`'parameterDerivatives'`, and :code:`'integratorParameters'`.
* :code:`periodic` (bool, optional): Whether the positions should be wrapped to the periodic box.  If None or not set, it will
  automatically decide whether to wrap positions based on whether the System uses
  periodic boundary conditions.


When the time comes for the next scheduled report, the :class:`Simulation` calls
:meth:`report()` to generate the report.  The arguments are the :class:`Simulation`
object, and a :class:`State` that is guaranteed to contain all the information that was
requested by :meth:`describeNextReport()`\ .  A State object contains a
snapshot of information about the simulation, such as forces or particle
positions.  We call :meth:`getForces()` to retrieve the forces and convert
them to the units we want to output (kJ/mole/nm).  Then we loop over each value
and write it to the file.  To keep the example simple, we just print the values
in text format, one line per particle.  In a real program, you might choose a
different output format.

Now that we have defined this class, we can use it exactly like any other
reporter.  For example,
::

    simulation.reporters.append(ForceReporter('forces.txt', 100))

will output forces to a file called “forces.txt” every 100 time steps.

Computing Energies
******************

This example illustrates a different sort of analysis.  Instead of running a
simulation, assume we have already identified a set of structures we are
interested in.  These structures are saved in a set of PDB files.  We want to
loop over all the files in a directory, load them in one at a time, and compute
the potential energy of each one.  Assume we have already created our :class:`System` and
:class:`Simulation`.  The following lines perform the analysis:

.. samepage::
    ::

        import os
        for file in os.listdir('structures'):
            pdb = PDBFile(os.path.join('structures', file))
            simulation.context.setPositions(pdb.positions)
            state = simulation.context.getState(energy=True)
            print(file, state.getPotentialEnergy())

    .. caption::

        :autonumber:`Example,computing energies`

We use Python’s :code:`listdir()` function to list all the files in the
directory.  We create a :class:`PDBFile` object for each one and call
:meth:`setPositions()` on the Context to specify the particle positions loaded
from the PDB file.  We then compute the energy by calling :meth:`getState()`
with the option :code:`energy=True`\ , and print it to the console along
with the name of the file.


"""
replicaexchangesampler.py: Performs multistate sampling with replica exchange

This is part of the OpenMM molecular simulation toolkit.
See https://openmm.org/development.

Portions copyright (c) 2026 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import openmm as mm
import openmm.unit as unit
from openmm.app.internal.multistatesampler import MultistateSampler
import math
import random

class ReplicaExchangeSampler(object):
    """
    ReplicaExchangeSampler uses replica exchange to simulate a system in a collection of different states.  It supports
    both temperature replica exchange, in which the states correspond to different temperatures, and Hamiltonian replica
    exchange, in which they correspond to different potential functions.  It also can combine them to vary temperature
    and potential function at the same time.

    You provide a Simulation describing the system to simulate and a list of states.  It creates multiple replicas of the
    system, one for each state.  It then simulates all the replicas at once while periodically swapping states between
    them in a way that ensures correct sampling.

    Replica exchange is used for a variety of purposes.  Here are some of the more common ones.

    - Temperature replica exchange is most often used to accelerate sampling of rare events.  A particular transition might
      only happen rarely at physiological temperature.  You therefore simulate both the temperature you care about and also
      higher temperatures at which the transition happens more easily.  Allowing each replica to move between temperatures
      produces accurate sampling at low temperature of the entire phase space.

    - Hamiltonian replica exchange may be used to sample an alchemical transition between two endpoints.  It produces
      efficient sampling at every point along the transition pathway from which a free energy difference can be computed.

    Each state (not to be confused with a State object) is represented by a dict containing property values.  All states
    must specify values for the same properties.  They describe the ways in which the states differ from each other.  The
    following properties are supported.

    - 'temperature': the simulation temperature
    - 'groups': a set containing the force groups to include when computing the energy, for example {0, 2}.  It also may
      be a weighted sum specified as a dict.  For example, {0:1.0, 2:0.5} means to compute the energy of group 0 plus
      half the energy of group 2.
    - Context parameters, specified by the parameter name
    - Global variables defined by a CustomIntegrator, specified by the variable name

    For example, a typical use of temperature replica exchange might define the states as

    >>> states = [{'temperature':t*kelvin} for t in np.geomspace(300.0, 500.0, 10)]

    A typical use of Hamiltonian replica exchange might include forces that depend on a parameter `lambda`.  The states
    could be defined as

    >>> states = [{'lambda':x} for x in np.linspace(0.0, 1.0, 11)]

    States may also specify multiple properties.  This example includes all combinations of two different parameters.

    >>> states = [{'lambda1':x, 'lambda2':y} for x in np.linspace(0.0, 1.0, 6) for y in np.linspace(0.0, 1.0, 6)]

    To run a replica exchange simulation, first create a Simulation object in the normal way.  Then create a
    ReplicaExchangeSampler:

    >>> sampler = ReplicaExchangeSampler(states, simulation, stepsPerIteration)

    The third argument is the number of time steps to integrate during each iteration.  Then perform a specified number
    of sampling iterations by calling

    >>> sampler.simulate(iterations)

    Because replica exchange involves simulating many distinct replicas, the standard reporters used to generate output
    from a simulation are often not appropriate.  ReplicaExchangeSampler instead provides its own reporting mechanism.
    Define a function that takes a ReplicaExchangeSampler as its only argument:

    >>> def report(sampler):
    >>>     # Generate whatever output you want here

    and then add it to the sampler:

    >>> sampler.reporters.append(report)

    The function will be called after every iteration.

    Attributes
    ----------
    states: list[dict]
        The states to sample
    simulation: openmm.app.Simulation
        The Simulation defining the System, Context, and Integrator to use
    stepsPerIteration: int
        The number of time steps to integrate for each replica on each iteration
    reinitializeVelocities: bool
        If true, a replica's velocities are reinitialized from a Boltzmann distribution every time its state changes.
        This may sometimes improve stability, but also decreases efficiency.
    exchangesPerIteration: int
        The number of state exchange attempts between pairs of replicas to perform on each iteration.  This is initialized
        to len(states)**2.  There is rarely a reason to change it.
    currentIteration: int
        The number of iterations that have been completed
    reporters: list
        Reporting functions to invoke after each iteration
    replicaStateIndex: list[int]
        The current assignment of states to replicas.  replicaStateIndex[i] is the state of the i'th replica, specified
        as an index into states.
    replicaConformation: list[openmm.State]
        The current conformation of each replica, specified as a State object
    replicaStateEnergy: list[numpy.ndarray]
        The current potential energy of each replica in each state.  replicaStateEnergy[i][j] is the energy of replica
        i in state j.  Note that this includes the energy of each replica in every state, not just the state it is
        currently in, because all of them are required for performing exchanges.
    """
    def __init__(self, states: list[dict], simulation: "mm.app.Simulation", stepsPerIteration: int, reinitializeVelocities: bool = False):
        """
        Create a ReplicaExchangeSampler to sample a set of states.

        Parameters
        ----------
        states: list[dict]
            The states to sample
        simulation: openmm.app.Simulation
            The Simulation defining the System, Context, and Integrator to use
        stepsPerIteration: int
            The number of time steps to integrate for each replica on each iteration
        reinitializeVelocities: bool
            If true, a replica's velocities are reinitialized from a Boltzmann distribution every time its state changes.
            This may sometimes improve stability, but also decreases efficiency.
        """
        self.states = states
        self.simulation = simulation
        self.stepsPerIteration = stepsPerIteration
        self.reinitializeVelocities = reinitializeVelocities
        self.exchangesPerIteration = len(states)**2
        self.currentIteration = 0
        self.reporters = []
        self.replicaStateIndex = list(range(len(states)))
        self._previousReplicaStateIndex = self.replicaStateIndex[:]
        self.replicaConformation = [simulation.context.getState(positions=True, velocities=True, parameters=True, integratorParameters=True)]*len(states)
        self.replicaStateEnergy = None
        self._sampler = MultistateSampler(states, simulation.context)
        if 'temperature' in states[0]:
            temperature = [s['temperature'] for s in states]
            for i, t in enumerate(temperature):
                if not unit.is_quantity(t):
                    temperature[i] = t*unit.kelvin
            self._kT = [unit.MOLAR_GAS_CONSTANT_R*t for t in temperature]
        else:
            self._kT = None
            if not hasattr(simulation.integrator, 'getTemperature'):
                raise ValueError('Cannot determine temperature because the integrator does not have a getTemperature() method. '
                                 'Specify the temperature in each state dict.')

    def simulate(self, iterations: int):
        """
        Run the simulation to perform sampling.  Each iteration consists of simulating each replica for stepsPerIteration
        time steps, then exchanging states between replicas.

        Parameters
        ----------
        iterations: int
            the number of iterations to perform
        """
        for _ in range(iterations):
            self.replicaStateEnergy = None
            energies = []
            for i in range(len(self.states)):
                self.simulateReplica(i)
                energies.append(self._sampler.computeAllEnergies())
            self.replicaStateEnergy = energies
            self.exchangeReplicas()
            self.currentIteration += 1
            for reporter in self.reporters:
                reporter(self)

    def simulateReplica(self, index: int, steps: int | None = None):
        """
        Simulate a single replica for a specified number of steps.  Most often this is called by simulate(), but you also
        can call it directly, for example to equilibrate each replica before beginning sampling.

        You can create subclasses that override this method to change the sampling behavior.  For example, a subclass
        might perform Monte Carlo moves instead of or in addition to running dynamics.

        Parameters
        ----------
        index: int
            the index of the replica to simulate
        steps: int | None
            the number of time steps to simulate.  If None (the default), stepsPerIteration steps are simulated.
        """
        self.simulation.context.setState(self.replicaConformation[index])
        if self.reinitializeVelocities:
            # Reinitialize velocities for replicas whose states have changed.

            if self._previousReplicaStateIndex[index] != self.replicaStateIndex[index]:
                state = self.states[self.replicaStateIndex[index]]
                if 'temperature' in state:
                    temperature = state['temperature']
                else:
                    temperature = self.simulation.integrator.getTemperature()
                self.simulation.context.setVelocitiesToTemperature(temperature)
        elif self._kT is not None:
            # Rescale velocities for replicas whose temperatures have changed.

            kT1 = self._kT[self._previousReplicaStateIndex[index]]
            kT2 = self._kT[self.replicaStateIndex[index]]
            if kT1 != kT2:
                scale = math.sqrt(kT2/kT1)
                self.simulation.context.setVelocities(self.replicaConformation[index].getVelocities(asNumpy=True)*scale)
        self._previousReplicaStateIndex[index] = self.replicaStateIndex[index]
        self._sampler.applyState(self.replicaStateIndex[index])
        if steps is None:
            steps = self.stepsPerIteration
        self.simulation.integrator.step(steps)
        self.replicaConformation[index] = self.simulation.context.getState(positions=True, velocities=True, parameters=True, integratorParameters=True)

    def exchangeReplicas(self):
        """
        Attempt to exchange states between replicas.  This is called by simulate(), and there is not normally a reason
        to call it directly.  It performs exchangesPerIteration exchange attempts, each between two randomly chosen
        replicas.

        You can create subclasses that override this method to perform exchanges in different ways.
        """
        n = len(self.states)
        if n < 2:
            return
        for _ in range(self.exchangesPerIteration):
            # Select a random pair of replicas.

            while True:
                i = random.randrange(n)
                j = random.randrange(n)
                if i != j:
                    break
            si = self.replicaStateIndex[i]
            sj = self.replicaStateIndex[j]

            # Compute the probability of exchanging their states.

            if self._kT is None:
                kT = unit.MOLAR_GAS_CONSTANT_R*self.simulation.integrator.getTemperature()
                exponent = (self.replicaStateEnergy[i][si]-self.replicaStateEnergy[i][sj]+self.replicaStateEnergy[j][sj]-self.replicaStateEnergy[j][si])/kT
            else:
                exponent = (self.replicaStateEnergy[i][si]-self.replicaStateEnergy[j][si])/self._kT[si] + (self.replicaStateEnergy[j][sj]-self.replicaStateEnergy[i][sj])/self._kT[sj]
            prob = min(1.0, math.exp(exponent))
            if random.random() <= prob:
                # Accept the exchange.

                self.replicaStateIndex[i] = sj
                self.replicaStateIndex[j] = si

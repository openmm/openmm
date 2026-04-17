from __future__ import print_function

"""
expandedensemblesampler.py: Performs multistate sampling with the expanded ensemble method

This is part of the OpenMM molecular simulation toolkit.
See https://openmm.org/development.

Portions copyright (c) 2015-2026 Stanford University and the Authors.
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
__author__ = "Peter Eastman"
__version__ = "1.0"

import openmm.unit as unit
from openmm.app.internal.multistatesampler import MultistateSampler
import math
import random

class ExpandedEnsembleSampler(object):
    """
    ExpandedEnsembleSampler uses the expanded ensemble method to simulate a system in a collection of thermodynamic
    states.  It supports both the temperature version, in which the states correspond to different temperatures, and the
    Hamiltonian version, in which they correspond to different potential functions.  It also can combine them to vary
    temperature and potential function at the same time.

    You provide a Simulation describing the system to simulate and a list of states.  It simulates the system while
    periodically moving between states in a way that ensures correct sampling.  The method is based on Gibbs sampling,
    which alternates sampling of conformations x and thermodynamic states s.  It first performs molecular dynamics to
    sample P(x|s), the distribution of conformations for a fixed state.  It then chooses a new state from the
    distribution P(s|x) for the current conformation.  For more information on this method, see
    https://doi.org/10.33011/livecoms.4.1.1583.

    Expanded ensemble sampling is used for a variety of purposes.  Here are some of the more common ones.

    - The temperature version is most often used to accelerate sampling of rare events.  A particular transition might
      only happen rarely at physiological temperature.  You therefore simulate both the temperature you care about and
      also higher temperatures at which the transition happens more easily.  Allowing it to move between temperatures
      produces accurate sampling at low temperature of the entire phase space.

    - The Hamiltonian version may be used to sample an alchemical transition between two endpoints.  It produces
      efficient sampling at every point along the transition pathway from which a free energy difference can be computed.

    Each thermodynamic state (not to be confused with a State object) is represented by a dict containing property
    values.  All states must specify values for the same properties.  They describe the ways in which the states differ
    from each other.  The following properties are supported.

    - 'temperature': the simulation temperature
    - 'groups': a set containing the force groups to include when computing the energy, for example {0, 2}.  It also may
      be a weighted sum specified as a dict.  For example, {0:1.0, 2:0.5} means to compute the energy of group 0 plus
      half the energy of group 2.
    - Context parameters, specified by the parameter name
    - Global variables defined by a CustomIntegrator, specified by the variable name

    For example, a typical use of temperature sampling might define the states as

    >>> states = [{'temperature':t*kelvin} for t in np.geomspace(300.0, 500.0, 10)]

    A typical use of Hamiltonian sampling might include forces that depend on a parameter `lambda`.  The states could
    be defined as

    >>> states = [{'lambda':x} for x in np.linspace(0.0, 1.0, 11)]

    States may also specify multiple properties.  This example includes all combinations of two different parameters.

    >>> states = [{'lambda1':x, 'lambda2':y} for x in np.linspace(0.0, 1.0, 6) for y in np.linspace(0.0, 1.0, 6)]

    To run an expanded ensemble simulation, first create a Simulation object in the normal way.  Then create a
    ExpandedEnsembleSampler:

    >>> sampler = ExpandedEnsembleSampler(states, simulation, stepsPerIteration)

    The third argument is the number of time steps to integrate between attempted state changes.  Then perform a
    simulation by calling simulation.step() in the usual way.  The ExpandedEnsembleSampler adds a reporter to the
    Simulation that manages state changes as it runs.

    The probability distribution P(s|x) is often very nonuniform, which can lead to the simulation spending much more
    time in some states than others.  This is addressed by including a weight factor for each state when computing the
    probability.  Ideally the weights should be chosen so it spends equal time in every state.  You can specify the
    weights yourself with a constructor argument, but usually it is easier to let the ExpandedEnsembleSampler choose
    them automatically.  This is done with a combination of the Wang-Landau (https://doi.org/10.1103/PhysRevLett.86.2050)
    and Self-Adjusted Mixture Sampling (http://dx.doi.org/10.1080/10618600.2015.1113975) algorithms.  You should monitor
    the weights and discard the initial part of the simulation where they are changing rapidly.  Until the weights have
    converged, the simulation will not sample the correct distribution.

    To analyze the results of a simulation, it is essential to know what state it was in at every point in time.  An
    ExpandedEnsembleSampler can act as a reporter, generating output to a file that reports the current state and
    weights at regular intervals.  The reporting interval should generally be the same as the one used for writing a
    trajectory so the two will be synchronized with each other.

    Attributes
    ----------
    states: list[dict]
        The states to sample
    simulation: openmm.app.Simulation
        The Simulation defining the System, Context, and Integrator to use
    stepsPerIteration: int
        The number of time steps to integrate between attempted state changes
    reinitializeVelocities: bool
        If true, the simulation's velocities are reinitialized from a Boltzmann distribution every time its state
        changes.  This may sometimes improve stability, but also decreases efficiency.
    currentIteration: int
        The number of iterations that have been completed, each consisting of stepsPerIteration time steps followed by
        an attempted state change
    currentStateIndex: int
        The current state of the simulation, specified as an index into states.
    """

    def __init__(self, states: list[dict], simulation: "mm.app.Simulation", stepsPerIteration: int,
                 reinitializeVelocities: bool = False, weights: list[float] | None = None, reportInterval: int = 1000,
                 logFile: str | object | None = None):
        """Create a new ExpandedEnsembleSampler.

        Parameters
        ----------
        states: list[dict]
            The states to sample
        simulation: openmm.app.Simulation
            The Simulation defining the System, Context, and Integrator to use
        stepsPerIteration: int
            The number of time steps to integrate between attempted state changes
        reinitializeVelocities: bool
            If true, the simulation's velocities are reinitialized from a Boltzmann distribution every time its state
            changes.  This may sometimes improve stability, but also decreases efficiency.
        weights: list[float] | None
            The weights to use for each state.  If None, weights are chosen automatically as the simulation runs.
        reportInterval: int
            The frequency at which to write output, measured in time steps
        logFile: str | object | None
            An optional file to write a log to.  This may be either a file-like object or a string containing the path
            to the file.
        """
        self.states = states
        self.simulation = simulation
        self.stepsPerIteration = stepsPerIteration
        self.reinitializeVelocities = reinitializeVelocities
        self.currentIteration = 0
        self._stage = 0
        self.currentStateIndex = 0
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
        self.reportInterval = reportInterval

        # Initialize the weights.

        if weights is None:
            self._weights = [0.0]*len(states)
            self._updateWeights = True
            self._weightUpdateFactor = 1.0
            self._histogram = [0]*len(states)
            self._hasMadeTransition = False
        else:
            self._weights = weights
            self._updateWeights = False

        # Apply the current state to the simulation.

        self._sampler.applyState(self.currentStateIndex)

        # Add a reporter to the simulation which will handle the updates and reports.

        class EEReporter(object):
            def __init__(self, sampler):
                self.sampler = sampler

            def describeNextReport(self, simulation):
                sampler = self.sampler
                steps1 = sampler.stepsPerIteration - simulation.currentStep%sampler.stepsPerIteration
                steps2 = sampler.reportInterval - simulation.currentStep%sampler.reportInterval
                steps = min(steps1, steps2)
                isUpdateAttempt = (steps1 == steps)
                if isUpdateAttempt and not sampler.reinitializeVelocities:
                    return {'steps': steps, 'periodic':None, 'include':['velocities']}
                else:
                    return {'steps': steps, 'periodic':None, 'include':[]}

            def report(self, simulation, state):
                sampler = self.sampler
                if simulation.currentStep%sampler.stepsPerIteration == 0:
                    sampler.attemptStateChange(state)
                if simulation.currentStep%sampler.reportInterval == 0:
                    sampler._writeReport()

        self._openedFile = isinstance(logFile, str)
        if self._openedFile:
            self._out = open(logFile, 'w', 1)
        else:
            self._out = logFile
        simulation.reporters.append(EEReporter(self))

        # Write out the header line.

        if self._out is not None:
            headers = ['Steps', 'Iteration', 'State'] + [f'Weight {i}' for i in range(len(self.states))]
            print('"%s"' % (',').join(headers), file=self._out)

    def __del__(self):
        if self._openedFile:
            self._out.close()

    @property
    def weights(self):
        return [x-self._weights[0] for x in self._weights]

    def step(self, steps):
        """Advance the simulation by integrating a specified number of time steps."""
        self.simulation.step(steps)

    def attemptStateChange(self, state):
        """Attempt to move to a different state.  This is called automatically during the simulation, and there is not
        normally a reason to call it directly.

        You can create subclasses that override this method to select states in different ways."""

        # Compute the probability for each state.  This is done in log space to avoid overflow.

        self.currentIteration += 1
        energies = self._sampler.computeAllEnergies()
        if self._kT is None:
            kT = [unit.MOLAR_GAS_CONSTANT_R*self.simulation.integrator.getTemperature()]*len(self.states)
        else:
            kT = self._kT
        logProbability = [(self._weights[i]-energies[i]/kT[i]) for i in range(len(self._weights))]
        maxLogProb = max(logProbability)
        offset = maxLogProb + math.log(sum(math.exp(x-maxLogProb) for x in logProbability))
        probability = [math.exp(x-offset) for x in logProbability]

        # Select a new state.

        prevState = self.currentStateIndex
        r = random.random()
        for j in range(len(probability)):
            if r < probability[j]:
                if j != self.currentStateIndex:
                    # Select this state.

                    self._hasMadeTransition = True
                    self.currentStateIndex = j
                if self._updateWeights:
                    # Update the weights.

                    updateFactor = 1/self.currentIteration
                    if self._stage == 0:
                        if self._weightUpdateFactor >= updateFactor:
                            updateFactor = self._weightUpdateFactor
                        else:
                            self._stage = 1
                    for k in range(len(probability)):
                        self._weights[k] -= updateFactor*probability[k]
                    if self._stage == 0:
                        # Early in the simulation we reduce the update factor when the histogram is sufficiently flat.
                        # Later we will switch over to just using 1/iteration.

                        self._histogram[j] += 1
                        minCounts = min(self._histogram)
                        if minCounts > 20 and minCounts >= 0.2*sum(self._histogram)/len(self._histogram):
                            # Reduce the weight update factor and reset the histogram.

                            self._weightUpdateFactor *= 0.5
                            self._histogram = [0]*len(self.states)
                            self._weights = [x-self._weights[0] for x in self._weights]
                        elif not self._hasMadeTransition and probability[self.currentStateIndex] > 0.99 and sum(self._histogram) >= 10:
                            # Rapidly increase the weight update factor at the start of the simulation to find
                            # a reasonable starting value.

                            self._weightUpdateFactor *= 2.0
                            self._histogram = [0]*len(self.states)
                        if self._weightUpdateFactor < 1/self.currentIteration:
                            self.stage = 1
                break
            r -= probability[j]

        # Apply the state.

        self._sampler.applyState(self.currentStateIndex)

        # Reinitialize or rescale the velocities.

        if self.reinitializeVelocities:
            if 'temperature' in self.states[self.currentStateIndex]:
                temperature = self.states[self.currentStateIndex]['temperature']
            else:
                temperature = self.simulation.integrator.getTemperature()
            self.simulation.context.setVelocitiesToTemperature(temperature)
        elif kT[prevState] != kT[self.currentStateIndex]:
            scale = math.sqrt(kT[self.currentStateIndex]/kT[prevState])
            self.simulation.context.setVelocities(scale*state.getVelocities(asNumpy=True))

    def _writeReport(self):
        """Write out a line to the report."""
        if self._out is not None:
            values = [self.currentStateIndex]+self.weights
            print(f'{self.simulation.currentStep},{self.currentIteration},' + ','.join('%g' % v for v in values), file=self._out)

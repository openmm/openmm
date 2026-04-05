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
    def __init__(self, states: list, simulation: "mm.app.Simulation", stepsPerIteration: int):
        self.states = states
        self.simulation = simulation
        self.stepsPerIteration = stepsPerIteration
        self.exchangesPerIteration = len(states)**3
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
        self.simulation.context.setState(self.replicaConformation[index])
        if self._kT is not None:
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
        n = len(self.states)
        if n < 2:
            return
        for _ in range(self.exchangesPerIteration):
            while True:
                i = random.randrange(n)
                j = random.randrange(n)
                if i != j:
                    break
            si = self.replicaStateIndex[i]
            sj = self.replicaStateIndex[j]
            if self._kT is None:
                kT = unit.MOLAR_GAS_CONSTANT_R*self.simulation.integrator.getTemperature()
                exponent = (self.replicaStateEnergy[i][si]-self.replicaStateEnergy[i][sj]+self.replicaStateEnergy[j][sj]-self.replicaStateEnergy[j][si])/kT
            else:
                exponent = (self.replicaStateEnergy[i][si]-self.replicaStateEnergy[j][si])/self._kT[si] + (self.replicaStateEnergy[j][sj]-self.replicaStateEnergy[i][sj])/self._kT[sj]
            prob = min(1.0, math.exp(exponent))
            if random.random() <= prob:
                self.replicaStateIndex[i] = sj
                self.replicaStateIndex[j] = si

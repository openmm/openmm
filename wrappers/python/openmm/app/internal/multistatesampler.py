"""
multistatesampler.py: Utilities for multistate sampling algorithms.

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

import openmm
from openmm.unit import kilojoules_per_mole

"""This is a utility for internal use by algorithms that perform multistate sampling.  Given the definitions of a set of
states, it can put a context into any state and efficiently evaluate the potential energies of states.

Each (not to be confused with a State object) is represented by a dict containing property values.  All states must
specify values for the same properties.  The following properties are supported.

- 'temperature': the simulation temperature.  It knows how to set the temperature for all standard integrators and forces.
- 'groups': a set containing the force groups to include when computing the energy, for example {0, 2}.  It also may be
  a weighted sum specified as a dict.  For example, {0:1.0, 2:0.5} means to return the energy of group 0 plus half the
  energy of group 2.
- Context parameters
- Global variables defined by a CustomIntegrator
"""
class MultistateSampler(object):
    def __init__(self, states, context):
        """Create a MultistateSampler.

        Parameters
        ----------
        states: list
            the states to sample.  Each entry should be a dict containing the property values for one state.
        context: Context
            the Context to apply the states to
        """
        self.states = states
        self.context = context
        keys = set(states[0].keys())
        for s in states:
            if set(s.keys()) != keys:
                raise ValueError('All states must include the same set of properties')

        # Process the keys and determine how to set the simulation to a state.

        self.parameters = [{} for _ in states]
        energy_parameters = [{} for _ in states]
        self.variables = [{} for _ in states]
        self.temperature = None
        self.groups = None
        if isinstance(context.getIntegrator(), openmm.CustomIntegrator):
            variables = set(context.getIntegrator().getGlobalVariableName(i) for i in range(context.getIntegrator().getNumGlobalVariables()))
        else:
            variables = set()
        for key in keys:
            if key == 'temperature':
                # Setting the temperature can involve call setTemperature() on an integrator and setting context parameters.

                if hasattr(context.getIntegrator(), 'setTemperature'):
                    self.temperature = [s['temperature'] for s in states]
                for force in context.getSystem().getForces():
                    if hasattr(type(force), 'Temperature'):
                        param = force.Temperature()
                        if param in context.getParameters():
                            for i, s in enumerate(states):
                                self.parameters[i][param] = s['temperature']
                if 'temperature' in variables:
                    for i, s in enumerate(states):
                        self.variables[i]['temperature'] = s['temperature']
            elif key == 'groups':
                # Force groups can be specified with either a set or a dict.

                self.groups = []
                for s in states:
                    groups = s['groups']
                    if isinstance(groups, set):
                        groups = {g:1.0 for g in groups}
                    self.groups.append(groups)
            elif key in context.getParameters():
                # A context parameter.

                for i, s in enumerate(states):
                    self.parameters[i][key] = s[key]
                    energy_parameters[i][key] = s[key]
            elif key in variables:
                # A global variable of a CustomIntegrator.

                for i, s in enumerate(states):
                    self.variables[i][key] = s[key]
            else:
                raise ValueError(f'Unknown property "{key}"')

        # We now identify subsets of states whose energies can be evaluated more efficiently.  Two states are in the
        # same subset if they have identical context parameters (ignoring parameters used only to set temperature).

        self.subsets = []
        for i, s in enumerate(states):
            match = None
            for subset in self.subsets:
                if energy_parameters[i] == energy_parameters[subset[0]]:
                    match = subset
                    break
            if match is None:
                self.subsets.append([i])
            else:
                match.append(i)

        # States within a subset may still vary in what force groups they include.  Figure out the most efficient way
        # of evaluating all states in a subset.

        if self.groups is not None:
            self.groups_of_groups = [[] for _ in range(len(self.subsets))]
            for i, subset in enumerate(self.subsets):
                # Find all force groups that appear in this subset.

                all_groups = set()
                for j in subset:
                    for k in self.groups[j].keys():
                        all_groups.add(k)

                # If two force groups always appear together with the same weight, they can be computed together.

                while len(all_groups) > 0:
                    first = next(iter(all_groups))
                    group_of_groups = {first}
                    for group in list(all_groups)[1:]:
                        match = True
                        for j in subset:
                            if self.groups[j].get(first) != self.groups[j].get(group):
                                match = False
                                break
                        if match:
                            group_of_groups.add(group)
                            all_groups.remove(group)
                    self.groups_of_groups[i].append(group_of_groups)
                    all_groups.remove(first)

    def applyState(self, index):
        """Modify the Context to match one of the states.

        Parameters
        ----------
        index: int
            the index of the state in self.states
        """
        if self.temperature is not None:
            self.context.getIntegrator().setTemperature(self.temperature[index])
        for name, value in self.parameters[index].items():
            self.context.setParameter(name, value)
        for name, value in self.variables[index].items():
            self.context.getIntegrator().setGlobalVariableByName(name, value)

    def computeEnergy(self, index):
        """Compute the potential energy of one of the states.  This method calls applyState(), so when it returns, the
        Context will be in the specified state.

        If you need to compute the energies of all states, computeAllEnergies() is much more efficient than calling this
        method for each one.

        Parameters
        ----------
        index: int
            the index of the state in self.states

        Returns
        -------
        the potential energy of the state
        """
        self.applyState(index)
        if self.groups is None:
            return self.context.getState(energy=True).getPotentialEnergy()
        weights = set(self.groups[index].values())
        energy = 0*kilojoules_per_mole
        for weight in weights:
            groups = set(g for g, w in self.groups[index].items() if w == weight)
            energy += weight*self.context.getState(energy=True, groups=groups).getPotentialEnergy()
        return energy

    def computeAllEnergies(self):
        """Compute the potential energies of states.  This is much more efficient than calling computeEnergy() for each
        state individually.

        If you only care about energy differences between states, not absolute energies, use computeRelativeEnergies()
        instead.  It can be even more efficient.

        When this method returns, it is undefined which state the Context will be in.

        Returns
        -------
        an array containing the potential energies of all states in the order they appear in self.states
        """
        energies = [0*kilojoules_per_mole for _ in self.states]
        for i, subset in enumerate(self.subsets):
            if self.groups is None:
                # States don't depend on force groups, so we can just evaluate the energy of one state.

                energy = self.computeEnergy(subset[0])
                for j in subset:
                    energies[j] = energy
            else:
                # Compute the energy of each set of groups, and add them with the proper weights for each state.

                self.applyState(subset[0])
                for groups in self.groups_of_groups[i]:
                    energy = self.context.getState(energy=True, groups=groups).getPotentialEnergy()
                    for j in subset:
                        g = next(iter(groups))
                        if g in self.groups[j]:
                            energies[j] += energy*self.groups[j][g]
        return energies

    def computeRelativeEnergies(self):
        """This is similar to computeAllEnergies(), but the energies are shifted by a constant to make the energy of the
        first state exactly zero.  The advantage is that if states differ only in temperature or other ways that do not
        affect energy, there is no need to compute any energies at all.

        When this method returns, it is undefined which state the Context will be in.

        Returns
        -------
        an array containing the shifted potential energies of all states in the order they appear in self.states
        """
        if len(self.subsets) == 1 and self.groups is None:
            return [0*kilojoules_per_mole]*len(self.states)
        energies = self.computeAllEnergies()
        return [e-energies[0] for e in energies]
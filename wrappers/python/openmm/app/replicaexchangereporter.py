"""
replicaexchangereporter.py: Writes output for replica exchange simulations

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
import openmm.app as app
import openmm.unit as unit
from openmm.app.internal import safesave
import os

class ReplicaExchangeReporter(object):
    """
    This is a reporter for use with ReplicaExchangeSampler.  It outputs the most important information required for
    nearly all replica exchange simulations.  It must be added to the sampler's list of reporters:

    >>> sampler.reporters.append(ReplicaExchangeReporter(directory, interval, sampler))

    As the simulation runs, it creates the following files in a user specified directory.

    - A CSV file containing a log of what state each replica was in at each iteration.
    - A set of XML files containing serialized State objects for each replica. These serve as checkpoints, allowing a
      simulation to be resumed later.
    - (optional) A trajectory file for each thermodynamic state.  The coordinates in these files jump discontinuously
      whenever an exchange happens.
    - (optional) A trajectory file for each replica.  The coordinates in these files are continuous, but the
      thermodynamic states they explore change with time.
    - (optional) A CSV file containing the reduced energy (potential energy divided by kT) of every replica in every
      state at each iteration.  This information is useful for calculating free energies.

    The reduced energies are written as a matrix in flattened order for each iteration.  You can load the file with
    NumPy using the command

    >>> energy = np.loadtxt('energy.csv', delimiter=',').reshape(-1, n_states, n_states)

    `energy[i][j][k]` is the reduced energy of state k for replica j in iteration i.

    To resume a simulation from the saved files, pass `resume=True` to the constructor.  It will load all necessary
    information and configure the ReplicaExchangeSampler correctly.
    """
    def __init__(self, directory: str, reportInterval: int, sampler: "app.ReplicaExchangeSampler",
                 trajectoryPerState: bool = True, trajectoryPerReplica: bool = False, trajectoryFormat: str = 'xtc',
                 enforcePeriodicBox: bool | None = None, energy: bool = False, resume: bool = False):
        """
        Create a ReplicaExchangeReporter.

        Parameters
        ----------
        directory: str
            The path to the directory where the files will be written.  The directory will be created if it does not
            already exist.  The directory must be empty unless resume is True, in which case it must contain the files
            from an earlier simulation.
        reportInterval: int
            The frequency at which to write output, measured in iterations
        sampler: openmm.app.ReplicaExchangeSampler
            The sampler this reporter will be added to
        trajectoryPerState: bool
            If True, a trajectory file will be saved for each thermodynamic state.
        trajectoryPerReplica: bool
            If True, a trajectory file will be saved for each replica.
        trajectoryFormat: str
            The format in which to save trajectories.  Supported options are 'xtc' and 'dcd'.
        enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
        energy: bool
            If True, a CSV file will be saved containing reduced energies.
        resume: bool
            Specifies whether to resume an earlier simulation.  If True, the checkpoint and log data will be loaded into
            the ReplicaExchangeSampler, and future output will be appended to the existing files.
        """
        trajectoryFormat = trajectoryFormat.lower()
        self.directory = directory
        self.reportInterval = reportInterval
        self.trajectoryPerState = trajectoryPerState
        self.trajectoryPerReplica = trajectoryPerReplica
        self.format = format
        self._log = None
        self._energy = None
        numStates = len(sampler.states)

        # Validate the inputs and create the output directory if necessary.

        if trajectoryFormat not in ('xtc', 'dcd'):
            raise ValueError(f'Unsupported trajectory format: {trajectoryFormat}.  Allowed values are "xtc" and "dcd".')
        if resume:
            if not os.path.isdir(directory):
                raise ValueError(f'Cannot resume because the directory does not exist: {directory}')

            def checkExists(filename):
                if not os.path.isfile(os.path.join(directory, filename)):
                    raise ValueError(f'Cannot resume because the file {filename} does not exist.')

            checkExists('log.csv')
            if energy:
                checkExists('energy.csv')
            for i in range(numStates):
                checkExists(f'checkpoint_{i}.xml')
                if trajectoryPerState:
                    checkExists(f'state_{i}.{trajectoryFormat}')
                if trajectoryPerReplica:
                    checkExists(f'replica_{i}.{trajectoryFormat}')
        else:
            if os.path.isdir(directory):
                if len(os.listdir(directory)) > 0:
                    raise ValueError(f'The directory {directory} is not empty.')
            else:
                os.makedirs(directory)

        # Load state assignments and checkpoints if resuming another simulation.

        if resume:
            for i in range(numStates):
                with open(os.path.join(directory, f'checkpoint_{i}.xml')) as input:
                    sampler.replicaConformation[i] = mm.XmlSerializer.deserialize(input.read())
            log = None
            with open(os.path.join(directory, 'log.csv')) as input:
                for line in input:
                    log = line
            if log is None:
                raise ValueError('Cannot resume because log file is empty.')
            fields = log.split(',')
            sampler.replicaStateIndex = [int(x) for x in fields[2:(numStates+2)]]
            sampler._previousReplicaStateIndex = sampler.replicaStateIndex[:]
            sampler.currentIteration = int(fields[0])

        # Create reporters and open files for writing output.

        def createReporter(filename):
            path = os.path.join(directory, filename)
            interval = reportInterval*sampler.stepsPerIteration
            if trajectoryFormat == 'xtc':
                return app.XTCReporter(path, reportInterval*sampler.stepsPerIteration, append=resume, enforcePeriodicBox=enforcePeriodicBox)
            return app.DCDReporter(path, interval, append=resume, enforcePeriodicBox=enforcePeriodicBox)

        self._stateReporters = []
        self._replicaReporters = []
        for i in range(numStates):
            if trajectoryPerState:
                self._stateReporters.append(createReporter(f'state_{i}.{trajectoryFormat}'))
            if trajectoryPerReplica:
                self._replicaReporters.append(createReporter(f'replica_{i}.{trajectoryFormat}'))
        self._log = open(os.path.join(directory, 'log.csv'), 'a' if resume else 'w')
        if not resume:
            print(','.join(['Iteration', 'Step']+[f'Replica_{i}_State' for i in range(numStates)]), file=self._log)
            self._log.flush()
        if energy:
            self._energy = open(os.path.join(directory, 'energy.csv'), 'a' if resume else 'w')

    def __del__(self):
        if self._log is not None:
            self._log.close()
        if self._energy is not None:
            self._energy.close()

    def __call__(self, sampler: "app.ReplicaExchangeSampler"):
        """This is invoked by the ReplicaExchangeSampler at the end of every iteration to generate output."""
        if sampler.currentIteration % self.reportInterval != 0:
            return
        logData = [sampler.currentIteration, sampler.simulation.currentStep]+sampler.replicaStateIndex
        print(','.join([str(x) for x in logData]), file=self._log)
        self._log.flush()
        if self._energy is not None:
            import numpy as np
            energy = np.array(sampler.replicaStateEnergy)
            if sampler._kT is None:
                energy /= unit.MOLAR_GAS_CONSTANT_R*sampler.simulation.integrator.getTemperature()
            else:
                energy /= sampler._kT
            print(','.join([str(x) for x in energy.flatten()]), file=self._energy)
            self._energy.flush()
        for i, conf in enumerate(sampler.replicaConformation):
            if len(self._stateReporters) > 0:
                self._stateReporters[sampler.replicaStateIndex[i]].report(sampler.simulation, conf)
            if len(self._replicaReporters) > 0:
                self._replicaReporters[i].report(sampler.simulation, conf)
            safesave.save(mm.XmlSerializer.serialize(conf), os.path.join(self.directory, f'checkpoint_{i}.xml'))

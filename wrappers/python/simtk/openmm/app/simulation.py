"""
simulation.py: Provides a simplified API for running simulations.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2015 Stanford University and the Authors.
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
from __future__ import absolute_import
__author__ = "Peter Eastman"
__version__ = "1.0"

import simtk.openmm as mm
import simtk.unit as unit
import sys
from datetime import datetime, timedelta
try:
    string_types = (unicode, str)
except NameError:
    string_types = (str,)

class Simulation(object):
    """Simulation provides a simplified API for running simulations with OpenMM and reporting results.

    A Simulation ties together various objects used for running a simulation: a Topology, System,
    Integrator, and Context.  To use it, you provide the Topology, System, and Integrator, and it
    creates the Context automatically.

    Simulation also maintains a list of "reporter" objects that record or analyze data as the simulation
    runs, such as writing coordinates to files or displaying structures on the screen.  For example,
    the following line will cause a file called "output.pdb" to be created, and a structure written to
    it every 1000 time steps:

    simulation.reporters.append(PDBReporter('output.pdb', 1000))
    """

    def __init__(self, topology, system, integrator, platform=None, platformProperties=None, state=None):
        """Create a Simulation.

        Parameters
        ----------
        topology : Topology
            A Topology describing the the system to simulate
        system : System or XML file name
            The OpenMM System object to simulate (or the name of an XML file
            with a serialized System)
        integrator : Integrator or XML file name
            The OpenMM Integrator to use for simulating the System (or the name
            of an XML file with a serialized System)
        platform : Platform=None
            If not None, the OpenMM Platform to use
        platformProperties : map=None
            If not None, a set of platform-specific properties to pass to the
            Context's constructor
        state : XML file name=None
            The name of an XML file containing a serialized State. If not None,
            the information stored in state will be transferred to the generated
            Simulation object.
        """
        self.topology = topology
        ## The System being simulated
        if isinstance(system, string_types):
            with open(system, 'r') as f:
                self.system = mm.XmlSerializer.deserialize(f.read())
        else:
            self.system = system
        ## The Integrator used to advance the simulation
        if isinstance(integrator, string_types):
            with open(integrator, 'r') as f:
                self.integrator = mm.XmlSerializer.deserialize(f.read())
        else:
            self.integrator = integrator
        ## The index of the current time step
        self.currentStep = 0
        ## A list of reporters to invoke during the simulation
        self.reporters = []
        if platform is None:
            ## The Context containing the current state of the simulation
            self.context = mm.Context(self.system, self.integrator)
        elif platformProperties is None:
            self.context = mm.Context(self.system, self.integrator, platform)
        else:
            self.context = mm.Context(self.system, self.integrator, platform, platformProperties)
        if state is not None:
            with open(state, 'r') as f:
                self.context.setState(mm.XmlSerializer.deserialize(f.read()))
        ## Determines whether or not we are using PBC. Try from the System first,
        ## fall back to Topology if that doesn't work
        try:
            self._usesPBC = self.system.usesPeriodicBoundaryConditions()
        except Exception: # OpenMM just raises Exception if it's not implemented everywhere
            self._usesPBC = topology.getUnitCellDimensions() is not None

    def minimizeEnergy(self, tolerance=10*unit.kilojoule/unit.mole, maxIterations=0):
        """Perform a local energy minimization on the system.

        Parameters
        ----------
        tolerance : energy=10*kilojoules/mole
            The energy tolerance to which the system should be minimized
        maxIterations : int=0
            The maximum number of iterations to perform.  If this is 0,
            minimization is continued until the results converge without regard
            to how many iterations it takes.
        """
        mm.LocalEnergyMinimizer.minimize(self.context, tolerance, maxIterations)

    def step(self, steps):
        """Advance the simulation by integrating a specified number of time steps."""
        self._simulate(endStep=self.currentStep+steps)

    def runForClockTime(self, time, checkpointFile=None, stateFile=None, checkpointInterval=None):
        """Advance the simulation by integrating time steps until a fixed amount of clock time has elapsed.

        This is useful when you have a limited amount of computer time available, and want to run the longest simulation
        possible in that time.  This method will continue taking time steps until the specified clock time has elapsed,
        then return.  It also can automatically write out a checkpoint and/or state file before returning, so you can
        later resume the simulation.  Another option allows it to write checkpoints or states at regular intervals, so
        you can resume even if the simulation is interrupted before the time limit is reached.

        Parameters
        ----------
        time : time
            the amount of time to run for.  If no units are specified, it is
            assumed to be a number of hours.
        checkpointFile : string or file=None
            if specified, a checkpoint file will be written at the end of the
            simulation (and optionally at regular intervals before then) by
            passing this to saveCheckpoint().
        stateFile : string or file=None
            if specified, a state file will be written at the end of the
            simulation (and optionally at regular intervals before then) by
            passing this to saveState().
        checkpointInterval : time=None
            if specified, checkpoints and/or states will be written at regular
            intervals during the simulation, in addition to writing a final
            version at the end.  If no units are specified, this is assumed to
            be in hours.
        """
        if unit.is_quantity(time):
            time = time.value_in_unit(unit.hours)
        if unit.is_quantity(checkpointInterval):
            checkpointInterval = checkpointInterval.value_in_unit(unit.hours)
        endTime = datetime.now()+timedelta(hours=time)
        while (datetime.now() < endTime):
            if checkpointInterval is None:
                nextTime = endTime
            else:
                nextTime = datetime.now()+timedelta(hours=checkpointInterval)
                if nextTime > endTime:
                    nextTime = endTime
            self._simulate(endTime=nextTime)
            if checkpointFile is not None:
                self.saveCheckpoint(checkpointFile)
            if stateFile is not None:
                self.saveState(stateFile)

    def _simulate(self, endStep=None, endTime=None):
        if endStep is None:
            endStep = sys.maxsize
        nextReport = [None]*len(self.reporters)
        while self.currentStep < endStep and (endTime is None or datetime.now() < endTime):
            nextSteps = endStep-self.currentStep
            
            # Find when the next report will happen.
            
            anyReport = False
            for i, reporter in enumerate(self.reporters):
                nextReport[i] = reporter.describeNextReport(self)
                if nextReport[i][0] > 0 and nextReport[i][0] <= nextSteps:
                    nextSteps = nextReport[i][0]
                    anyReport = True
            stepsToGo = nextSteps
            while stepsToGo > 10:
                self.integrator.step(10) # Only take 10 steps at a time, to give Python more chances to respond to a control-c.
                self.currentStep += 10
                stepsToGo -= 10
                if endTime is not None and datetime.now() >= endTime:
                    break
            else:
                self.integrator.step(stepsToGo)
                self.currentStep += stepsToGo
                stepsToGo = 0
            if anyReport and stepsToGo == 0:
                # One or more reporters are ready to generate reports.  Organize them into three
                # groups: ones that want wrapped positions, ones that want unwrapped positions,
                # and ones that don't care about positions.
                
                wrapped = []
                unwrapped = []
                either = []
                for reporter, report in zip(self.reporters, nextReport):
                    if report[0] == nextSteps:
                        if len(report) > 5:
                            wantWrap = report[5]
                            if wantWrap is None:
                                wantWrap = self._usesPBC
                        else:
                            wantWrap = self._usesPBC
                        if not report[1]:
                            either.append((reporter, report))
                        elif wantWrap:
                            wrapped.append((reporter, report))
                        else:
                            unwrapped.append((reporter, report))
                if len(wrapped) > len(unwrapped):
                    wrapped += either
                else:
                    unwrapped += either
                
                # Generate the reports.

                if len(wrapped) > 0:
                    self._generate_reports(wrapped, True)
                if len(unwrapped) > 0:
                    self._generate_reports(unwrapped, False)
    
    def _generate_reports(self, reports, periodic):
        getPositions = False
        getVelocities = False
        getForces = False
        getEnergy = False
        for reporter, next in reports:
            if next[1]:
                getPositions = True
            if next[2]:
                getVelocities = True
            if next[3]:
                getForces = True
            if next[4]:
                getEnergy = True
        state = self.context.getState(getPositions=getPositions, getVelocities=getVelocities, getForces=getForces,
                                      getEnergy=getEnergy, getParameters=True, enforcePeriodicBox=periodic)
        for reporter, next in reports:
            reporter.report(self, state)

    def saveCheckpoint(self, file):
        """Save a checkpoint of the simulation to a file.

        The output is a binary file that contains a complete representation of the current state of the Simulation.
        It includes both publicly visible data such as the particle positions and velocities, and also internal data
        such as the states of random number generators.  Reloading the checkpoint will put the Simulation back into
        precisely the same state it had before, so it can be exactly continued.

        A checkpoint file is highly specific to the Simulation it was created from.  It can only be loaded into
        another Simulation that has an identical System, uses the same Platform and OpenMM version, and is running on
        identical hardware.  If you need a more portable way to resume simulations, consider using saveState() instead.

        Parameters
        ----------
        file : string or file
            a File-like object to write the checkpoint to, or alternatively a
            filename
        """
        if isinstance(file, str):
            with open(file, 'wb') as f:
                f.write(self.context.createCheckpoint())
        else:
            file.write(self.context.createCheckpoint())

    def loadCheckpoint(self, file):
        """Load a checkpoint file that was created with saveCheckpoint().

        Parameters
        ----------
        file : string or file
            a File-like object to load the checkpoint from, or alternatively a
            filename
        """
        if isinstance(file, str):
            with open(file, 'rb') as f:
                self.context.loadCheckpoint(f.read())
        else:
            self.context.loadCheckpoint(file.read())

    def saveState(self, file):
        """Save the current state of the simulation to a file.

        The output is an XML file containing a serialized State object.  It includes all publicly visible data,
        including positions, velocities, and parameters.  Reloading the State will put the Simulation back into
        approximately the same state it had before.

        Unlike saveCheckpoint(), this does not store internal data such as the states of random number generators.
        Therefore, you should not expect the following trajectory to be identical to what would have been produced
        with the original Simulation.  On the other hand, this means it is portable across different Platforms or
        hardware.

        Parameters
        ----------
        file : string or file
            a File-like object to write the state to, or alternatively a
            filename
        """
        state = self.context.getState(getPositions=True, getVelocities=True, getParameters=True)
        xml = mm.XmlSerializer.serialize(state)
        if isinstance(file, str):
            with open(file, 'w') as f:
                f.write(xml)
        else:
            file.write(xml)

    def loadState(self, file):
        """Load a State file that was created with saveState().

        Parameters
        ----------
        file : string or file
            a File-like object to load the state from, or alternatively a
            filename
        """
        if isinstance(file, str):
            with open(file, 'r') as f:
                xml = f.read()
        else:
            xml = file.read()
        self.context.setState(mm.XmlSerializer.deserialize(xml))

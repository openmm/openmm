"""
progressreporter: Human-readable output about the progress of a simulation

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2013 Stanford University and the Authors.
Authors: Robert McGibbon
Contributors: Lee-Ping Wang

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

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import math
import time

from simtk import unit
import simtk.openmm as mm
from simtk.openmm.app import StateDataReporter

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

def prettyTime(secs):
    """Format the time in a pretty way"""

    if math.isnan(secs):
        return "??"
    secs = int(secs)

    days = secs // 86400
    secs -= 86400*days

    hrs = secs // 3600
    secs -= 3600*hrs

    mins = secs // 60
    secs -= 60*mins

    if days > 0:
        s = "%d:%d:%02d:%02d" % (days, hrs, mins, secs)
    elif hrs > 0:
        s = "%d:%02d:%02d" % (hrs, mins, secs)
    elif mins > 0:
        s = "%d:%02d" % (mins, secs)
    else:
        s = "0:%02d" % secs

    return s

class ProgressReporter(StateDataReporter):
    """ProgressReporter outputs information about the progress of a simulation
    in a human readable, ASCII format

    To use it, create a ProgressReporter, then add it to the Simulation's list
    of reporters.
    """
    def __init__(self, reportInterval, totalSteps, file=sys.stdout):
        """Create a ProgressReporter

        Parameters:
         - reportInterval (int) The interval (in time steps) at which to write frames
         - totalSteps (int) The total number of steps in the simulations, for
         calculating the estimated time it will finish.
         - file (string or file) The file to write to, specified as a file name or file object
        """
        
        super(ProgressReporter, self).__init__(file, reportInterval, step=False, time=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
            temperature=True)

        self._totalSteps = totalSteps

    def _initializeConstants(self, simulation, state):
        if simulation.topology.getUnitCellDimensions() is not None:
            self._volume = True
            self._density = True

        platform = mm.context.getPlatform()
        print >>self._out, 'Number of available platforms: %d' % mm.Platform.getNumPlatforms()
        print >>self._out, 'Selected Platform: %s' % platform.getName()
        for key in platform.getPropertyNames():
            print >>self._out, '%s = %s', key, platform.getPropertyValue(simulation.context, key)

        # this needs to come after _density and _volume are set so
        # that the mass gets computed, if needed
        super(ProgressReporter, self)._initializeConstants(simulation)

        # initialize these as late as possible, so that as little initialization
        # code gets counted in the elapsed walltime. When that happens, it
        # makes it look like the simulation is getting faster and faster,
        # since that time is being amoritized out.
        self._initialWallTime = time.time()
        self._initialStep = simulation.currentStep
        self._initialSimTime = state.getTime()

    def _constructReportValues(self, simulation, state):
        progressPercent = 100 * float(simulation.currentStep - self._initialStep) / self._totalSteps
        if progressPercent > 0:
            timeLeft = (time.time() - self._initialWallTime) * (100.0 - progressPercent) / progressPercent

            elapsedSim = (state.getTime() - self._initialSimTime).value_in_unit(unit.nanoseconds)
            walltime = ((time.time() - self._initialWallTime)*unit.seconds).value_in_unit(unit.days)
            rate = elapsedSim / walltime
        else:
            timeLeft = float('nan')
            rate = 0

        values = [progressPercent, prettyTime(timeLeft), rate] + \
                 super(ProgressReporter, self)._constructReportValues(simulation, state)
        return values

    def _constructHeaders(self):
        headers = [('Progress', '(%)'),
                   ('WallTime Left', '(d:h:m:s)'),
                   ('Speed', '(ns/day)'),
                   ('Time',  '(ps)'),
                   ('P.E.', '(kJ/mol)'),
                   ('K.E.', '(kJ/mol)'),
                   ('Total E.', '(kJ/mol)'),
                   ('Temp', '(K)'),
                  ]

        widths =  [8,          15,      10,        13,       15,
                   15,         15,      13]
        formats = ['%7.3f%%', '%15s', '%10.2f', '%13.5f', '%15.5f',
                   '%15.5f', '%15.5f', '%13.5f']

        if self._volume:
            headers.append(('Vol', '(nm^3)'))
            formats.append('%10.4f')
            widths.append(10)
        if self._density:
            headers.append(('Rho', '(g/mL)'))
            formats.append('%10.4f')
            widths.append(10)

        self._formats = formats

        row1, row2 = zip(*headers)
        headerwidths = ['%{w}s'.format(w=w) for w in widths]
        print >>self._out, ' '.join(f % e for f, e in zip(headerwidths, row1))
        print >>self._out, ' '.join(f % e for f, e in zip(headerwidths, row2))


    def report(self, simulation, state):
        """Generate a report.

        Parameters:
         - simulation (Simulation) The Simulation to generate a report for
         - state (State) The current state of the simulation
        """
        if not self._hasInitialized:
            self._initializeConstants(simulation, state)
            self._constructHeaders()
            self._hasInitialized = True

        # Check for errors.
        self._checkForErrors(simulation, state)

        # Query for the values
        values = self._constructReportValues(simulation, state)

        print >>self._out, ' '.join(f % v for f, v in zip(self._formats, values))

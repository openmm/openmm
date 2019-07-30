from __future__ import print_function

"""
simulatedtempering.py: Implements simulated tempering

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2015 Stanford University and the Authors.
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

import simtk.unit as unit
import math
import random
from sys import stdout

try:
    import bz2
    have_bz2 = True
except: have_bz2 = False

try:
    import gzip
    have_gzip = True
except: have_gzip = False

class SimulatedTempering(object):
    """SimulatedTempering implements the simulated tempering algorithm for accelerated sampling.
    
    It runs a simulation while allowing the temperature to vary.  At high temperatures, it can more easily cross
    energy barriers to explore a wider area of conformation space.  At low temperatures, it can thoroughly
    explore each local region.  For details, see Marinari, E. and Parisi, G., Europhys. Lett. 19(6). pp. 451-458 (1992).
    
    The set of temperatures to sample can be specified in two ways.  First, you can explicitly provide a list
    of temperatures by using the "temperatures" argument.  Alternatively, you can specify the minimum and
    maximum temperatures, and the total number of temperatures to use.  The temperatures are chosen spaced
    exponentially between the two extremes.  For example,
    
    st = SimulatedTempering(simulation, numTemperatures=15, minTemperature=300*kelvin, maxTemperature=450*kelvin)
    
    After creating the SimulatedTempering object, call step() on it to run the simulation.
    
    Transitions between temperatures are performed at regular intervals, as specified by the "tempChangeInterval"
    argument.  For each transition, a new temperature is selected using the independence sampling method, as
    described in Chodera, J. and Shirts, M., J. Chem. Phys. 135, 194110 (2011).
    
    Simulated tempering requires a "weight factor" for each temperature.  Ideally, these should be chosen so
    the simulation spends equal time at every temperature.  You can specify the list of weights to use with the
    optional "weights" argument.  If this is omitted, weights are selected automatically using the Wang-Landau
    algorithm as described in Wang, F. and Landau, D. P., Phys. Rev. Lett. 86(10), pp. 2050-2053 (2001).
    
    To properly analyze the results of the simulation, it is important to know the temperature and weight factors
    at every point in time.  The SimulatedTempering object functions as a reporter, writing this information
    to a file or stdout at regular intervals (which should match the interval at which you save frames from the
    simulation).  You can specify the output file and reporting interval with the "reportFile" and "reportInterval"
    arguments.
    """

    def __init__(self, simulation, temperatures=None, numTemperatures=None, minTemperature=None, maxTemperature=None, weights=None, tempChangeInterval=25, reportInterval=1000, reportFile=stdout):
        """Create a new SimulatedTempering.
        
        Parameters
        ----------
        simulation: Simulation
            The Simulation defining the System, Context, and Integrator to use
        temperatures: list
            The list of temperatures to use for tempering, in increasing order
        numTemperatures: int
            The number of temperatures to use for tempering.  If temperatures is not None, this is ignored.
        minTemperature: temperature
            The minimum temperature to use for tempering.  If temperatures is not None, this is ignored.
        maxTemperature: temperature
            The maximum temperature to use for tempering.  If temperatures is not None, this is ignored.
        weights: list
            The weight factor for each temperature.  If none, weights are selected automatically.
        tempChangeInterval: int
            The interval (in time steps) at which to attempt transitions between temperatures
        reportInterval: int
            The interval (in time steps) at which to write information to the report file
        reportFile: string or file
            The file to write reporting information to, specified as a file name or file object
        """        
        self.simulation = simulation
        if temperatures is None:
            if unit.is_quantity(minTemperature):
                minTemperature = minTemperature.value_in_unit(unit.kelvin)
            if unit.is_quantity(maxTemperature):
                maxTemperature = maxTemperature.value_in_unit(unit.kelvin)
            self.temperatures = [minTemperature*((float(maxTemperature)/minTemperature)**(i/float(numTemperatures-1))) for i in range(numTemperatures)]*unit.kelvin
        else:
            numTemperatures = len(temperatures)
            self.temperatures = [(t.value_in_unit(unit.kelvin) if unit.is_quantity(t) else t)*unit.kelvin for t in temperatures]
            if any(self.temperatures[i] >= self.temperatures[i+1] for i in range(numTemperatures-1)):
                raise ValueError('The temperatures must be in strictly increasing order')
        self.tempChangeInterval = tempChangeInterval
        self.reportInterval = reportInterval
        self.inverseTemperatures = [1.0/(unit.MOLAR_GAS_CONSTANT_R*t) for t in self.temperatures]

        # If necessary, open the file we will write reports to.

        self._openedFile = isinstance(reportFile, str)
        if self._openedFile:
            # Detect the desired compression scheme from the filename extension
            # and open all files unbuffered
            if reportFile.endswith('.gz'):
                if not have_gzip:
                    raise RuntimeError("Cannot write .gz file because Python could not import gzip library")
                self._out = gzip.GzipFile(fileobj=open(reportFile, 'wb', 0))
            elif reportFile.endswith('.bz2'):
                if not have_bz2:
                    raise RuntimeError("Cannot write .bz2 file because Python could not import bz2 library")
                self._out = bz2.BZ2File(reportFile, 'w', 0)
            else:
                self._out = open(reportFile, 'w', 1)
        else:
            self._out = reportFile
        
        # Initialize the weights.
        
        if weights is None:
            self._weights = [0.0]*numTemperatures
            self._updateWeights = True
            self._weightUpdateFactor = 1.0
            self._histogram = [0]*numTemperatures
            self._hasMadeTransition = False
        else:
            self._weights = weights
            self._updateWeights = False

        # Select the initial temperature.
        
        self.currentTemperature = 0
        self.simulation.integrator.setTemperature(self.temperatures[self.currentTemperature])
        
        # Add a reporter to the simulation which will handle the updates and reports.
        
        class STReporter(object):
            def __init__(self, st):
                self.st = st

            def describeNextReport(self, simulation):
                st = self.st
                steps1 = st.tempChangeInterval - simulation.currentStep%st.tempChangeInterval
                steps2 = st.reportInterval - simulation.currentStep%st.reportInterval
                steps = min(steps1, steps2)
                isUpdateAttempt = (steps1 == steps)
                return (steps, False, isUpdateAttempt, False, isUpdateAttempt)

            def report(self, simulation, state):
                st = self.st
                if simulation.currentStep%st.tempChangeInterval == 0:
                    st._attemptTemperatureChange(state)
                if simulation.currentStep%st.reportInterval == 0:
                    st._writeReport()
        
        simulation.reporters.append(STReporter(self))
        
        # Write out the header line.
        
        headers = ['Steps', 'Temperature (K)']
        for t in self.temperatures:
            headers.append('%gK Weight' % t.value_in_unit(unit.kelvin))
        print('#"%s"' % ('"\t"').join(headers), file=self._out)

    def __del__(self):
        if self._openedFile:
            self._out.close()
    
    @property
    def weights(self):
        return [x-self._weights[0] for x in self._weights]

    def step(self, steps):
        """Advance the simulation by integrating a specified number of time steps."""
        self.simulation.step(steps)
    
    def _attemptTemperatureChange(self, state):
        """Attempt to move to a different temperature."""
        
        # Compute the probability for each temperature.  This is done in log space to avoid overflow.

        logProbability = [(self._weights[i]-self.inverseTemperatures[i]*state.getPotentialEnergy()) for i in range(len(self._weights))]
        maxLogProb = max(logProbability)
        offset = maxLogProb + math.log(sum(math.exp(x-maxLogProb) for x in logProbability))
        probability = [math.exp(x-offset) for x in logProbability]
        r = random.random()
        for j in range(len(probability)):
            if r < probability[j]:
                if j != self.currentTemperature:
                    # Rescale the velocities.
                    
                    scale = math.sqrt(self.temperatures[j]/self.temperatures[self.currentTemperature])
                    velocities = [v*scale for v in state.getVelocities().value_in_unit(unit.nanometers/unit.picoseconds)]
                    self.simulation.context.setVelocities(velocities)

                    # Select this temperature.
                    
                    self._hasMadeTransition = True
                    self.currentTemperature = j
                    self.simulation.integrator.setTemperature(self.temperatures[j])
                if self._updateWeights:
                    # Update the weight factors.
                    
                    self._weights[j] -= self._weightUpdateFactor
                    self._histogram[j] += 1
                    minCounts = min(self._histogram)
                    if minCounts > 20 and minCounts >= 0.2*sum(self._histogram)/len(self._histogram):
                        # Reduce the weight update factor and reset the histogram.
                        
                        self._weightUpdateFactor *= 0.5
                        self._histogram = [0]*len(self.temperatures)
                        self._weights = [x-self._weights[0] for x in self._weights]
                    elif not self._hasMadeTransition and probability[self.currentTemperature] > 0.99:
                        # Rapidly increase the weight update factor at the start of the simulation to find
                        # a reasonable starting value.
                        
                        self._weightUpdateFactor *= 2.0
                        self._histogram = [0]*len(self.temperatures)
                return
            r -= probability[j]

    def _writeReport(self):
        """Write out a line to the report."""
        temperature = self.temperatures[self.currentTemperature].value_in_unit(unit.kelvin)
        values = [temperature]+self.weights
        print(('%d\t' % self.simulation.currentStep) + '\t'.join('%g' % v for v in values), file=self._out)

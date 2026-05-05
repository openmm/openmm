from __future__ import print_function

"""
simulatedtempering.py: Implements simulated tempering

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
from sys import stdout


class SimulatedTempering(object):
    """SimulatedTempering implements the simulated tempering algorithm for accelerated sampling.

    @deprecated This class exists mainly for historical reasons.  It now is just a thin wrapper around ExpandedEnsembleSampler.  It is still supported for backward compatibility, but using ExpandedEnsembleSampler directly is recommended since it is much more flexible.

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
        if temperatures is None:
            if unit.is_quantity(minTemperature):
                minTemperature = minTemperature.value_in_unit(unit.kelvin)
            if unit.is_quantity(maxTemperature):
                maxTemperature = maxTemperature.value_in_unit(unit.kelvin)
            temperatures = [minTemperature*((float(maxTemperature)/minTemperature)**(i/float(numTemperatures-1))) for i in range(numTemperatures)]*unit.kelvin
        self.temperatures = [(t.value_in_unit(unit.kelvin) if unit.is_quantity(t) else t)*unit.kelvin for t in temperatures]
        states = [{'temperature':t} for t in temperatures]
        from . import ExpandedEnsembleSampler
        self._sampler = ExpandedEnsembleSampler(states, simulation, tempChangeInterval, weights=weights,
                                                reportInterval=reportInterval, logFile=reportFile)

    @property
    def weights(self):
        return self._sampler.weights

    @property
    def currentTemperature(self):
        return self._sampler.currentStateIndex

    def step(self, steps):
        """Advance the simulation by integrating a specified number of time steps."""
        self._sampler.step(steps)

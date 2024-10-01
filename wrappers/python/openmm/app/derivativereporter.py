"""
derivativereporter.py: Reports the derivative of the energy with regards to a global parameter

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2024 the Authors.
Authors: Marc Schuh
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
__author__ = "Marc Schuh"
__version__ = "1.0"

from openmm import unit

__all__ = ['DerivativeReporter']


class DerivativeReporter(object):
    """DerivativeReporter reports the derivative of the energy with regards to one or more global parameters.

    To use this reporter, the forces for which the energy derivative should be calculated need to have the option `addEnergyParameterDerivative` set for the corresponding global variable.
    
    Based on openMM source code for StateDataReporter, and openMM user guide http://docs.openmm.org/latest/userguide/application/04_advanced_sim_examples.html#extracting-and-reporting-forces-and-other-data

    """
    def __init__(self, file, reportInterval: int, derivativeNames: list, append: bool = False, timeUnit: unit = unit.picosecond):
        """Create a DerivativeReporter.

        Parameters
        ----------
        file : string or file
            The file to write to, specified as a file name or file object
        reportInterval : int
            The interval (in time steps) at which to report derivatives
        derivativeNames : list
            A list of global parameter names for which to report energy derivatives for
        append : bool=False
            If true, append to an existing file.  This has two effects.  First,
            the file is opened in append mode.  Second, the header line is not
            written, since there is assumed to already be a header line at the
            start of the file.
        timeUnit : unit=unit.picosecond
            The unit of the simulation time in the output.
        """        
        self._append = append
        self._derivativeNames = derivativeNames
        self._reportInterval = reportInterval
        self._timeUnit = timeUnit

        # Open the file if it is a file, otherwise write directly to it (probably stdout)
        self._openedFile = isinstance(file, str)
        if self._openedFile:
            self._out = open(file, 'a' if append else 'w')
        else:
            self._out = file

        self._hasInitialized = False

    def __del__(self):
        if self._openedFile:
            self._out.close()

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        dict
            A dictionary describing the required information for the next report
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return {'steps':steps, 'periodic':None, 'include':['parameterDerivatives']}

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        # print header if it is the first time this object has been called and is not set to append.
        if not self._hasInitialized:
            self._hasInitialized = True
            if not self._append:
                header = ','.join(self._derivativeNames)
                header = "time," + header

                print(header, file=self._out)
                try:
                    self._out.flush()
                except AttributeError:
                    pass

        derivatives = [str(state.getEnergyParameterDerivatives()[derivativeName]) for derivativeName in self._derivativeNames]

        values = ','.join(derivatives)
        values = f"{state.getTime().value_in_unit(self._timeUnit): .3f},{values}"

        print(values, file=self._out)
        try:
            self._out.flush()
        except AttributeError:
            pass
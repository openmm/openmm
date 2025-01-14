"""
debuggingreporter.py: Used for debugging hard to reproduce errors

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2024 Stanford University and the Authors.
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

class DebuggingReporter(object):
    """DebuggingReporter is a useful tool for debugging simulations that blow up
    after running for a while.  It records a checkpoint after every time step, keeping
    a list of the most recent ones in memory.  It also checks the positions,
    velocities, and forces to make sure they are still finite.  If it discovers that 
    anything has become NaN or infinite, it writes its saved checkpoints to disk and
    throws an exception to immediately end the simulation.  This allows you to
    inspect the behavior in the steps leading up to the error and retrace how it
    happened.

    To use this class, create a DebuggingReporter and add it to the Simulation's list of reporters.
    """

    def __init__(self, numCheckpoints=10, prefix='debug'):
        """Create a DebuggingReporter.

        Parameters
        ----------
        states : int
            The number of States to keep in memory and write to disk when an
            error happens.
        prefix : string
            The prefix to use for the saved State files.
        """
        self._numCheckpoints = numCheckpoints
        self._prefix = prefix
        self._checkpoints = []

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
        return {'steps':1, 'periodic':False, 'include':['positions', 'velocities', 'forces']}

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if len(self._checkpoints) >= self._numCheckpoints:
            del self._checkpoints[0]
        self._checkpoints.append(simulation.context.createCheckpoint())
        import numpy as np
        if not np.all(np.isfinite(state.getPositions(asNumpy=True)._value)) or \
           not np.all(np.isfinite(state.getVelocities(asNumpy=True)._value)) or \
           not np.all(np.isfinite(state.getForces(asNumpy=True)._value)):
            for i, c in enumerate(self._checkpoints):
                with open(f'{self._prefix}{i+1}.chk', 'wb') as f:
                    f.write(c)
            raise Exception('Debugging reporter detected error')

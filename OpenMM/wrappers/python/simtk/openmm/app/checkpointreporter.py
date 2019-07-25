"""
checkpointreporter.py: Saves checkpoint files for a simulation

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2014-2016 Stanford University and the Authors.
Authors: Robert McGibbon
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
__author__ = "Robert McGibbon"
__version__ = "1.0"

import simtk.openmm as mm
import os
import os.path

__all__ = ['CheckpointReporter']


class CheckpointReporter(object):
    """CheckpointReporter saves periodic checkpoints of a simulation.
    The checkpoints will overwrite one another -- only the last checkpoint
    will be saved in the file.

    To use it, create a CheckpointReporter, then add it to the Simulation's
    list of reporters. To load a checkpoint file and continue a simulation,
    use the following recipe:

    >>> with open('checkput.chk', 'rb') as f:
    >>>     simulation.context.loadCheckpoint(f.read())

    Notes:
    A checkpoint contains not only publicly visible data such as the particle
    positions and velocities, but also internal data such as the states of
    random number generators.  Ideally, loading a checkpoint should restore the
    Context to an identical state to when it was written, such that continuing
    the simulation will produce an identical trajectory.  This is not strictly
    guaranteed to be true, however, and should not be relied on.  For most
    purposes, however, the internal state should be close enough to be
    reasonably considered equivalent.

    A checkpoint contains data that is highly specific to the Context from
    which it was created. It depends on the details of the System, the
    Platform being used, and the hardware and software of the computer it was
    created on.  If you try to load it on a computer with different hardware,
    or for a System that is different in any way, loading is likely to fail.
    Checkpoints created with different versions of OpenMM are also often
    incompatible.  If a checkpoint cannot be loaded, that is signaled by
    throwing an exception.

    """
    def __init__(self, file, reportInterval):
        """Create a CheckpointReporter.

        Parameters
        ----------
        file : string or open file object
            The file to write to. Any current contents will be overwritten.
        reportInterval : int
            The interval (in time steps) at which to write checkpoints.
        """

        self._reportInterval = reportInterval
        self._file = file

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A five element tuple. The first element is the number of steps
            until the next report. The remaining elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, False, False, False)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if isinstance(self._file, str):
            # Do a safe save.

            tempFilename1 = self._file+".backup1"
            tempFilename2 = self._file+".backup2"
            with open(tempFilename1, 'w+b', 0) as out:
                out.write(simulation.context.createCheckpoint())
            exists = os.path.exists(self._file)
            if exists:
                os.rename(self._file, tempFilename2)
            os.rename(tempFilename1, self._file)
            if exists:
                os.remove(tempFilename2)
        else:
            # Replace the contents of the file.

            self._file.seek(0)
            chk = simulation.context.createCheckpoint()
            self._file.write(chk)
            self._file.truncate()
            self._file.flush()

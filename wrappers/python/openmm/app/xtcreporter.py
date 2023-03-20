__author__ = "Raul P. Pelaez"
from openmm.app import XTCFile
from openmm.unit import nanometers, femtoseconds
import numpy as np


class XTCReporter(object):
    """Report simulation progress in the Gromacs XTC format.

    Parameters
    ----------
    file : str
        The file to write to
    reportInterval : int
        The interval (in time steps) at which to write frames
    enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
    """

    def __init__(self, file, reportInterval, enforcePeriodicBox=None):
        self._file = XTCFile(file)
        self._reportInterval = reportInterval
        self._enforcePeriodicBox = enforcePeriodicBox

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A six element tuple. The first element is the number of steps
            until the next report. The next four elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.  The final element specifies whether
            positions should be wrapped to lie in a single periodic box.
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False, self._enforcePeriodicBox)

    def report(self, simulation, state):
        """Generate a report.
        Stores the current positions, box vectors, time, and step number in the XTC file.
        Positions are written in nanometers, box vectors in nanometers, time in femtoseconds, and step number as an integer.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        step = simulation.currentStep
        coords = state.getPositions().value_in_unit(nanometers)
        box = state.getPeriodicBoxVectors().value_in_unit(nanometers)
        if (
            box[0][1] != 0
            or box[0][2] != 0
            or box[1][0] != 0
            or box[1][2] != 0
            or box[2][0] != 0
            or box[2][1] != 0
        ):
            raise Exception("XTCReporter does not support triclinic boxes")
        if box is None:
            box = np.array([0, 0, 0], dtype=np.float32)
        else:
            box = np.array([box[0][0], box[1][1], box[2][2]], dtype=np.float32)
        time = state.getTime().value_in_unit(femtoseconds)
        self._file.writeFrame(coords, box, time, step)

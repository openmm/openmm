from sys import _current_frames
from openmm.app import XTCFile, Topology
from openmm.unit import nanometers, femtoseconds



class XTCReporter(object):
    def __init__(self, file, reportInterval, enforcePeriodicBox=None):
        self._file = XTCFile(file)
        self._reportInterval = reportInterval
        self._enforcePeriodicBox = enforcePeriodicBox
        self._step = 0

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
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, True, False, False, False, self._enforcePeriodicBox)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        self._step += 1
        coords = state.getPositions().value_in_unit(nanometers)
        box = state.getPeriodicBoxVectors().value_in_unit(nanometers)
        time = state.getTime().value_in_unit(femtoseconds)
        step = self._step
        self._file.write_frame(coords, box, time, step)


    def __del__(self):
        self._file.close()

__author__ = "Raul P. Pelaez"
from openmm.app import XTCFile

class XTCReporter(object):
    """XTCReporter outputs a series of frames from a Simulation to a XTC file.

    To use it, create a XTCReporter, then add it to the Simulation's list of reporters.
    """

    def __init__(self, file, reportInterval, append=False, enforcePeriodicBox=None):
        """Create a XTCReporter.

        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        append : bool=False
            If True, open an existing XTC file to append to.  If False, create a new file.
        enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
        """
        self._reportInterval = reportInterval
        self._append = append
        self._enforcePeriodicBox = enforcePeriodicBox
        self._fileName = file
        self._xtc = None
        if not append:
            open(file, 'wb').close()

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
        return {'steps':steps, 'periodic':self._enforcePeriodicBox, 'include':['positions']}

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """

        if self._xtc is None:
            self._xtc = XTCFile(
                self._fileName,
                simulation.topology,
                simulation.integrator.getStepSize(),
                simulation.currentStep,
                self._reportInterval,
                self._append,
            )
        self._xtc.writeModel(
            state.getPositions(asNumpy=True), periodicBoxVectors=state.getPeriodicBoxVectors()
        )

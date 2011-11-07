"""
pdbreporter.py: Outputs simulation trajectories in PDB format
"""
__author__ = "Peter Eastman"
__version__ = "1.0"

import simtk.openmm as mm
from simtk.openmm.app import PDBFile
    
class PDBReporter(object):
    """PDBReporter outputs a series of frames from a Simulation to a PDB file.
    
    To use it, create a PDBReporter, than add it to the Simulation's list of reporters.
    """
    
    def __init__(self, file, reportInterval):
        """Create a PDBReporter.
    
        Parameters:
         - file (string) The file to write to
         - reportInterval (int) The interval (in time steps) at which to write frames
        """
        self._reportInterval = reportInterval
        self._out = open(file, 'w')
        self._topology = None
        self._nextModel = 0
    
    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.
        
        Parameters:
         - simulation (Simulation) The Simulation to generate a report for
        Returns: A five element tuple.  The first element is the number of steps until the
        next report.  The remaining elements specify whether that report will require
        positions, velocities, forces, and energies respectively.
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, True, False, False, False)
    
    def report(self, simulation, state):
        """Generate a report.
        
        Parameters:
         - simulation (Simulation) The Simulation to generate a report for
         - state (State) The current state of the simulation
        """
        if self._nextModel == 0:
            PDBFile.writeHeader(simulation.topology, self._out)
            self._topology = simulation.topology
            self._nextModel += 1
        PDBFile.writeModel(simulation.topology, state.getPositions(), self._out, self._nextModel)
        self._nextModel += 1
    
    def __del__(self):
        PDBFile.writeFooter(self._topology, self._out)
        self._out.close()

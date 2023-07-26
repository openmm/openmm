"""
reweightingreporter.py: Outputs reweighting factors evaluated within simulation using a LangevinSplittingGirsanov integrator implemented in openmmtools.
Author: Joana-Lysiane Schaefer, Freie Universität Berlin
"""

import openmm as mm
import openmm.unit as unit
import math
import time

## ToDo. enable more pertubations
##       include remaining time see State reporter script

class ReweightingReporter(object):
    """ReweightingReporter outputs the energy of the perturbation U.
    
    To use it, create a PertubationEnergyReporter, then add it to the Simulation's list of reporters.  The set of
    data to write is configurable using boolean flags passed to the constructor.  By default the data is
    written in comma-separated-value (CSV) format, but you can specify a different separator to use.
    """

    def __init__(self, file, reportInterval, integrator, factorM=True, unperturebed=False, firtsPertubation=False, secondPertubation=False, thirdPertubation=False, fourthPertubation=False, separator=','):
        """Create a PertubationEnergyReporter.
        To print the energy at a simulation step simulation.context.getState() is used. with argument groups one of the force
        groups can be selected. Depending on the pertubation defined a number of groups are of interesst and can be choosen with
        the number sequence (after b) inputed to groups : "groups = 0b00000000000000000000000000000000".
        
        Parameters
        ----------
        file : string or file
            The file to write to, specified as a file name or file object
        reportInterval : int
            The interval (in time steps) at which to write frames
        unperturebed : bool=False
            Whether to write the unperturbed energy to the file, groups argument : "0b00000000000000000000000000000001"
        firtsPertubation : bool=False
            Whether to write the first pertubation energy to the file, groups argument : "0b00000000000000000000000000000010"
        secondPertubation : bool=False
            Whether to write the second pertubation energy to the file, groups argument : "0b00000000000000000000000000000100"
        thirdPertubation : bool=False
            Whether to write the third pertubation energy to the file, groups argument : "0b00000000000000000000000000001000"
        fourthPertubation : bool=False
            Whether to write the fourth pertubation energy to the file, groups argument : "0b00000000000000000000000000010000"
        separator : string=','
            The separator to use between columns in the file
        """
        self._reportInterval = reportInterval
        self._integrator = integrator
        self._out = open(file, 'w')
        
        
        self._factorM = factorM
        self._unperturebed = unperturebed
        self._firtsPertubation = firtsPertubation
        self._secondPertubation = secondPertubation
        self._thirdPertubation = thirdPertubation
        self._fourthPertubation = fourthPertubation
        self._separator = separator
        self._hasInitialized = False
        self._needsPositions = False
        self._needsVelocities = False
        self._needsForces = True
        self._reportForces = False
        self._needEnergy = True

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
        return (steps, self._needsPositions, self._needsVelocities, self._needsForces, self._needEnergy)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if not self._hasInitialized:
            headers = self._constructHeaders()
            print('#"%s"' % ('"'+self._separator+'"').join(headers), file=self._out)
            try:
                self._out.flush()
            except AttributeError:
                pass
            self._initialClockTime = time.time()
            self._initialSimulationTime = state.getTime()
            self._initialSteps = simulation.currentStep
            self._hasInitialized = True

        # Check for errors.
        self._checkForErrors(simulation, state)

        # Query for the values
        values = self._constructReportValues(simulation, state)

        # Write the values.
        print(self._separator.join(str(v) for v in values), file=self._out)
        try:
            self._out.flush()
        except AttributeError:
            pass

    def _constructReportValues(self, simulation, state):
        """Query the simulation for the current state of our observables of interest.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation

        Returns
        -------
        A list of values summarizing the current state of
        the simulation, to be printed or saved. Each element in the list
        corresponds to one of the columns in the resulting CSV file.
        """
        values = []
        clockTime = time.time()
        if self._factorM:
            #values.append(self._integrator.getPerDofVariableByName("TESTPDOF"))  ## DELETE ME
            values.append(self._integrator.getGlobalVariableByName("M"))
        if self._unperturebed:
            state  = simulation.context.getState(getPositions = self._needsPositions, getVelocities = False,
                                                getForces = False, getEnergy = True,
                                                getParameters = False, enforcePeriodicBox = False,
                                                groups = 0b00000000000000000000000000000000) 
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            values.append(energy)
        if self._firtsPertubation:
            state  = simulation.context.getState(getPositions = self._needsPositions, getVelocities = False,
                                                getForces = False, getEnergy = True,
                                                getParameters = False, enforcePeriodicBox = False,
                                                groups = 0b00000000000000000000000000000010) 
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            values.append(energy)            
        if self._secondPertubation:
            state  = simulation.context.getState(getPositions = self._needsPositions, getVelocities = False,
                                                getForces = False, getEnergy = True,
                                                getParameters = False, enforcePeriodicBox = False,
                                                groups = 0b00000000000000000000000000000100) 
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            values.append(energy)        
        if self._thirdPertubation:
            state  = simulation.context.getState(getPositions = self._needsPositions, getVelocities = False,
                                                getForces = False, getEnergy = True,
                                                getParameters = False, enforcePeriodicBox = False,
                                                groups = 0b00000000000000000000000000001000) 
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            values.append(energy)
        if self._fourthPertubation:
            state  = simulation.context.getState(getPositions = self._needsPositions, getVelocities = False,
                                                getForces = False, getEnergy = True,
                                                getParameters = False, enforcePeriodicBox = False,
                                                groups = 0b00000000000000000000000000001000) 
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            values.append(energy)            
        return values


    def _constructHeaders(self):
        """Construct the headers for the CSV output

        Returns: a list of strings giving the title of each observable being reported on.
        """
        headers = []
        if self._factorM:
            headers.append('Reweighting Factor M (no exp())')
        if self._unperturebed:
            headers.append('Unperturebed Energy (kJ/mole)')
        if self._firtsPertubation:
            headers.append('First Pertubation (kJ/mole)')
        if self._secondPertubation:
            headers.append('Second Pertubation (kJ/mole)')
        if self._thirdPertubation:
            headers.append('Third Pertubation (kJ/mole)')
        if self._fourthPertubation:
            headers.append('Fourth Pertubation (kJ/mole)')
        return headers

    def _checkForErrors(self, simulation, state):
        """Check for errors in the current state of the simulation

        Parameters
         - simulation (Simulation) The Simulation to generate a report for
         - state (State) The current state of the simulation
        """
        if self._needEnergy:
            energy = (state.getKineticEnergy()+state.getPotentialEnergy()).value_in_unit(unit.kilojoules_per_mole)
            if math.isnan(energy):
                raise ValueError('Energy is NaN')
            if math.isinf(energy):
                raise ValueError('Energy is infinite')

  

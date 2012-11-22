"""
statedatareporter.py: Outputs data about a simulation

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
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

import simtk.openmm as mm
import simtk.unit as unit
from simtk.openmm.app import PDBFile
import math
    
class StateDataReporter(object):
    """StateDataReporter outputs information about a simulation, such as energy and temperature, to a file.
    
    To use it, create a StateDataReporter, then add it to the Simulation's list of reporters.  The set of
    data to write is configurable using boolean flags passed to the constructor.  By default the data is
    written in comma-separated-value (CSV) format, but you can specify a different separator to use.
    """
    
    def __init__(self, file, reportInterval, step=False, time=False, potentialEnergy=False, kineticEnergy=False, totalEnergy=False, temperature=False, volume=False, density=False, separator=',', systemMass=None):
        """Create a StateDataReporter.
    
        Parameters:
         - file (string or file) The file to write to, specified as a file name or file object
         - reportInterval (int) The interval (in time steps) at which to write frames
         - step (boolean=False) Whether to write the current step index to the file
         - time (boolean=False) Whether to write the current time to the file
         - potentialEnergy (boolean=False) Whether to write the potential energy to the file
         - kineticEnergy (boolean=False) Whether to write the kinetic energy to the file
         - totalEnergy (boolean=False) Whether to write the total energy to the file
         - temperature (boolean=False) Whether to write the instantaneous temperature to the file
         - volume (boolean=False) Whether to write the periodic box volume to the file
         - density (boolean=False) Whether to write the system density to the file
         - separator (string=',') The separator to use between columns in the file
         - systemMass (mass=None) The total mass to use for the system when reporting density.  If this is
           None (the default), the system mass is computed by summing the masses of all particles.  This
           parameter is useful when the particle masses do not reflect their actual physical mass, such as
           when some particles have had their masses set to 0 to immobilize them.
        """
        self._reportInterval = reportInterval
        self._openedFile = isinstance(file, str)
        if self._openedFile:
            self._out = open(file, 'w')
        else:
            self._out = file
        self._step = step
        self._time = time
        self._potentialEnergy = potentialEnergy
        self._kineticEnergy = kineticEnergy
        self._totalEnergy = totalEnergy
        self._temperature = temperature
        self._volume = volume
        self._density = density
        self._separator = separator
        self._needEnergy = potentialEnergy or kineticEnergy or totalEnergy or temperature
        self._totalMass = systemMass
        self._hasInitialized = False
    
    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.
        
        Parameters:
         - simulation (Simulation) The Simulation to generate a report for
        Returns: A five element tuple.  The first element is the number of steps until the
        next report.  The remaining elements specify whether that report will require
        positions, velocities, forces, and energies respectively.
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, False, False, self._needEnergy)
    
    def report(self, simulation, state):
        """Generate a report.
        
        Parameters:
         - simulation (Simulation) The Simulation to generate a report for
         - state (State) The current state of the simulation
        """
        if not self._hasInitialized:
            system = simulation.system
            if self._temperature:
                # Compute the number of degrees of freedom.
                dof = 0
                for i in range(system.getNumParticles()):
                    if system.getParticleMass(i) > 0*unit.dalton:
                        dof += 3
                dof -= system.getNumConstraints()
                if any(type(system.getForce(i)) == mm.CMMotionRemover for i in range(system.getNumForces())):
                    dof -= 3
                self._dof = dof
            if self._density:
                if self._totalMass is None:
                    # Compute the total system mass.
                    self._totalMass = 0*unit.dalton
                    for i in range(system.getNumParticles()):
                        self._totalMass += system.getParticleMass(i)
                elif not unit.is_quantity(self._totalMass):
                    self._totalMass = self._totalMass*unit.dalton
            
            # Write the headers.
            
            headers = []
            if self._step:
                headers.append('Step')
            if self._time:
                headers.append('Time (ps)')
            if self._potentialEnergy:
                headers.append('Potential Energy (kJ/mole)')
            if self._kineticEnergy:
                headers.append('Kinetic Energy (kJ/mole)')
            if self._totalEnergy:
                headers.append('Total Energy (kJ/mole)')
            if self._temperature:
                headers.append('Temperature (K)')
            if self._volume:
                headers.append('Box Volume (nm^3)')
            if self._density:
                headers.append('Density (g/mL)')
            print >>self._out, '#"%s"' % ('"'+self._separator+'"').join(headers)
            self._hasInitialized = True

        # Check for errors.
        
        if self._needEnergy:
            energy = (state.getKineticEnergy()+state.getPotentialEnergy()).value_in_unit(unit.kilojoules_per_mole)
            if math.isnan(energy):
                raise ValueError('Energy is NaN')
            if math.isinf(energy):
                raise ValueError('Energy is infinite')
                
        # Write the values.
        
        values = []
        box = state.getPeriodicBoxVectors()
        volume = box[0][0]*box[1][1]*box[2][2]
        if self._step:
            values.append(simulation.currentStep)
        if self._time:
            values.append(state.getTime().value_in_unit(unit.picosecond))
        if self._potentialEnergy:
            values.append(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
        if self._kineticEnergy:
            values.append(state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole))
        if self._totalEnergy:
            values.append((state.getKineticEnergy()+state.getPotentialEnergy()).value_in_unit(unit.kilojoules_per_mole))
        if self._temperature:
            values.append((2*state.getKineticEnergy()/(self._dof*0.00831451)).value_in_unit(unit.kilojoules_per_mole))
        if self._volume:
            values.append(volume.value_in_unit(unit.nanometer**3))
        if self._density:
            values.append((self._totalMass/volume).value_in_unit(unit.gram/unit.item/unit.milliliter))
        print >>self._out, self._separator.join(str(v) for v in values)
        
    def __del__(self):
        if self._openedFile:
            self._out.close()

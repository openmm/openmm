# -*- coding: UTF-8 -*-
"""
nhcintegrator.py: Implements the Nosé-Hoover chain thermostated integration method.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

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
from __future__ import absolute_import, print_function
__author__ = "Andrew C. Simmonett"
__version__ = "1.0"

from simtk.openmm import CustomIntegrator
from simtk.unit import Quantity, kilojoules_per_mole, nanometers, picoseconds, kelvin, amu
from simtk.unit.constants import MOLAR_GAS_CONSTANT_R
import math

class NHCIntegrator(CustomIntegrator):
    """ Implements the Nosé-Hoover thermostat, using the reversible multi time step velocity
        Verlet algorithm detailed in this paper:

            "Explicit reversible integrators for extended systems dynamics",
            G. J. Martyna, M. E. Tuckerman, D. J. Tobias, and Michael L. Klein,
            Mol. Phys., 87, 1117 (1996)

        More details about the method can be found here:

            "Nosé–Hoover chains: The canonical ensemble via continuous dynamics
            G. J. Martyna, M. L. Klein, and M. Tuckerman
            J. Chem. Phys., 97, 2635 (1992)


        Parameters:
        -----------

        timestep: unit.Quantity compatible with femtoseconds, default=1*simtk.unit.femtoseconds
            The integration timestep for particles.

        temperature: unit.Quantity compatible with kelvin, default=298*simtk.unit.kelvin
            The target temperature for the thermostat.

        chainlength: integer, default=5
            The number of thermostat particles in the Nosé-Hoover chain.  Increasing
            this parameter will affect computational cost, but will make the simulation
            less sensitive to thermostat parameters, particularly in stiff systems; see
            the 1992 paper referenced above for more details.

        oscillatorperiod: unit.Quantity compatible with picoseconds, default=0.5*simtk.unit.picoseconds
            The period of motion of the thermostat particles.  A very large value will result
            in a distribution approaching the microcanonical ensemble, while a small value will
            cause rapid fluctuations in the temperature before convergence.

        n_mts: integer, default=5
            The number of timesteps used in the multi-timestep procedure to update the thermostat
            positions.  A higher value will increase the stability of the dynamics, but will also
            increase the compuational cost of the integration.

        n_yoshidasuzuki: integer, default=5
            The number of Yoshida-Suzuki steps used to subdivide each of the multi-timesteps used
            to update the thermostat positions.  A higher value will increase the stability of the 
            dynamics, but will also increase the computational cost of the integration; only certain
            values (currently 1,3, or 5) are supported, because weights must be defined.
    """

    YSWeights = {
        1 : [ 1.0000000000000000 ],
        3 : [ 0.8289815435887510, -0.6579630871775020,  0.8289815435887510 ],
        5 : [ 0.2967324292201065,  0.2967324292201065, -0.1869297168804260, 0.2967324292201065, 0.2967324292201065 ]
    }

    def get_ndf(self):
        """ Returns the number of degrees of freedom, 3*natom, ignoring constraints. """
        return int(self.getGlobalVariable(0))

    def get_ke(self):
        """ Returns the particles' kinetic energy. """
        return self.getGlobalVariable(1) * kilojoules_per_mole

    def get_pe(self):
        """ Returns the particles' potential energy. """
        return self.getGlobalVariable(2) * kilojoules_per_mole

    def get_bath_ke(self):
        """ Returns the thermostat bath kinetic energy. """
        return self.getGlobalVariable(3) * kilojoules_per_mole

    def get_bath_pe(self):
        """ Returns the thermostat bath potential energy. """
        return self.getGlobalVariable(4) * kilojoules_per_mole

    def get_conserved_energy(self):
        """ Returns the NHC conserved energy, which is the sum
            of all particle and thermostat energies. """
        return self.getGlobalVariable(5) * kilojoules_per_mole

    def get_current_temperature(self):
        """ Returns the instantaneous temperature of the system. """
        return self.getGlobalVariable(6) * kelvin

    def __init__(self, timestep=0.001*picoseconds, temperature=298*kelvin, chainlength=5,
                 oscillatorperiod=0.5*picoseconds, n_mts=5, n_yoshidasuzuki=5):
        CustomIntegrator.__init__(self, timestep)
        self.n_c         = n_mts
        self.n_ys        = n_yoshidasuzuki
        try:
            self.weights     = self.YSWeights[self.n_ys]
        except KeyError:
            raise Exception("Invalid Yoshida-Suzuki value. Allowed values are: %s"%
                             ",".join(map(str,self.YSWeights.keys())))
        if chainlength < 0:
            raise Exception("Nosé-Hoover chain length must be at least 0")
        if chainlength == 0:
            print("WARNING: Nosé-Hoover chain length is 0; falling back to regular velocity verlet algorithm.")
        self.M           = chainlength
        self.R           = MOLAR_GAS_CONSTANT_R.value_in_unit(kilojoules_per_mole/kelvin)
        self.RT          = self.R * temperature.value_in_unit(kelvin)

        # Define the "mass" of the thermostat particles (multiply by ndf for particle 0)
        frequency = (2*math.pi/oscillatorperiod).value_in_unit(picoseconds**-1)
        Q = self.RT/frequency**2

        #
        # Define global variables
        #
        self.addGlobalVariable("ndf", 0)           # 0 number of degrees of freedom
        self.addGlobalVariable("PKE", 0.0)         # 1 Particle kinetic energy
        self.addGlobalVariable("PPE", 0.0)         # 2 Particle potential energy
        self.addGlobalVariable("BKE", 0.0)         # 3 Thermostat bath kinetic energy
        self.addGlobalVariable("BPE", 0.0)         # 4 Thermostat bath potential energy
        self.addGlobalVariable("conserved", 0.0)   # 5 Conserved energy
        self.addGlobalVariable("temperature", 0.0) # 6 Instantaneous temperature
        self.addGlobalVariable("RT", self.RT)      # 7 RT 
        self.addGlobalVariable("R", self.R)        # 8 R
        self.addGlobalVariable("KE2", 0.0)         # 9 Twice the kinetic energy
        self.addGlobalVariable("Q", Q)
        self.addGlobalVariable("scale", 1.0)     
        self.addGlobalVariable("aa", 0.0)
        self.addGlobalVariable("wdt", 0.0)
        for w in range(self.n_ys):
            self.addGlobalVariable("w%d"%w, self.weights[w])

        #
        # Initialize thermostat parameters
        #
        for i in range(self.M):
            self.addGlobalVariable("xi%d"%i, 0)         # Thermostat particle
            self.addGlobalVariable("vxi%d"%i, 0)        # Thermostat particle velocities in ps^-1
            self.addGlobalVariable("G%d"%i, -self.RT/Q) # Forces on thermostat particles in ps^-2
            self.addGlobalVariable("Q%d"%i, 0)          # Thermostat "masses" in ps^2 kJ/mol
        # The masses need the number of degrees of freedom, which is approximated here.  Need a
        # better solution eventually, to properly account for constraints, translations, etc.
        self.addPerDofVariable("ones", 1.0)
        self.addPerDofVariable("x1", 0);
        self.addComputeSum("ndf", "ones")
        if self.M:
            self.addComputeGlobal("Q0", "ndf*Q")
            for i in range(1, self.M):
                self.addComputeGlobal("Q%d"%i, "Q")

        #
        # Take a step
        #
        if self.M:
            self.propagateNHC()

        # Velocity verlet algorithm to propagate particles
        self.addUpdateContextState()
        self.addComputePerDof("v", "v+0.5*dt*f/m")
        self.addComputePerDof("x", "x+dt*v")
        self.addComputePerDof("x1", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
        self.addConstrainVelocities()

        if self.M:
            self.propagateNHC()

        #
        # Compute statistics
        #
        self.compute_energies()

    def propagateNHC(self):
        """ Propagate the Nosé-Hoover chain """
        self.addComputeGlobal("scale", "1.0")
        self.addComputeSum("KE2", "m*v^2")
        self.addComputeGlobal("G0", "(KE2 - ndf*RT)/Q0")
        for ncval in range(self.n_c):
            for nysval in range(self.n_ys):
                self.addComputeGlobal("wdt", "w%d*dt/%d"%(nysval, self.n_c))
                self.addComputeGlobal("vxi%d"%(self.M-1), "vxi%d + 0.25*wdt*G%d"%(self.M-1, self.M-1))
                for j in range(self.M-2, -1, -1):
                    self.addComputeGlobal("aa", "exp(-0.125*wdt*vxi%d)"%(j+1))
                    self.addComputeGlobal("vxi%d"%j, "aa*(aa*vxi%d + 0.25*wdt*G%d)"%(j,j))
                # update particle velocities
                self.addComputeGlobal("aa", "exp(-0.5*wdt*vxi0)")
                self.addComputeGlobal("scale", "scale*aa")
                # update the thermostat positions
                for j in range(self.M):
                    self.addComputeGlobal("xi%d"%j, "xi%d + 0.5*wdt*vxi%d"%(j,j))
                # update the forces
                self.addComputeGlobal("G0", "(scale*scale*KE2 - ndf*RT)/Q0")
                # update thermostat velocities
                for j in range(self.M-1):
                    self.addComputeGlobal("aa", "exp(-0.125*wdt*vxi%d)"%(j+1))
                    self.addComputeGlobal("vxi%d"%j, "aa*(aa*vxi%d + 0.25*wdt*G%d)"%(j,j))
                    self.addComputeGlobal("G%d"%(j+1), "(Q%d*vxi%d*vxi%d - RT)/Q%d"%(j,j,j,j+1))
                self.addComputeGlobal("vxi%d"%(self.M-1), "vxi%d + 0.25*wdt*G%d"%(self.M-1, self.M-1))
        # update particle velocities
        self.addComputePerDof("v", "scale*v")

    def compute_energies(self):
        """ Computes kinetic and potential energies for particles and thermostat
            particles, as well as the conserved energy and current temperature """
        # Particle kinetic energy
        self.addComputeSum("PKE", "0.5*m*v^2")
        # Particle potential energy
        self.addComputeGlobal("PPE", "energy")
        # Bath kinetic energy
        self.addComputeGlobal("BKE", "0.0")
        for i in range(self.M):
            self.addComputeGlobal("BKE", "BKE + 0.5*Q%d*vxi%d^2"%(i,i))
        # Bath potential energy
        self.addComputeGlobal("BPE", "ndf*xi0")
        for i in range(1,self.M):
            self.addComputeGlobal("BPE", "BPE + xi%d"%i)
        self.addComputeGlobal("BPE", "RT*BPE")
        # Conserved energy is the sum of particle and bath energies
        self.addComputeGlobal("conserved", "PKE + PPE + BKE + BPE")
        # Instantaneous temperature
        self.addComputeGlobal("temperature", "2*PKE/(ndf*R)")



if __name__ == "__main__":
    import os
    from simtk.unit import angstrom, elementary_charge
    from simtk.openmm import System, Vec3, NonbondedForce, Context, LocalEnergyMinimizer, VerletIntegrator
    import numpy as np

    # This example was "borrowed" from John Chodera's really neat Argon free energy example

    nparticles      = 120                       # number of argon atoms
    mass            = 39.9 * amu                # mass
    sigma           = 3.4 * angstrom            # Lennard-Jones sigma
    epsilon         = 1.0 * kilojoules_per_mole # Lennard-Jones well-depth
    charge          = 0.0 * elementary_charge   # argon model has no charge
    reduced_density = 0.85                      # reduced density rho*sigma^3
    temperature     = 300 * kelvin              # temperature

    volume = nparticles*(sigma**3)/reduced_density
    box_edge = volume**(1.0/3.0)
    cutoff = min(box_edge*0.49, 2.5*sigma) # Compute cutoff, obeying minimum image

    print()
    print("sigma    = %.2f A" % sigma.value_in_unit(angstrom))
    print("box_edge = %.2f A" % box_edge.value_in_unit(angstrom))
    print("cutoff   = %.2f A" % cutoff.value_in_unit(angstrom))
    print()

    system = System()
    system.setDefaultPeriodicBoxVectors(Vec3(box_edge, 0, 0), Vec3(0, box_edge, 0), Vec3(0, 0, box_edge))
    force = NonbondedForce()
    force.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
    force.setCutoffDistance(cutoff) 
    for particle_index in range(nparticles):
        system.addParticle(mass)
        force.addParticle(charge, sigma, epsilon)
    system.addForce(force)

    # Create Integrator and Context.
    integrator = NHCIntegrator(temperature = 300*kelvin,
                               timestep=0.001*picoseconds,
                               chainlength=10,
                               oscillatorperiod=0.1*picoseconds,
                               n_mts=5,
                               n_yoshidasuzuki=5)

    context = Context(system, integrator)

    positions = Quantity(np.random.uniform(high=box_edge/angstrom, size=[nparticles,3]), angstrom)
    context.setPositions(positions)
    LocalEnergyMinimizer.minimize(context, tolerance=100*kilojoules_per_mole/nanometers)

    print("  Step          KE             PE          bath KE        bath PE      Temperature  Conserved Energy")
    print("=====================================================================================================")
    for iteration in range(20000):
        integrator.step(1)
        if iteration%20 == 0:
            PKE = integrator.get_ke().value_in_unit(kilojoules_per_mole)
            PPE = integrator.get_pe().value_in_unit(kilojoules_per_mole)
            BKE = integrator.get_bath_ke().value_in_unit(kilojoules_per_mole)
            BPE = integrator.get_bath_pe().value_in_unit(kilojoules_per_mole)
            conserved = PKE + PPE + BKE + BPE
            temp = integrator.get_current_temperature().value_in_unit(kelvin)
            print("%7d   %12.6f   %12.6f   %12.6f   %12.6f   %12.6f   %12.6f"%
                      (iteration, PKE, PPE, BKE, BPE, temp, conserved))

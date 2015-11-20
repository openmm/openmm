"""
amd.py: Implements the aMD integration method.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Peter Eastman, Steffen Lindert
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
__author__ = "Peter Eastman"
__version__ = "1.0"

from simtk.openmm import CustomIntegrator
from simtk.unit import kilojoules_per_mole, is_quantity

class AMDIntegrator(CustomIntegrator):
    """AMDIntegrator implements the aMD integration algorithm.

    The system is integrated based on a modified potential.  Whenever the energy V(r) is less than a
    cutoff value E, the following effective potential is used:

    V*(r) = V(r) + (E-V(r))^2 / (alpha+E-V(r))

    For details, see Hamelberg et al., J. Chem. Phys. 127, 155102 (2007).
    """

    def __init__(self, dt, alpha, E):
        """Create an AMDIntegrator.

        Parameters
        ----------
        dt : time
            The integration time step to use
        alpha : energy
            The alpha parameter to use
        E : energy
            The energy cutoff to use
        """
        CustomIntegrator.__init__(self, dt)
        self.addGlobalVariable("alpha", alpha)
        self.addGlobalVariable("E", E)
        self.addPerDofVariable("oldx", 0)
        self.addUpdateContextState();
        self.addComputePerDof("v", "v+dt*fprime/m; fprime=f*((1-modify) + modify*(alpha/(alpha+E-energy))^2); modify=step(E-energy)")
        self.addComputePerDof("oldx", "x")
        self.addComputePerDof("x", "x+dt*v")
        self.addConstrainPositions()
        self.addComputePerDof("v", "(x-oldx)/dt")

    def getAlpha(self):
        """Get the value of alpha for the integrator."""
        return self.getGlobalVariable(0)*kilojoules_per_mole

    def setAlpha(self, alpha):
        """Set the value of alpha for the integrator."""
        self.setGlobalVariable(0, alpha)

    def getE(self):
        """Get the energy threshold E for the integrator."""
        return self.getGlobalVariable(1)*kilojoules_per_mole

    def setE(self, E):
        """Set the energy threshold E for the integrator."""
        self.setGlobalVariable(1, E)

    def getEffectiveEnergy(self, energy):
        """Given the actual potential energy of the system, return the value of the effective potential."""
        alpha = self.getAlpha()
        E = self.getE()
        if not is_quantity(energy):
            energy = energy*kilojoules_per_mole # Assume kJ/mole
        if (energy > E):
            return energy
        return energy+(E-energy)*(E-energy)/(alpha+E-energy)


class AMDForceGroupIntegrator(CustomIntegrator):
    """AMDForceGroupIntegrator implements a single boost aMD integration algorithm.

    This is similar to AMDIntegrator, but is applied based on the energy of a single force group
    (typically representing torsions).

    For details, see Hamelberg et al., J. Chem. Phys. 127, 155102 (2007).
    """

    def __init__(self, dt, group, alphaGroup, EGroup):
        """Create a AMDForceGroupIntegrator.

        Parameters
        ----------
        dt : time
            The integration time step to use
        group : int
            The force group to apply the boost to
        alphaGroup : energy
            The alpha parameter to use for the boosted force group
        EGroup : energy
            The energy cutoff to use for the boosted force group
        """
        CustomIntegrator.__init__(self, dt)
        self.addGlobalVariable("alphaGroup", alphaGroup)
        self.addGlobalVariable("EGroup", EGroup)
        self.addGlobalVariable("groupEnergy", 0)
        self.addPerDofVariable("oldx", 0)
        self.addPerDofVariable("fg", 0)
        self.addUpdateContextState();
        self.addComputeGlobal("groupEnergy", "energy"+str(group))
        self.addComputePerDof("fg", "f"+str(group))
        self.addComputePerDof("v", "v+dt*fprime/m; fprime=fother + fg*((1-modify) + modify*(alphaGroup/(alphaGroup+EGroup-groupEnergy))^2); fother=f-fg; modify=step(EGroup-groupEnergy)")
        self.addComputePerDof("oldx", "x")
        self.addComputePerDof("x", "x+dt*v")
        self.addConstrainPositions()
        self.addComputePerDof("v", "(x-oldx)/dt")

    def getAlphaGroup(self):
        """Get the value of alpha for the boosted force group."""
        return self.getGlobalVariable(0)*kilojoules_per_mole

    def setAlphaGroup(self, alpha):
        """Set the value of alpha for the boosted force group."""
        self.setGlobalVariable(0, alpha)

    def getEGroup(self):
        """Get the energy threshold E for the boosted force group."""
        return self.getGlobalVariable(1)*kilojoules_per_mole

    def setEGroup(self, E):
        """Set the energy threshold E for the boosted force group."""
        self.setGlobalVariable(1, E)

    def getEffectiveEnergy(self, groupEnergy):
        """Given the actual group energy of the system, return the value of the effective potential.

        Parameters
        ----------
        groupEnergy : energy
            the actual potential energy of the boosted force group

        Returns
        -------
        energy
            the value of the effective potential
        """
        alphaGroup = self.getAlphaGroup()
        EGroup = self.getEGroup()
        if not is_quantity(groupEnergy):
            groupEnergy = groupEnergy*kilojoules_per_mole # Assume kJ/mole
        dE = 0.0*kilojoules_per_mole
        if (groupEnergy < EGroup):
            dE = dE + (EGroup-groupEnergy)*(EGroup-groupEnergy)/(alphaGroup+EGroup-groupEnergy)
        return groupEnergy+dE



class DualAMDIntegrator(CustomIntegrator):
    """DualAMDIntegrator implements a dual boost aMD integration algorithm.

    This is similar to AMDIntegrator, but two different boosts are applied to the potential:
    one based on the total energy, and one based on the energy of a single force group
    (typically representing torsions).

    For details, see Hamelberg et al., J. Chem. Phys. 127, 155102 (2007).
    """

    def __init__(self, dt, group, alphaTotal, ETotal, alphaGroup, EGroup):
        """Create a DualAMDIntegrator.

        Parameters
        ----------
        dt : time
            The integration time step to use
        group : int
            The force group to apply the second boost to
        alphaTotal : energy
            The alpha parameter to use for the total energy
        ETotal : energy
            The energy cutoff to use for the total energy
        alphaGroup : energy
            The alpha parameter to use for the boosted force group
        EGroup : energy
            The energy cutoff to use for the boosted force group
        """
        CustomIntegrator.__init__(self, dt)
        self.addGlobalVariable("alphaTotal", alphaTotal)
        self.addGlobalVariable("ETotal", ETotal)
        self.addGlobalVariable("alphaGroup", alphaGroup)
        self.addGlobalVariable("EGroup", EGroup)
        self.addGlobalVariable("groupEnergy", 0)
        self.addPerDofVariable("oldx", 0)
        self.addPerDofVariable("fg", 0)
        self.addUpdateContextState();
        self.addComputeGlobal("groupEnergy", "energy"+str(group))
        self.addComputePerDof("fg", "f"+str(group))
        self.addComputePerDof("v", """v+dt*fprime/m;
                                      fprime=fprime1 + fprime2;
                                      fprime2=fg*((1-modifyGroup) + modifyGroup*(alphaGroup/(alphaGroup+EGroup-groupEnergy))^2);
                                      fprime1=fother*((1-modifyTotal) + modifyTotal*(alphaTotal/(alphaTotal+ETotal-energy))^2);
                                      fother=f-fg;
                                      modifyTotal=step(ETotal-energy); modifyGroup=step(EGroup-groupEnergy)""")
        self.addComputePerDof("oldx", "x")
        self.addComputePerDof("x", "x+dt*v")
        self.addConstrainPositions()
        self.addComputePerDof("v", "(x-oldx)/dt")

    def getAlphaTotal(self):
        """Get the value of alpha for the total energy."""
        return self.getGlobalVariable(0)*kilojoules_per_mole

    def setAlphaTotal(self, alpha):
        """Set the value of alpha for the total energy."""
        self.setGlobalVariable(0, alpha)

    def getETotal(self):
        """Get the energy threshold E for the total energy."""
        return self.getGlobalVariable(1)*kilojoules_per_mole

    def setETotal(self, E):
        """Set the energy threshold E for the total energy."""
        self.setGlobalVariable(1, E)

    def getAlphaGroup(self):
        """Get the value of alpha for the boosted force group."""
        return self.getGlobalVariable(2)*kilojoules_per_mole

    def setAlphaGroup(self, alpha):
        """Set the value of alpha for the boosted force group."""
        self.setGlobalVariable(2, alpha)

    def getEGroup(self):
        """Get the energy threshold E for the boosted force group."""
        return self.getGlobalVariable(3)*kilojoules_per_mole

    def setEGroup(self, E):
        """Set the energy threshold E for the boosted force group."""
        self.setGlobalVariable(3, E)

    def getEffectiveEnergy(self, totalEnergy, groupEnergy):
        """Given the actual potential energy of the system, return the value of the effective potential.

        Parameters
        ----------
        totalEnergy : energy
            the actual potential energy of the whole system
        groupEnergy : energy
            the actual potential energy of the boosted force group

        Returns
        -------
        energy
            the value of the effective potential
        """
        alphaTotal = self.getAlphaTotal()
        ETotal = self.getETotal()
        alphaGroup = self.getAlphaGroup()
        EGroup = self.getEGroup()
        if not is_quantity(totalEnergy):
            totalEnergy = totalEnergy*kilojoules_per_mole # Assume kJ/mole
        if not is_quantity(groupEnergy):
            groupEnergy = groupEnergy*kilojoules_per_mole # Assume kJ/mole
        dE = 0.0*kilojoules_per_mole
        if (totalEnergy < ETotal):
            dE = dE + (ETotal-totalEnergy)*(ETotal-totalEnergy)/(alphaTotal+ETotal-totalEnergy)
        if (groupEnergy < EGroup):
            dE = dE + (EGroup-groupEnergy)*(EGroup-groupEnergy)/(alphaGroup+EGroup-groupEnergy)
        return totalEnergy+dE

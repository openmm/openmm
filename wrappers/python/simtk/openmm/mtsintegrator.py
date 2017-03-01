"""
respa.py: Implements the rRESPA multiple time step integration method.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2013-2015 Stanford University and the Authors.
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
from __future__ import absolute_import, division
__author__ = "Peter Eastman"
__version__ = "1.0"

from simtk.openmm import CustomIntegrator

class MTSIntegrator(CustomIntegrator):
    """MTSIntegrator implements the rRESPA multiple time step integration algorithm.

    This integrator allows different forces to be evaluated at different frequencies,
    for example to evaluate the expensive, slowly changing forces less frequently than
    the inexpensive, quickly changing forces.

    To use it, you must first divide your forces into two or more groups (by calling
    setForceGroup() on them) that should be evaluated at different frequencies.  When
    you create the integrator, you provide a tuple for each group specifying the index
    of the force group and the frequency (as a fraction of the outermost time step) at
    which to evaluate it.  For example::

        integrator = MTSIntegrator(4*femtoseconds, [(0,1), (1,2), (2,8)])

    This specifies that the outermost time step is 4 fs, so each step of the integrator
    will advance time by that much.  It also says that force group 0 should be evaluated
    once per time step, force group 1 should be evaluated twice per time step (every 2 fs),
    and force group 2 should be evaluated eight times per time step (every 0.5 fs).

    A common use of this algorithm is to evaluate reciprocal space nonbonded interactions
    less often than the bonded and direct space nonbonded interactions.  The following
    example looks up the NonbondedForce, sets the reciprocal space interactions to their
    own force group, and then creates an integrator that evaluates them once every 4 fs,
    but all other interactions every 2 fs::

        nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
        nonbonded.setReciprocalSpaceForceGroup(1)
        integrator = MTSIntegrator(4*femtoseconds, [(1,1), (0,2)])

    For details, see Tuckerman et al., J. Chem. Phys. 97(3) pp. 1990-2001 (1992).
    """

    def __init__(self, dt, groups):
        """Create an MTSIntegrator.

        Parameters
        ----------
        dt : time
            The largest (outermost) integration time step to use
        groups : list
            A list of tuples defining the force groups.  The first element of
            each tuple is the force group index, and the second element is the
            number of times that force group should be evaluated in one time step.
        """
        if len(groups) == 0:
            raise ValueError("No force groups specified")
        groups = sorted(groups, key=lambda x: x[1])
        CustomIntegrator.__init__(self, dt)
        self.addPerDofVariable("x1", 0)
        self.addUpdateContextState();
        self._createSubsteps(1, groups)
        self.addConstrainVelocities();

    def _createSubsteps(self, parentSubsteps, groups):
        group, substeps = groups[0]
        stepsPerParentStep = substeps / parentSubsteps
        if stepsPerParentStep < 1 or stepsPerParentStep != int(stepsPerParentStep):
            raise ValueError("The number for substeps for each group must be a multiple of the number for the previous group")
        stepsPerParentStep = int(stepsPerParentStep)
        if group < 0 or group > 31:
            raise ValueError("Force group must be between 0 and 31")
        for i in range(stepsPerParentStep):
            self.addComputePerDof("v", "v+0.5*(dt/"+str(substeps)+")*f"+str(group)+"/m")
            if len(groups) == 1:
                self.addComputePerDof("x", "x+(dt/"+str(substeps)+")*v")
                self.addComputePerDof("x1", "x")
                self.addConstrainPositions();
                self.addComputePerDof("v", "v+(x-x1)/(dt/"+str(substeps)+")");
                self.addConstrainVelocities()
            else:
                self._createSubsteps(substeps, groups[1:])
            self.addComputePerDof("v", "v+0.5*(dt/"+str(substeps)+")*f"+str(group)+"/m")

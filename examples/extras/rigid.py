"""
rigid.py: Implements rigid bodies

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970.

Portions copyright (c) 2016-2025 Stanford University and the Authors.
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

import openmm as mm
import openmm.unit as unit
import numpy as np
import numpy.linalg as lin
from itertools import combinations

def createRigidBodies(system, positions, bodies):
    """Modify a System to turn specified sets of particles into rigid bodies.
    
    For every rigid body, four particles are selected as "real" particles whose positions are integrated.
    Constraints are added between them to make them move as a rigid body.  All other particles in the body
    are then turned into virtual sites whose positions are computed based on the "real" particles.
    
    Because virtual sites are massless, the mass properties of the rigid bodies will be slightly different
    from the corresponding sets of particles in the original system.  The masses of the non-virtual particles
    are chosen to guarantee that the total mass and center of mass of each rigid body exactly match those of
    the original particles.  The moment of inertia will be similar to that of the original particles, but
    not identical.
    
    Care is needed when using constraints, since virtual particles cannot participate in constraints.  If the
    input system includes any constraints, this function will automatically remove ones that connect two
    particles in the same rigid body.  But if there is a constraint beween a particle in a rigid body and
    another particle not in that body, it will likely lead to an exception when you try to create a context.
    
    Parameters:
     - system (System) the System to modify
     - positions (list) the positions of all particles in the system
     - bodies (list) each element of this list defines one rigid body.  Each element should itself be a list
       of the indices of all particles that make up that rigid body.
    """
    # Remove any constraints involving particles in rigid bodies.
    
    for i in range(system.getNumConstraints()-1, -1, -1):
        p1, p2, distance = system.getConstraintParameters(i)
        if (any(p1 in body and p2 in body for body in bodies)):
            system.removeConstraint(i)
    
    # Loop over rigid bodies and process them.
    
    for particles in bodies:
        if len(particles) < 5:
            # All the particles will be "real" particles.
            
            realParticles = particles
            realParticleMasses = [system.getParticleMass(i) for i in particles]
        else:
            # Select four particles to use as the "real" particles.  All others will be virtual sites.
            
            pos = [positions[i] for i in particles]
            mass = [system.getParticleMass(i) for i in particles]
            cm = unit.sum([p*m for p, m in zip(pos, mass)])/unit.sum(mass)
            r = [p-cm for p in pos]
            avgR = unit.sqrt(unit.sum([unit.dot(x, x) for x in r])/len(particles))
            rank = sorted(range(len(particles)), key=lambda i: abs(unit.norm(r[i])-avgR))
            for p in combinations(rank, 4):
                # Select masses for the "real" particles.  If any is negative, reject this set of particles
                # and keep going.
                
                matrix = np.zeros((4, 4))
                for i in range(4):
                    particleR = r[p[i]].value_in_unit(unit.nanometers)
                    matrix[0][i] = particleR[0]
                    matrix[1][i] = particleR[1]
                    matrix[2][i] = particleR[2]
                    matrix[3][i] = 1.0
                rhs = np.array([0.0, 0.0, 0.0, unit.sum(mass).value_in_unit(unit.amu)])
                weights = lin.solve(matrix, rhs)
                if all(w > 0.0 for w in weights):
                    # We have a good set of particles.
                    
                    realParticles = [particles[i] for i in p]
                    realParticleMasses = [float(w) for w in weights]*unit.amu
                    break
        
        # Set particle masses.
        
        for i, m in zip(realParticles, realParticleMasses):
            system.setParticleMass(i, m)
        
        # Add constraints between the real particles.
        
        for p1, p2 in combinations(realParticles, 2):
            distance = unit.norm(positions[p1]-positions[p2])
            key = (min(p1, p2), max(p1, p2))
            system.addConstraint(p1, p2, distance)
        
        # Select which three particles to use for defining virtual sites.
        
        bestNorm = 0
        for p1, p2, p3 in combinations(realParticles, 3):
            d12 = (positions[p2]-positions[p1]).value_in_unit(unit.nanometer)
            d13 = (positions[p3]-positions[p1]).value_in_unit(unit.nanometer)
            crossNorm = unit.norm((d12[1]*d13[2]-d12[2]*d13[1], d12[2]*d13[0]-d12[0]*d13[2], d12[0]*d13[1]-d12[1]*d13[0]))
            if crossNorm > bestNorm:
                bestNorm = crossNorm
                vsiteParticles = (p1, p2, p3)
        
        # Create virtual sites.
        
        d12 = (positions[vsiteParticles[1]]-positions[vsiteParticles[0]]).value_in_unit(unit.nanometer)
        d13 = (positions[vsiteParticles[2]]-positions[vsiteParticles[0]]).value_in_unit(unit.nanometer)
        cross = mm.Vec3(d12[1]*d13[2]-d12[2]*d13[1], d12[2]*d13[0]-d12[0]*d13[2], d12[0]*d13[1]-d12[1]*d13[0])
        matrix = np.zeros((3, 3))
        for i in range(3):
            matrix[i][0] = d12[i]
            matrix[i][1] = d13[i]
            matrix[i][2] = cross[i]
        for i in particles:
            if i not in realParticles:
                system.setParticleMass(i, 0)
                rhs = np.array((positions[i]-positions[vsiteParticles[0]]).value_in_unit(unit.nanometer))
                weights = lin.solve(matrix, rhs)
                system.setVirtualSite(i, mm.OutOfPlaneSite(vsiteParticles[0], vsiteParticles[1], vsiteParticles[2], weights[0], weights[1], weights[2]))

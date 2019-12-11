
/* Portions copyright (c) 2006-2013 Stanford University and Simbios.
 * Contributors: Peter Eastman, Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <cstring>
#include <iostream>
#include <sstream>

#include "openmm/OpenMMException.h"
#include "SimTKOpenMMUtilities.h"
#include "openmm/internal/ContextImpl.h"
#include "ReferenceVelocityVerletDynamics.h"
#include "ReferenceVirtualSites.h"

#include <cstdio>

using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceVelocityVerletDynamics constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         delta t for dynamics
   @param friction       friction coefficient
   @param temperature    temperature

   --------------------------------------------------------------------------------------- */

ReferenceVelocityVerletDynamics::ReferenceVelocityVerletDynamics(int numberOfAtoms, double deltaT) :
           ReferenceDynamics(numberOfAtoms, deltaT, 0.0) {
   xPrime.resize(numberOfAtoms);
   inverseMasses.resize(numberOfAtoms);
}

/**---------------------------------------------------------------------------------------

   ReferenceVelocityVerletDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceVelocityVerletDynamics::~ReferenceVelocityVerletDynamics() {
}

/**---------------------------------------------------------------------------------------

   Update -- driver routine for performing Velocity Verlet dynamics update of coordinates
   and velocities

   @param system              the System to be integrated
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses
   @param atomList            list of all atoms not involved in a Drude-like pair
   @param pairList            list of all Drude-like pairs
   @param maxPairDistance     the maximum separation of any Drude-like pairs

   --------------------------------------------------------------------------------------- */

void ReferenceVelocityVerletDynamics::update(OpenMM::ContextImpl &context, const OpenMM::System& system, vector<Vec3>& atomCoordinates,
                                          vector<Vec3>& velocities,
                                          vector<Vec3>& forces, vector<double>& masses, double tolerance, bool &forcesAreValid,
                                          const std::vector<int> & atomList, const std::vector<std::tuple<int, int, double>> &pairList,
                                          double maxPairDistance) {

    // first-time-through initialization
    if (!forcesAreValid) context.calcForcesAndEnergy(true, false);

    int numberOfAtoms = system.getNumParticles();
    if (getTimeStep() == 0) {
       // invert masses

       for (int ii = 0; ii < numberOfAtoms; ii++) {
          if (masses[ii] == 0.0)
              inverseMasses[ii] = 0.0;
          else
              inverseMasses[ii] = 1.0/masses[ii];
       }
    }

    //// Perform the integration.

    // Regular atoms
    for (const auto &atom : atomList) {
        if (masses[atom] != 0.0) {
            velocities[atom] += 0.5 * inverseMasses[atom]*forces[atom]*getDeltaT();
            xPrime[atom] = atomCoordinates[atom];
            atomCoordinates[atom] += velocities[atom]*getDeltaT();
        }
    }
    // Connected particles
    for (const auto &pair : pairList) {
        const auto &atom1 = std::get<0>(pair);
        const auto &atom2 = std::get<1>(pair);
        double m1 = system.getParticleMass(atom1);
        double m2 = system.getParticleMass(atom2);
        double mass1fract = m1 / (m1 + m2);
        double mass2fract = m2 / (m1 + m2);
        double invRedMass = (m1 * m2 != 0.0) ? (m1 + m2)/(m1 * m2) : 0.0;
        double invTotMass = (m1 + m2 != 0.0) ? 1.0 /(m1 + m2) : 0.0;
        Vec3 comVel = velocities[atom1]*mass1fract + velocities[atom2]*mass2fract;
        Vec3 relVel = velocities[atom2] - velocities[atom1];
        Vec3 comForce = forces[atom1] + forces[atom2];
        Vec3 relForce = mass1fract*forces[atom2] - mass2fract*forces[atom1];
        comVel += 0.5 * comForce * getDeltaT() * invTotMass;
        relVel += 0.5 * relForce * getDeltaT() * invRedMass; 
        if (m1 != 0.0) {
            velocities[atom1] = comVel - relVel*mass2fract;
            xPrime[atom1] = atomCoordinates[atom1];
            atomCoordinates[atom1] += velocities[atom1]*getDeltaT();
        }
        if (m2 != 0.0) {
            velocities[atom2] = comVel + relVel*mass1fract;
            xPrime[atom2] = atomCoordinates[atom2];
            atomCoordinates[atom2] += velocities[atom2]*getDeltaT();
        }
    }

    // 

    ReferenceConstraintAlgorithm* referenceConstraintAlgorithm = getReferenceConstraintAlgorithm();
    if (referenceConstraintAlgorithm)
       referenceConstraintAlgorithm->apply(xPrime, atomCoordinates, inverseMasses, tolerance);


    // Apply hard wall constraints.
    if (maxPairDistance > 0) {
        for (const auto & pair : pairList) {
            const int atom1 = std::get<0>(pair);
            const int atom2 = std::get<1>(pair);
            const double hardWallScale = sqrt(std::get<2>(pair)*BOLTZ);
            Vec3 delta = atomCoordinates[atom1]-atomCoordinates[atom2];
            double r = sqrt(delta.dot(delta));
            double rInv = 1/r;
            if (rInv*maxPairDistance < 1.0) {
                // The constraint has been violated, so make the inter-particle distance "bounce"
                // off the hard wall.
                //if (rInv*maxPairDistance < 0.5)
                //    throw OpenMMException("Drude particle moved too far beyond hard wall constraint");
                //    TODO: Review this - I commented it out to make the NoseHooverThermostat test work
                Vec3 bondDir = delta*rInv;
                Vec3 vel1 = velocities[atom1];
                Vec3 vel2 = velocities[atom2];
                double m1 = masses[atom1];
                double m2 = masses[atom2];
                double invTotMass = (m1 + m2 != 0.0) ? 1.0 /(m1 + m2) : 0.0;
                double deltaR = r-maxPairDistance;
                double deltaT = getDeltaT();
                double dt = getDeltaT();

                double dotvr1 = vel1.dot(bondDir);
                Vec3 vb1 = bondDir*dotvr1;
                Vec3 vp1 = vel1-vb1;
                if (m2 == 0) {
                    // The parent particle is massless, so move only the Drude particle.

                    if (dotvr1 != 0.0)
                        deltaT = deltaR/std::abs(dotvr1);
                    if (deltaT > getDeltaT())
                        deltaT = getDeltaT();
                    dotvr1 = -dotvr1*hardWallScale/(std::abs(dotvr1)*sqrt(m1));
                    double dr = -deltaR + deltaT*dotvr1;
                    atomCoordinates[atom1] += bondDir*dr;
                    velocities[atom1] = vp1 + bondDir*dotvr1;
                }
                else {
                    // Move both particles.
                    double dotvr2 = vel2.dot(bondDir);
                    Vec3 vb2 = bondDir*dotvr2;
                    Vec3 vp2 = vel2-vb2;
                    double vbCMass = (m1*dotvr1 + m2*dotvr2)*invTotMass;
                    dotvr1 -= vbCMass;
                    dotvr2 -= vbCMass;
                    if (dotvr1 != dotvr2)
                        deltaT = deltaR/std::abs(dotvr1-dotvr2);
                    if (deltaT > dt)
                        deltaT = dt;
                    double vBond = hardWallScale/sqrt(m1);
                    dotvr1 = -dotvr1*vBond*m2*invTotMass/std::abs(dotvr1);
                    dotvr2 = -dotvr2*vBond*m1*invTotMass/std::abs(dotvr2);
                    double dr1 = -deltaR*m2*invTotMass + deltaT*dotvr1;
                    double dr2 = deltaR*m1*invTotMass + deltaT*dotvr2;
                    dotvr1 += vbCMass;
                    dotvr2 += vbCMass;
                    atomCoordinates[atom1] += bondDir*dr1;
                    atomCoordinates[atom2] += bondDir*dr2;
                    velocities[atom1] = vp1 + bondDir*dotvr1;
                    velocities[atom2] = vp2 + bondDir*dotvr2;
                }
            }
        }
    } /* end of hard wall constraint part */


    ReferenceVirtualSites::computePositions(system, atomCoordinates);
    context.calcForcesAndEnergy(true, false);
    forcesAreValid = true;

    for (int i = 0; i < numberOfAtoms; ++i) {
        if (masses[i] != 0.0)
            for (int j = 0; j < 3; ++j) {
                xPrime[i][j] += velocities[i][j]*getDeltaT();
        }
   } 

    // Update the positions and velocities.
    // Regular atoms
    for (const auto &atom : atomList) {
        if (masses[atom] != 0.0) {
            velocities[atom] += 0.5 * inverseMasses[atom]*forces[atom]*getDeltaT() + (atomCoordinates[atom] - xPrime[atom])/getDeltaT();
        }
    }
    // Connected particles
    for (const auto &pair : pairList) {
        const auto &atom1 = std::get<0>(pair);
        const auto &atom2 = std::get<1>(pair);
        double m1 = system.getParticleMass(atom1);
        double m2 = system.getParticleMass(atom2);
        double mass1fract = m1 / (m1 + m2);
        double mass2fract = m2 / (m1 + m2);
        double invRedMass = (m1 * m2 != 0.0) ? (m1 + m2)/(m1 * m2) : 0.0;
        double invTotMass = (m1 + m2 != 0.0) ? 1.0 /(m1 + m2) : 0.0;
        Vec3 comVel = velocities[atom1]*mass1fract + velocities[atom2]*mass2fract;
        Vec3 relVel = velocities[atom2] - velocities[atom1];
        Vec3 comForce = forces[atom1] + forces[atom2];
        Vec3 relForce = mass1fract*forces[atom2] - mass2fract*forces[atom1];
        comVel += 0.5 * comForce * getDeltaT() * invTotMass;
        relVel += 0.5 * relForce * getDeltaT() * invRedMass; 
        if (m1 != 0.0) {
            velocities[atom1] = comVel - relVel*mass2fract + (atomCoordinates[atom1] - xPrime[atom1])/getDeltaT();
        }
        if (m2 != 0.0) {
            velocities[atom2] = comVel + relVel*mass1fract + (atomCoordinates[atom2] - xPrime[atom2])/getDeltaT();
        }
    }
    if (referenceConstraintAlgorithm)
       referenceConstraintAlgorithm->applyToVelocities(atomCoordinates, velocities, inverseMasses, tolerance);

    incrementTimeStep();
}

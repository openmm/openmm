/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2023 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "sfmt/SFMT.h"
#include "SimTKOpenMMRealType.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/DrudeForce.h"
#include "openmm/DrudeLangevinIntegrator.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include "openmm/OpenMMException.h"

#include <set>

using namespace std;

namespace OpenMM {

/**
 * Identify normal particles (not part of a pair) and Drude particle pairs.
 */
void findParticlesAndPairs(const System &system, vector<int>& normalParticles, vector<pair<int, int> >& pairParticles) {
    // Find the underlying Drude force object
    const DrudeForce* drudeForce = NULL;
    for (int i = 0; i < system.getNumForces(); i++)
        if (dynamic_cast<const DrudeForce*>(&system.getForce(i)) != NULL) {
            if (drudeForce == NULL)
                drudeForce = dynamic_cast<const DrudeForce*>(&system.getForce(i));
            else
                throw OpenMMException("The System contains multiple DrudeForces");
        }
    if (drudeForce == NULL)
        throw OpenMMException("The System does not contain a DrudeForce");

    // Figure out which particles are individual and which are Drude pairs
    set<int> particles;
    for (int i = 0; i < system.getNumParticles(); i++)
        particles.insert(i);
    for (int i = 0; i < drudeForce->getNumParticles(); i++) {
        int p, p1, p2, p3, p4;
        double charge, polarizability, aniso12, aniso34;
        drudeForce->getParticleParameters(i, p, p1, p2, p3, p4, charge, polarizability, aniso12, aniso34);
        particles.erase(p);
        particles.erase(p1);
        pairParticles.emplace_back(p, p1);
    }
    normalParticles.insert(normalParticles.begin(), particles.begin(), particles.end());
}

vector<Vec3> assignDrudeVelocities(const System &system, double temperature, double drudeTemperature, int randomSeed) {
    vector<int> normalParticles;
    vector<pair<int, int> > pairParticles;
    findParticlesAndPairs(system, normalParticles, pairParticles);

    // Generate the list of Gaussian random numbers.
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(randomSeed, sfmt);
    vector<double> randoms;
    while (randoms.size() < system.getNumParticles()*3) {
        double x, y, r2;
        do {
            x = 2.0*genrand_real2(sfmt)-1.0;
            y = 2.0*genrand_real2(sfmt)-1.0;
            r2 = x*x + y*y;
        } while (r2 >= 1.0 || r2 == 0.0);
        double multiplier = sqrt((-2.0*log(r2))/r2);
        randoms.push_back(x*multiplier);
        randoms.push_back(y*multiplier);
    }

    // Assign the velocities.
    vector<Vec3> velocities(system.getNumParticles(), Vec3());
    int nextRandom = 0;
    // First the indivitual atoms
    for (const auto &atom : normalParticles ) {
        double mass = system.getParticleMass(atom);
        if (mass != 0) {
            double velocityScale = sqrt(BOLTZ*temperature/mass);
            velocities[atom] = Vec3(randoms[nextRandom++], randoms[nextRandom++], randoms[nextRandom++])*velocityScale;
        }
    }
    // Now the particle-Drude pairs
    for (const auto &pair : pairParticles ) {
        const auto atom1 = pair.first;
        const auto atom2 = pair.second;
        double mass1 = system.getParticleMass(atom1);
        double mass2 = system.getParticleMass(atom2);
        if (mass1 != 0 && mass2 != 0) {
            double invMass = 1.0 / (mass1 + mass2);
            double redMass = mass1 * mass2 * invMass;
            double fracM1 = mass1 * invMass;
            double fracM2 = mass2 * invMass;
            Vec3 comVelocity = Vec3(randoms[nextRandom++], randoms[nextRandom++], randoms[nextRandom++])*sqrt(BOLTZ*temperature*invMass);
            Vec3 relVelocity = Vec3(randoms[nextRandom++], randoms[nextRandom++], randoms[nextRandom++])*sqrt(BOLTZ*drudeTemperature/redMass);
            velocities[atom1] = comVelocity - fracM2 * relVelocity;
            velocities[atom2] = comVelocity + fracM1 * relVelocity;
        }
        else if (mass2 != 0) {
            double velocityScale = sqrt(BOLTZ*temperature/mass2);
            velocities[atom2] = Vec3(randoms[nextRandom++], randoms[nextRandom++], randoms[nextRandom++])*velocityScale;
        }
    }
    return velocities;
}


/**
 * Computes the instantaneous temperatures of the system and the internal Drude motion and returns a pair (T_system, T_drude)
 */
pair<double, double> computeTemperaturesFromVelocities(const System& system, const vector<Vec3>& velocities) {
    vector<int> normalParticles;
    vector<pair<int, int> > pairParticles;
    findParticlesAndPairs(system, normalParticles, pairParticles);
    double energy = 0.0, drudeEnergy = 0;
    int dof = 0, drudeDof = 0;
    
    // Kinetic energy and degrees of freedom from normal particles.
    
    for (int i : normalParticles) {
        double mass = system.getParticleMass(i);
        if (mass > 0) {
            energy += mass*(velocities[i].dot(velocities[i]));
            dof += 3;
        }
    }
    
    // Kinetic energy and degrees of freedom from Drude particle pairs.
    
    for (auto pair : pairParticles) {
        int p1 = pair.first;
        int p2 = pair.second;
        double mass1 = system.getParticleMass(p1);
        double mass2 = system.getParticleMass(p2);
	double reducedMass = (mass1*mass2)/(mass1+mass2);
        if (mass1 != 0 || mass2 != 0) {
            Vec3 momentum = mass1*velocities[p1] + mass2*velocities[p2];
            energy += momentum.dot(momentum)/(mass1+mass2);
	    Vec3 drudeVelocity = velocities[p1] - velocities[p2];
            drudeEnergy += reducedMass*drudeVelocity.dot(drudeVelocity);
	    dof += 3;
	    drudeDof += 3;
        }
    }
    
    // Reduce degrees of freedom for constraints.

    for (int i = 0; i < system.getNumConstraints(); i++) {
        int p1, p2;
        double distance;
        system.getConstraintParameters(i, p1, p2, distance);
        if (system.getParticleMass(p1) > 0 || system.getParticleMass(p2) > 0)
            dof--;
    }

    // Reduce degrees of freedom if there is a CMMotionRemover.

    for (int i = 0; i < system.getNumForces(); i++)
        if (dynamic_cast<const CMMotionRemover*>(&system.getForce(i)) != NULL) {
            dof -= 3;
            break;
        }
    energy *= 0.5;
    drudeEnergy *= 0.5;
    if (drudeDof == 0)
        drudeDof = 1; // so that the drude temperature is reported as 0
    return make_pair<double, double>(2*energy/(dof*BOLTZ), 2*drudeEnergy/(drudeDof*BOLTZ));
}

}

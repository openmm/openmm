/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/Context.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/AndersenThermostat.h"
#include "openmm/MonteCarloBarostat.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>
#include <stdlib.h>

using namespace OpenMM;
using namespace std;

void testSerialization() {

    // Create a System.
	const int numParticles=50;
    System system;
	system.setDefaultPeriodicBoxVectors(Vec3(6.2, 0, 0), Vec3(0, 6.2, 0), Vec3(0, 0, 6.2 ));
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::Ewald);
    nonbonded->setCutoffDistance(0.8);
    nonbonded->setEwaldErrorTolerance(0.01);
    for (int i = 0; i < numParticles/2; i++)
        system.addParticle(22.99);
    for (int i = 0; i < numParticles/2; i++)
        system.addParticle(35.45);
    for (int i = 0; i < numParticles/2; i++)
        nonbonded->addParticle(1.0, 1.0,0.0);
    for (int i = 0; i < numParticles/2; i++)
        nonbonded->addParticle(-1.0, 1.0,0.0);
    system.addForce(nonbonded);
	system.addForce(new AndersenThermostat(393.3, 19.3));
	system.addForce(new MonteCarloBarostat(25, 393.3, 25));
	LangevinIntegrator intg(300,79,0.002);
	Context context(system, intg);
	
	// Set positions, velocities, forces
	vector<Vec3> positions;
	for (int i = 0; i < numParticles; i++) {
		positions.push_back(Vec3( ((float) rand()/(float) RAND_MAX)*6.2, ((float) rand()/(float) RAND_MAX)*6.2, ((float) rand()/(float) RAND_MAX)*6.2));
	}
	vector<Vec3> velocities;
	for (int i = 0; i < numParticles; i++) {
		velocities.push_back(Vec3( ((float) rand()/(float) RAND_MAX)*6.2, ((float) rand()/(float) RAND_MAX)*6.2, ((float) rand()/(float) RAND_MAX)*6.2));
	}

	context.setPositions(positions);
	context.setVelocities(velocities);

	// Serialize and then deserialize it.
	State s1 = context.getState(State::Positions | State::Velocities | State::Forces | State::Energy | State::Parameters);

	stringstream buffer;
    XmlSerializer::serialize<State>(&s1, "State", buffer);
    State* copy = XmlSerializer::deserialize<State>(buffer);
	State& s2 = *copy;

    // Compare the two states to see if they are identical.
	vector<Vec3> pos1 = s1.getPositions();
	vector<Vec3> pos2 = s2.getPositions();
	ASSERT_EQUAL(pos1.size(), pos2.size());
	ASSERT_EQUAL(pos1.size(), positions.size());
	for (int i = 0; i < (int) pos1.size(); i++) {
		ASSERT_EQUAL_VEC(pos1[i],pos2[i],0);
	}
	vector<Vec3> vel1 = s1.getVelocities();
	vector<Vec3> vel2 = s2.getVelocities();
	ASSERT_EQUAL(vel1.size(), vel2.size());
	for (int i = 0; i < (int) pos1.size(); i++) {
		ASSERT_EQUAL_VEC(vel1[i],vel2[i],0);
	}
	vector<Vec3> forces1 = s1.getForces();
	vector<Vec3> forces2 = s2.getForces();
	ASSERT_EQUAL(forces1.size(), forces2.size());
	for (int i = 0; i < (int) pos1.size(); i++) {
		ASSERT_EQUAL_VEC(forces1[i],forces2[i],0);
	}
	Vec3 a1,a2,a3,b1,b2,b3;
	s1.getPeriodicBoxVectors(a1,a2,a3);
	s2.getPeriodicBoxVectors(b1,b2,b3);
	ASSERT_EQUAL_VEC(a1,b1,0);
	ASSERT_EQUAL_VEC(a2,b2,0);
	ASSERT_EQUAL_VEC(a3,b3,0);

	ASSERT_EQUAL(s1.getPotentialEnergy(), s2.getPotentialEnergy());
	ASSERT_EQUAL(s1.getKineticEnergy(), s2.getKineticEnergy());
	ASSERT_EQUAL(s1.getTime(), s2.getTime());

	map<string, double> p1 = s1.getParameters();
	map<string, double> p2 = s2.getParameters();

	ASSERT_EQUAL(p1.size(), p2.size());
	map<string, double>::const_iterator it1=p1.begin();
	map<string, double>::const_iterator it2=p2.begin();
	//maps are ordered, so iterators should be in the same order. 
	for (it1 = p1.begin(); it1 != p1.end(); it1++, it2++) {
		assert((it1->first).compare(it2->first) == 0);
		ASSERT_EQUAL(it1->second, it2->second);
	}
    delete copy;

    // Now create a series of States that include only one type of information.  Verify
    // that serialization works correctly for them.

    for (int types = 1; types <= 16; types *= 2) {
        State s3 = context.getState(types);
        stringstream buffer2;
        XmlSerializer::serialize<State>(&s3, "State", buffer2);
        copy = XmlSerializer::deserialize<State>(buffer2);
        int foundTypes = 0;
        try {
            copy->getPositions();
            foundTypes += State::Positions;
        }
        catch (...) {
            // Ignore
        }
        try {
            copy->getVelocities();
            foundTypes += State::Velocities;
        }
        catch (...) {
            // Ignore
        }
        try {
            copy->getForces();
            foundTypes += State::Forces;
        }
        catch (...) {
            // Ignore
        }
        try {
            copy->getPotentialEnergy();
            foundTypes += State::Energy;
        }
        catch (...) {
            // Ignore
        }
        try {
            copy->getParameters();
            foundTypes += State::Parameters;
        }
        catch (...) {
            // Ignore
        }
        delete copy;
        ASSERT_EQUAL(types, foundTypes);
    }
}

int main() {
    try {
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;  
		return 1;
    }
    cout << "Done" << endl;
    return 0;
}



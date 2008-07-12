/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

/**
 * This tests the findExclusions() method of StandardMMForceFieldImpl, which identifies pairs of atoms
 * whose nonbonded atoms are either excluded or decreased.  The test system is a chain with branches:
 * 
 * 1  3  5  7  9  11  13  15  17  19
 * |  |  |  |  |  |   |   |   |   |
 * 0--2--4--6--8--10--12--14--16--18
 */

#include "AssertionUtilities.h"
#include "Kernel.h"
#include "KernelFactory.h"
#include "OpenMMContext.h"
#include "Platform.h"
#include "StandardMMForceField.h"
#include "Stream.h"
#include "StreamFactory.h"
#include "System.h"
#include "VerletIntegrator.h"
#include "kernels.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <set>
#include <vector>

using namespace OpenMM;
using namespace std;

static const int NUM_ATOMS = 20;

/**
 * Add a pair of atoms to the list of exclusions.
 */

void addAtomsToExclusions(int atom1, int atom2, vector<set<int> >& exclusions) {
    if (atom2 < NUM_ATOMS) {
        exclusions[atom1].insert(atom2);
        exclusions[atom2].insert(atom1);
    }
}

/**
 * Verify that the exclusions are what we expect.
 */

void verifyExclusions(const vector<set<int> >& exclusions) {
    vector<set<int> > expected(NUM_ATOMS);
    for (int i = 0; i < NUM_ATOMS; i += 2) {
        addAtomsToExclusions(i, i+1, expected);
        addAtomsToExclusions(i, i+2, expected);
        addAtomsToExclusions(i, i+3, expected);
        addAtomsToExclusions(i, i+4, expected);
        addAtomsToExclusions(i, i+5, expected);
        addAtomsToExclusions(i, i+6, expected);
        addAtomsToExclusions(i+1, i+2, expected);
        addAtomsToExclusions(i+1, i+3, expected);
        addAtomsToExclusions(i+1, i+4, expected);
    }
    ASSERT_EQUAL(expected.size(), exclusions.size());
    for (int i = 0; i < NUM_ATOMS; ++i) {
        ASSERT_EQUAL(expected[i].size(), exclusions[i].size());
        vector<int> intersection(0);
        insert_iterator<vector<int> > inserter(intersection, intersection.begin());
        set_intersection(exclusions[i].begin(), exclusions[i].end(), expected[i].begin(), expected[i].end(), inserter);
        ASSERT_EQUAL(expected[i].size(), intersection.size());
    }
}

/**
 * Add a pair of atoms to the list of 1-4 pairs.
 */

void addAtomsTo14List(int atom1, int atom2, set<pair<int, int> >& bonded14Indices) {
    if (atom2 < NUM_ATOMS)
        bonded14Indices.insert(pair<int, int>(atom1, atom2));
}

/**
 * Verify that the 1-4 pairs are what we expect.
 */

void verify14(const vector<vector<int> >& bonded14Indices) {
    set<pair<int, int> > expected, found;
    for (int i = 0; i < NUM_ATOMS; i += 2) {
        addAtomsTo14List(i, i+5, expected);
        addAtomsTo14List(i, i+6, expected);
        addAtomsTo14List(i+1, i+3, expected);
        addAtomsTo14List(i+1, i+4, expected);
    }
    ASSERT_EQUAL(expected.size(), bonded14Indices.size());
    for (size_t i = 0; i < bonded14Indices.size(); ++i) {
        int atom1 = bonded14Indices[i][0];
        int atom2 = bonded14Indices[i][1];
        found.insert(pair<int, int>(min(atom1, atom2), max(atom1, atom2)));
    }
    vector<pair<int, int> > intersection(0);
    insert_iterator<vector<pair<int, int> > > inserter(intersection, intersection.begin());
    set_intersection(expected.begin(), expected.end(), found.begin(), found.end(), inserter);
    ASSERT_EQUAL(expected.size(), intersection.size());
}

/**
 * The following classes define a Platform whose job is to check whether the correct values were passed
 * to the initialize() methods.
 */

class DummyForceKernel : public CalcStandardMMForceFieldKernel {
public:
    DummyForceKernel(string name, const Platform& platform) : CalcStandardMMForceFieldKernel(name, platform) {
    }
    void initialize(const vector<vector<int> >& bondIndices, const vector<vector<double> >& bondParameters,
            const vector<vector<int> >& angleIndices, const vector<vector<double> >& angleParameters,
            const vector<vector<int> >& periodicTorsionIndices, const vector<vector<double> >& periodicTorsionParameters,
            const vector<vector<int> >& rbTorsionIndices, const vector<vector<double> >& rbTorsionParameters,
            const vector<vector<int> >& bonded14Indices, double lj14Scale, double coulomb14Scale,
            const vector<set<int> >& exclusions, const vector<vector<double> >& nonbondedParameters,
            NonbondedMethod nonbondedMethod, double nonbondedCutoff, double periodicBoxSize[3]) {
        verifyExclusions(exclusions);
        verify14(bonded14Indices);
    }
    void executeForces(const Stream& positions, Stream& forces) {
    }
    double executeEnergy(const Stream& positions) {
		return 0.0;
    }
};

class DummyIntegratorKernel : public IntegrateVerletStepKernel {
public:
    DummyIntegratorKernel(string name, const Platform& platform) : IntegrateVerletStepKernel(name, platform) {
    }
    void initialize(const vector<double>& masses, const vector<vector<int> >& constraintIndices, const vector<double>& constraintLengths) {
    }
    void execute(Stream& positions, Stream& velocities, const Stream& forces, double stepSize) {
    }
};

class DummyKEKernel : public CalcKineticEnergyKernel {
public:
    DummyKEKernel(string name, const Platform& platform) : CalcKineticEnergyKernel(name, platform) {
    }
    void initialize(const vector<double>& masses) {
    }
    double execute(const Stream& positions) {
        return 0.0;
    }
};

class DummyStreamImpl : public StreamImpl {
public:
    DummyStreamImpl(string name, int size, Stream::DataType type, const Platform& platform) : StreamImpl(name, size, type, platform) {
    }
    void loadFromArray(const void* array) {
    }
    void saveToArray(void* array) {
    }
    void fillWithValue(void* value) {
    }
};

class DummyKernelFactory : public KernelFactory {
public:
    KernelImpl* createKernelImpl(string name, const Platform& platform, OpenMMContextImpl& context) const {
        if (name == CalcStandardMMForceFieldKernel::Name())
            return new DummyForceKernel(name, platform);
        if (name == IntegrateVerletStepKernel::Name())
            return new DummyIntegratorKernel(name, platform);
        if (name == CalcKineticEnergyKernel::Name())
            return new DummyKEKernel(name, platform);
        return 0;
    }
};

class DummyStreamFactory : public StreamFactory {
public:
    StreamImpl* createStreamImpl(string name, int size, Stream::DataType type, const Platform& platform, OpenMMContextImpl& context) const {
        return new DummyStreamImpl(name, size, type, platform);
    }
};

class DummyPlatform : public Platform {
public:
    DummyPlatform() {
        registerKernelFactory(CalcStandardMMForceFieldKernel::Name(), new DummyKernelFactory());
        registerKernelFactory(IntegrateVerletStepKernel::Name(), new DummyKernelFactory());
        registerKernelFactory(CalcKineticEnergyKernel::Name(), new DummyKernelFactory());
    }
    string getName() const {
        return "Dummy";
    }
    double getSpeed() const {
        return 1.0;
    }
    bool supportsDoublePrecision() const {
        return true;
    }
    const StreamFactory& getDefaultStreamFactory() const {
        return streamFactory;
    }
private:
    DummyStreamFactory streamFactory;
};

int main() {
    try {
        DummyPlatform platform;
        System system(NUM_ATOMS, 0);
        VerletIntegrator integrator(0.01);
        StandardMMForceField* forces = new StandardMMForceField(NUM_ATOMS, NUM_ATOMS-1, 0, 0, 0);
        
        // loop over all main-chain atoms (even numbered atoms)
        for (int i = 0; i < NUM_ATOMS-1; i += 2) 
        {
        	// side-chain bonds
        	forces->setBondParameters(i, i, i+1, 1.0, 1.0);

        	// main-chain bonds
            if (i < NUM_ATOMS-2) // penultimate atom (NUM_ATOMS-2) has no subsequent main-chain atom
                forces->setBondParameters(i+1, i, i+2, 1.0, 1.0);
        }
        
        system.addForce(forces);
        OpenMMContext context(system, integrator, platform);
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

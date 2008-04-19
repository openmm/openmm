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

#include "ReferenceKernels.h"

using namespace OpenMM;
using namespace std;

int** allocateIntArray(int length, int width) {
    int** array = new int*[length];
    for (int i = 0; i < length; ++i)
        array[i] = new int[width];
    return array;
}

RealOpenMM** allocateRealArray(int length, int width) {
    RealOpenMM** array = new RealOpenMM*[length];
    for (int i = 0; i < length; ++i)
        array[i] = new RealOpenMM[width];
    return array;
}

int** copyToArray(const vector<vector<int> > vec) {
    if (vec.size() == 0)
        return new int*[0];
    int** array = allocateIntArray(vec.size(), vec[0].size());
    for (int i = 0; i < vec.size(); ++i)
        for (int j = 0; j < vec[i].size(); ++j)
            array[i][j] = vec[i][j];
    return array;
}

RealOpenMM** copyToArray(const vector<vector<double> > vec) {
    if (vec.size() == 0)
        return new RealOpenMM*[0];
    RealOpenMM** array = allocateRealArray(vec.size(), vec[0].size());
    for (int i = 0; i < vec.size(); ++i)
        for (int j = 0; j < vec[i].size(); ++j)
            array[i][j] = vec[i][j];
    return array;
}

void disposeIntArray(int** array, int size) {
    if (array) {
        for (int i = 0; i < size; ++i)
            delete[] array[i];
        delete[] array;
    }
}

void disposeRealArray(RealOpenMM** array, int size) {
    if (array) {
        for (int i = 0; i < size; ++i)
            delete[] array[i];
        delete[] array;
    }
}

ReferenceCalcStandardMMForcesKernel::~ReferenceCalcStandardMMForcesKernel() {
    disposeIntArray(bondIndexArray, numBonds);
    disposeRealArray(bondParamArray, numBonds);
    disposeIntArray(angleIndexArray, numAngles);
    disposeRealArray(angleParamArray, numAngles);
    disposeIntArray(periodicTorsionIndexArray, numPeriodicTorsions);
    disposeRealArray(periodicTorsionParamArray, numPeriodicTorsions);
    disposeIntArray(rbTorsionIndexArray, numRBTorsions);
    disposeRealArray(rbTorsionParamArray, numRBTorsions);
}

void ReferenceCalcStandardMMForcesKernel::initialize(const vector<vector<int> >& bondIndices, const vector<vector<double> >& bondParameters,
        const vector<vector<int> >& angleIndices, const vector<vector<double> >& angleParameters,
        const vector<vector<int> >& periodicTorsionIndices, const vector<vector<double> >& periodicTorsionParameters,
        const vector<vector<int> >& rbTorsionIndices, const vector<vector<double> >& rbTorsionParameters,
        const vector<vector<int> >& bonded14Indices, const vector<set<int> >& exclusions,
        const vector<vector<double> >& nonbondedParameters) {
    numBonds = bondIndices.size();
    numAngles = angleIndices.size();
    numPeriodicTorsions = periodicTorsionIndices.size();
    numRBTorsions = rbTorsionIndices.size();
    bondIndexArray = copyToArray(bondIndices);
    bondParamArray = copyToArray(bondParameters);
    angleIndexArray = copyToArray(angleIndices);
    angleParamArray = copyToArray(angleParameters);
    periodicTorsionIndexArray = copyToArray(periodicTorsionIndices);
    periodicTorsionParamArray = copyToArray(periodicTorsionParameters);
    rbTorsionIndexArray = copyToArray(rbTorsionIndices);
    rbTorsionParamArray = copyToArray(rbTorsionParameters);
}

void ReferenceCalcStandardMMForcesKernel::execute(const Stream& positions, Stream& forces) {
    
}

void ReferenceCalcStandardMMEnergyKernel::initialize(const vector<vector<int> >& bondIndices, const vector<vector<double> >& bondParameters,
        const vector<vector<int> >& angleIndices, const vector<vector<double> >& angleParameters,
        const vector<vector<int> >& periodicTorsionIndices, const vector<vector<double> >& periodicTorsionParameters,
        const vector<vector<int> >& rbTorsionIndices, const vector<vector<double> >& rbTorsionParameters,
        const vector<vector<int> >& bonded14Indices, const vector<set<int> >& exclusions,
        const vector<vector<double> >& nonbondedParameters) {
    
}

double ReferenceCalcStandardMMEnergyKernel::execute(const Stream& positions) {
    return 0.0; // TODO implement correctly
}

void ReferenceCalcGBSAOBCForcesKernel::initialize(const vector<double>& bornRadii, const vector<vector<double> >& atomParameters,
        double solventDielectric, double soluteDielectric) {
    
}

void ReferenceCalcGBSAOBCForcesKernel::execute(const Stream& positions, Stream& forces) {
    
}

void ReferenceCalcGBSAOBCEnergyKernel::initialize(const vector<double>& bornRadii, const vector<vector<double> >& atomParameters,
        double solventDielectric, double soluteDielectric) {
    
}

double ReferenceCalcGBSAOBCEnergyKernel::execute(const Stream& positions) {
    return 0.0; // TODO implement correctly
}

void ReferenceIntegrateVerletStepKernel::initialize(const vector<double>& masses, const vector<vector<int> >& constraintIndices,
        const vector<double>& constraintLengths) {
    
}

void ReferenceIntegrateVerletStepKernel::execute(Stream& positions, Stream& velocities, const Stream& forces, double stepSize) {
    
}

void ReferenceIntegrateLangevinStepKernel::initialize(const vector<double>& masses, const vector<vector<int> >& constraintIndices,
        const vector<double>& constraintLengths) {
    
}

void ReferenceIntegrateLangevinStepKernel::execute(Stream& positions, Stream& velocities, const Stream& forces, double temperature, double friction, double stepSize) {
    
}

void ReferenceIntegrateBrownianStepKernel::initialize(const vector<double>& masses, const vector<vector<int> >& constraintIndices,
        const vector<double>& constraintLengths) {
    
}

void ReferenceIntegrateBrownianStepKernel::execute(Stream& positions, Stream& velocities, const Stream& forces, double temperature, double friction, double stepSize) {
    
}

void ReferenceApplyAndersenThermostatKernel::initialize(const vector<double>& masses) {
    
}

void ReferenceApplyAndersenThermostatKernel::execute(Stream& velocities, double temperature, double collisionFrequency, double stepSize) {
    
}

void ReferenceCalcKineticEnergyKernel::initialize(const vector<double>& masses) {
    
}

double ReferenceCalcKineticEnergyKernel::execute(const Stream& positions) {
    return 0.0; // TODO implement correctly
}

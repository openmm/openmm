/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                             *
 * Portions copyright (c) 2025 Stanford University and the Authors.            *
 * Authors: Muhammad Hasyim                                                    *
 * Contributors:                                                               *
 *                                                                             *
 * Permission is hereby granted, free of charge, to any person obtaining a     *
 * copy of this software and associated documentation files (the "Software"),  *
 * to deal in the Software without restriction, including without limitation   *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,    *
 * and/or sell copies of the Software, and to permit persons to whom the       *
 * Software is furnished to do so, subject to the following conditions:        *
 *                                                                             *
 * The above copyright notice and this permission notice shall be included in  *
 * all copies or substantial portions of the Software.                         *
 *                                                                             *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL     *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,     *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR       *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE   *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                      *
 * -------------------------------------------------------------------------- */

#include "openmm/UMAForce.h"
#include "openmm/internal/UMAForceImpl.h"

using namespace OpenMM;
using std::string;
using std::vector;
using std::map;

UMAForce::UMAForce(const string& modelName, TaskType taskType) :
        modelName(modelName), taskType(taskType), charge(0), spin(1), 
        deviceIndex(0), usePeriodic(false), cutoffRadius(0.6), maxNeighbors(50),
        inferenceSettings("default") {
}

void UMAForce::setModelPath(const string& path) {
    modelPath = path;
}

void UMAForce::setInferenceSettings(const string& settings) {
    inferenceSettings = settings;
}

void UMAForce::setDevice(int deviceIndex) {
    this->deviceIndex = deviceIndex;
}

void UMAForce::setCharge(int charge) {
    this->charge = charge;
}

void UMAForce::setSpin(int spin) {
    this->spin = spin;
}

void UMAForce::setAtomSubset(const vector<int>& atoms) {
    atomSubset = atoms;
}

void UMAForce::setAtomicReferences(const map<int, double>& refs) {
    atomicRefs = refs;
}

void UMAForce::setCutoffRadius(double cutoff) {
    cutoffRadius = cutoff;
}

void UMAForce::setMaxNeighbors(int maxNeighbors) {
    this->maxNeighbors = maxNeighbors;
}

bool UMAForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}

void UMAForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

ForceImpl* UMAForce::createImpl() const {
    return new UMAForceImpl(*this);
}

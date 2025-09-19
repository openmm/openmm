/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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

#include "openmm/OrientationRestraintForce.h"
#include "openmm/internal/OrientationRestraintForceImpl.h"

using namespace OpenMM;
using namespace std;

OrientationRestraintForce::OrientationRestraintForce(double k, const vector<Vec3>& referencePositions, const vector<int>& particles) :
        k(k), referencePositions(referencePositions), particles(particles) {
}

void OrientationRestraintForce::setK(double k) {
    this->k = k;
}

void OrientationRestraintForce::setReferencePositions(const std::vector<Vec3>& positions) {
    referencePositions = positions;
}

void OrientationRestraintForce::setParticles(const std::vector<int>& particles) {
    this->particles = particles;
}

void OrientationRestraintForce::updateParametersInContext(Context& context) {
    dynamic_cast<OrientationRestraintForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

ForceImpl* OrientationRestraintForce::createImpl() const {
    return new OrientationRestraintForceImpl(*this);
}

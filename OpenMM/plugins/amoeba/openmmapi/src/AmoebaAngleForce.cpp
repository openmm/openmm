/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
 * Authors:                                                                   *
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

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/AmoebaAngleForce.h"
#include "openmm/internal/AmoebaAngleForceImpl.h"

using namespace OpenMM;

AmoebaAngleForce::AmoebaAngleForce() : usePeriodic(false) {
    _globalCubicK = _globalQuarticK = _globalPenticK = _globalSexticK = 0.0;
}

int AmoebaAngleForce::addAngle(int particle1, int particle2, int particle3,  double length, double quadraticK) {
    angles.push_back(AngleInfo(particle1, particle2, particle3, length, quadraticK));
    return angles.size()-1;
}

void AmoebaAngleForce::getAngleParameters(int index, int& particle1, int& particle2, int& particle3,
                                                  double& length, double&  quadraticK) const {
    particle1       = angles[index].particle1;
    particle2       = angles[index].particle2;
    particle3       = angles[index].particle3;
    length          = angles[index].length;
    quadraticK      = angles[index].quadraticK;
}

void AmoebaAngleForce::setAngleParameters(int index, int particle1, int particle2, int particle3, 
                                                  double length, double quadraticK) {
    angles[index].particle1  = particle1;
    angles[index].particle2  = particle2;
    angles[index].particle3  = particle3;
    angles[index].length     = length;
    angles[index].quadraticK = quadraticK;
}

double AmoebaAngleForce::getAmoebaGlobalAngleCubic() const {
    return _globalCubicK;
}

void AmoebaAngleForce::setAmoebaGlobalAngleCubic(double cubicK) {
    _globalCubicK           = cubicK;
}

double AmoebaAngleForce::getAmoebaGlobalAngleQuartic() const {
    return _globalQuarticK;
}

void AmoebaAngleForce::setAmoebaGlobalAngleQuartic(double quarticK) {
    _globalQuarticK         = quarticK;
}

double AmoebaAngleForce::getAmoebaGlobalAnglePentic() const {
    return _globalPenticK;
}

void AmoebaAngleForce::setAmoebaGlobalAnglePentic(double penticK) {
    _globalPenticK           = penticK;
}

double AmoebaAngleForce::getAmoebaGlobalAngleSextic() const {
    return _globalSexticK;
}

void AmoebaAngleForce::setAmoebaGlobalAngleSextic(double sexticK) {
    _globalSexticK         = sexticK;
}

ForceImpl* AmoebaAngleForce::createImpl() const {
    return new AmoebaAngleForceImpl(*this);
}

void AmoebaAngleForce::updateParametersInContext(Context& context) {
    dynamic_cast<AmoebaAngleForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

void AmoebaAngleForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool AmoebaAngleForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}

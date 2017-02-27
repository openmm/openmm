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
#include "openmm/AmoebaInPlaneAngleForce.h"
#include "openmm/internal/AmoebaInPlaneAngleForceImpl.h"

using namespace OpenMM;

AmoebaInPlaneAngleForce::AmoebaInPlaneAngleForce() : usePeriodic(false) {
    _globalCubicK = _globalQuarticK = _globalPenticK = _globalSexticK = 0.0;
}

int AmoebaInPlaneAngleForce::addAngle(int particle1, int particle2, int particle3, int particle4,  double length, double quadraticK) {
    angles.push_back(AngleInfo(particle1, particle2, particle3, particle4, length, quadraticK));
    return angles.size()-1;
}

void AmoebaInPlaneAngleForce::getAngleParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4,
                                                  double& length, double&  quadraticK) const {
    particle1       = angles[index].particle1;
    particle2       = angles[index].particle2;
    particle3       = angles[index].particle3;
    particle4       = angles[index].particle4;
    length          = angles[index].length;
    quadraticK      = angles[index].quadraticK;
}

void AmoebaInPlaneAngleForce::setAngleParameters(int index, int particle1, int particle2, int particle3, int particle4,
                                                    double length, double quadraticK) {
    angles[index].particle1  = particle1;
    angles[index].particle2  = particle2;
    angles[index].particle3  = particle3;
    angles[index].particle4  = particle4;
    angles[index].length     = length;
    angles[index].quadraticK = quadraticK;
}

void AmoebaInPlaneAngleForce::setAmoebaGlobalInPlaneAngleCubic(double cubicK) {
    _globalCubicK           = cubicK;
}

void AmoebaInPlaneAngleForce::setAmoebaGlobalInPlaneAngleQuartic(double quarticK) {
    _globalQuarticK         = quarticK;
}

double AmoebaInPlaneAngleForce::getAmoebaGlobalInPlaneAngleCubic() const {
    return _globalCubicK;
}

double AmoebaInPlaneAngleForce::getAmoebaGlobalInPlaneAngleQuartic() const {
    return _globalQuarticK;
}

void AmoebaInPlaneAngleForce::setAmoebaGlobalInPlaneAnglePentic(double cubicK) {
    _globalPenticK           = cubicK;
}

void AmoebaInPlaneAngleForce::setAmoebaGlobalInPlaneAngleSextic(double quarticK) {
    _globalSexticK         = quarticK;
}

double AmoebaInPlaneAngleForce::getAmoebaGlobalInPlaneAnglePentic() const {
    return _globalPenticK;
}

double AmoebaInPlaneAngleForce::getAmoebaGlobalInPlaneAngleSextic() const {
    return _globalSexticK;
}

ForceImpl* AmoebaInPlaneAngleForce::createImpl() const {
    return new AmoebaInPlaneAngleForceImpl(*this);
}

void AmoebaInPlaneAngleForce::updateParametersInContext(Context& context) {
    dynamic_cast<AmoebaInPlaneAngleForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

void AmoebaInPlaneAngleForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool AmoebaInPlaneAngleForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}

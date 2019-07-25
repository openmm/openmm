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
#include "openmm/AmoebaBondForce.h"
#include "openmm/internal/AmoebaBondForceImpl.h"

using namespace OpenMM;

AmoebaBondForce::AmoebaBondForce() : usePeriodic(false) {
   _globalCubicK = _globalQuarticK = 0.0;
}

int AmoebaBondForce::addBond(int particle1, int particle2, double length, double quadraticK) {
    bonds.push_back(BondInfo(particle1, particle2, length, quadraticK));
    return bonds.size()-1;
}

void AmoebaBondForce::getBondParameters(int index, int& particle1, int& particle2, double& length, double&  quadraticK) const {
    particle1       = bonds[index].particle1;
    particle2       = bonds[index].particle2;
    length          = bonds[index].length;
    quadraticK      = bonds[index].quadraticK;
}

void AmoebaBondForce::setBondParameters(int index, int particle1, int particle2, double length, double quadraticK) {
    bonds[index].particle1  = particle1;
    bonds[index].particle2  = particle2;
    bonds[index].length     = length;
    bonds[index].quadraticK = quadraticK;
}

void AmoebaBondForce::setAmoebaGlobalBondCubic(double cubicK) {
    _globalCubicK           = cubicK;
}

void AmoebaBondForce::setAmoebaGlobalBondQuartic(double quarticK) {
    _globalQuarticK         = quarticK;
}

double AmoebaBondForce::getAmoebaGlobalBondCubic() const {
    return _globalCubicK;
}

double AmoebaBondForce::getAmoebaGlobalBondQuartic() const {
    return _globalQuarticK;
}

ForceImpl* AmoebaBondForce::createImpl() const {
    return new AmoebaBondForceImpl(*this);
}

void AmoebaBondForce::updateParametersInContext(Context& context) {
    dynamic_cast<AmoebaBondForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

void AmoebaBondForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool AmoebaBondForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}

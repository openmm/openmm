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
#include "openmm/AmoebaPiTorsionForce.h"
#include "openmm/internal/AmoebaPiTorsionForceImpl.h"

using namespace OpenMM;

AmoebaPiTorsionForce::AmoebaPiTorsionForce() : usePeriodic(false) {
}

int AmoebaPiTorsionForce::addPiTorsion(int particle1, int particle2, int particle3, int particle4, int particle5, int particle6,  double k) {
    piTorsions.push_back(PiTorsionInfo(particle1, particle2, particle3, particle4, particle5, particle6, k));
    return piTorsions.size()-1;
}

void AmoebaPiTorsionForce::getPiTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, int& particle5, int& particle6, double& k) const {
    particle1       = piTorsions[index].particle1;
    particle2       = piTorsions[index].particle2;
    particle3       = piTorsions[index].particle3;
    particle4       = piTorsions[index].particle4;
    particle5       = piTorsions[index].particle5;
    particle6       = piTorsions[index].particle6;
    k               = piTorsions[index].k;
}

void AmoebaPiTorsionForce::setPiTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, int particle5, int particle6,  double k) {
    piTorsions[index].particle1  = particle1;
    piTorsions[index].particle2  = particle2;
    piTorsions[index].particle3  = particle3;
    piTorsions[index].particle4  = particle4;
    piTorsions[index].particle5  = particle5;
    piTorsions[index].particle6  = particle6;
    piTorsions[index].k          = k;
}

ForceImpl* AmoebaPiTorsionForce::createImpl() const {
    return new AmoebaPiTorsionForceImpl(*this);
}

void AmoebaPiTorsionForce::updateParametersInContext(Context& context) {
    dynamic_cast<AmoebaPiTorsionForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

void AmoebaPiTorsionForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool AmoebaPiTorsionForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}

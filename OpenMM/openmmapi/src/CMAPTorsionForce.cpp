/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2016 Stanford University and the Authors.      *
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

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/CMAPTorsionForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/CMAPTorsionForceImpl.h"

using namespace OpenMM;

CMAPTorsionForce::CMAPTorsionForce() : usePeriodic(false) {
}

int CMAPTorsionForce::addMap(int size, const std::vector<double>& energy) {
    if (energy.size() != size*size)
        throw OpenMMException("CMAPTorsionForce: incorrect number of energy values");
    maps.push_back(MapInfo(size, energy));
    return maps.size()-1;
}

void CMAPTorsionForce::getMapParameters(int index, int& size, std::vector<double>& energy) const {
    ASSERT_VALID_INDEX(index, maps);
    size = maps[index].size;
    energy = maps[index].energy;
}

void CMAPTorsionForce::setMapParameters(int index, int size, const std::vector<double>& energy) {
    ASSERT_VALID_INDEX(index, maps);
    if (energy.size() != size*size)
        throw OpenMMException("CMAPTorsionForce: incorrect number of energy values");
    maps[index].size = size;
    maps[index].energy = energy;
}

int CMAPTorsionForce::addTorsion(int map, int a1, int a2, int a3, int a4, int b1, int b2, int b3, int b4) {
    torsions.push_back(CMAPTorsionInfo(map, a1, a2, a3, a4, b1, b2, b3, b4));
    return torsions.size()-1;
}

void CMAPTorsionForce::getTorsionParameters(int index, int& map, int& a1, int& a2, int& a3, int& a4, int& b1, int& b2, int& b3, int& b4) const {
    ASSERT_VALID_INDEX(index, torsions);
    map = torsions[index].map;
    a1 = torsions[index].a1;
    a2 = torsions[index].a2;
    a3 = torsions[index].a3;
    a4 = torsions[index].a4;
    b1 = torsions[index].b1;
    b2 = torsions[index].b2;
    b3 = torsions[index].b3;
    b4 = torsions[index].b4;
}

void CMAPTorsionForce::setTorsionParameters(int index, int map, int a1, int a2, int a3, int a4, int b1, int b2, int b3, int b4) {
    ASSERT_VALID_INDEX(index, torsions);
    torsions[index].map = map;
    torsions[index].a1 = a1;
    torsions[index].a2 = a2;
    torsions[index].a3 = a3;
    torsions[index].a4 = a4;
    torsions[index].b1 = b1;
    torsions[index].b2 = b2;
    torsions[index].b3 = b3;
    torsions[index].b4 = b4;
}

ForceImpl* CMAPTorsionForce::createImpl() const {
    return new CMAPTorsionForceImpl(*this);
}

void CMAPTorsionForce::updateParametersInContext(Context& context) {
    dynamic_cast<CMAPTorsionForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

void CMAPTorsionForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool CMAPTorsionForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}

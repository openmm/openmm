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

#include "Force.h"
#include "OpenMMException.h"
#include "PeriodicTorsionForce.h"
#include "internal/PeriodicTorsionForceImpl.h"

using namespace OpenMM;

PeriodicTorsionForce::PeriodicTorsionForce(int numTorsions) : periodicTorsions(numTorsions) {
}

void PeriodicTorsionForce::getTorsionParameters(int index, int& atom1, int& atom2, int& atom3, int& atom4, int& periodicity, double& phase, double& k) const {
    atom1 = periodicTorsions[index].atom1;
    atom2 = periodicTorsions[index].atom2;
    atom3 = periodicTorsions[index].atom3;
    atom4 = periodicTorsions[index].atom4;
    periodicity = periodicTorsions[index].periodicity;
    phase = periodicTorsions[index].phase;
    k = periodicTorsions[index].k;
}

void PeriodicTorsionForce::setTorsionParameters(int index, int atom1, int atom2, int atom3, int atom4, int periodicity, double phase, double k) {
    periodicTorsions[index].atom1 = atom1;
    periodicTorsions[index].atom2 = atom2;
    periodicTorsions[index].atom3 = atom3;
    periodicTorsions[index].atom4 = atom4;
    periodicTorsions[index].periodicity = periodicity;
    periodicTorsions[index].phase = phase;
    periodicTorsions[index].k = k;
}

ForceImpl* PeriodicTorsionForce::createImpl() {
    return new PeriodicTorsionForceImpl(*this);
}

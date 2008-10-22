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
#include "HarmonicAngleForce.h"
#include "internal/HarmonicAngleForceImpl.h"

using namespace OpenMM;

HarmonicAngleForce::HarmonicAngleForce(int numAngles) : angles(numAngles) {
}

void HarmonicAngleForce::getAngleParameters(int index, int& atom1, int& atom2, int& atom3, double& angle, double& k) const {
    atom1 = angles[index].atom1;
    atom2 = angles[index].atom2;
    atom3 = angles[index].atom3;
    angle = angles[index].angle;
    k = angles[index].k;
}

void HarmonicAngleForce::setAngleParameters(int index, int atom1, int atom2, int atom3, double angle, double k) {
    angles[index].atom1 = atom1;
    angles[index].atom2 = atom2;
    angles[index].atom3 = atom3;
    angles[index].angle = angle;
    angles[index].k = k;
}

ForceImpl* HarmonicAngleForce::createImpl() {
    return new HarmonicAngleForceImpl(*this);
}

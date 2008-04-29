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
#include "GBSAOBCForceField.h"
#include "internal/GBSAOBCForceFieldImpl.h"

using namespace OpenMM;

GBSAOBCForceField::GBSAOBCForceField(int numAtoms) : atoms(numAtoms), solventDielectric(78.3), soluteDielectric(1.0) {
}

void GBSAOBCForceField::getAtomParameters(int index, double& charge, double& radius, double& scalingFactor) const {
    charge = atoms[index].charge;
    radius = atoms[index].radius;
    scalingFactor = atoms[index].scalingFactor;
}

void GBSAOBCForceField::setAtomParameters(int index, double charge, double radius, double scalingFactor) {
    atoms[index].charge = charge;
    atoms[index].radius = radius;
    atoms[index].scalingFactor = scalingFactor;
}

ForceImpl* GBSAOBCForceField::createImpl(OpenMMContextImpl& context) {
    return new GBSAOBCForceFieldImpl(*this, context);
}

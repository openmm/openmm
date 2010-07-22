/* -------------------------------------------------------------------------- *
 *                                AmoebaOpenMM                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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
#include "AmoebaHarmonicBondForce.h"
#include "internal/AmoebaHarmonicBondForceImpl.h"

using namespace OpenMM;

const double AmoebaHarmonicBondForce::DEFAULT_GLOBAL_K = -1234.567891234;

AmoebaHarmonicBondForce::AmoebaHarmonicBondForce() {
   _globalCubicK = _globalQuarticK = 0.0;
}

int AmoebaHarmonicBondForce::addBond(int particle1, int particle2, double length, double quadraticK, double cubicK, double quarticK) {
    bonds.push_back(BondInfo(particle1, particle2, length, quadraticK, cubicK, quarticK));
    return bonds.size()-1;
}

void AmoebaHarmonicBondForce::getBondParameters(int index, int& particle1, int& particle2, double& length, double&  quadraticK, double& cubicK, double& quarticK ) const {
    particle1       = bonds[index].particle1;
    particle2       = bonds[index].particle2;
    length          = bonds[index].length;
    quadraticK      = bonds[index].quadraticK;
    cubicK          = bonds[index].cubicK   == DEFAULT_GLOBAL_K ? _globalCubicK   : bonds[index].cubicK;
    quarticK        = bonds[index].quarticK == DEFAULT_GLOBAL_K ? _globalQuarticK : bonds[index].quarticK;
}

void AmoebaHarmonicBondForce::setBondParameters(int index, int particle1, int particle2, double length, double quadraticK, double cubicK, double quarticK ) {
    bonds[index].particle1  = particle1;
    bonds[index].particle2  = particle2;
    bonds[index].length     = length;
    bonds[index].quadraticK = quadraticK;
    bonds[index].cubicK     = cubicK;
    bonds[index].quarticK   = quarticK;
}

void AmoebaHarmonicBondForce::setAmoebaGlobalHarmonicBondCubic(double cubicK ) {
    _globalCubicK           = cubicK;
}

void AmoebaHarmonicBondForce::setAmoebaGlobalHarmonicBondQuartic(double quarticK ) {
    _globalQuarticK         = quarticK;
}

double AmoebaHarmonicBondForce::getAmoebaGlobalHarmonicBondCubic( void ) const {
    return _globalCubicK;
}

double AmoebaHarmonicBondForce::getAmoebaGlobalHarmonicBondQuartic( void ) const {
    return _globalQuarticK;
}

ForceImpl* AmoebaHarmonicBondForce::createImpl() {
    return new AmoebaHarmonicBondForceImpl(*this);
}

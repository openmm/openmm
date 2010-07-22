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
#include "AmoebaHarmonicInPlaneAngleForce.h"
#include "internal/AmoebaHarmonicInPlaneAngleForceImpl.h"

using namespace OpenMM;

AmoebaHarmonicInPlaneAngleForce::AmoebaHarmonicInPlaneAngleForce() {
    _globalCubicK = _globalQuarticK = _globalPenticK = _globalSexticK = 0.0;
}

int AmoebaHarmonicInPlaneAngleForce::addAngle(int particle1, int particle2, int particle3, int particle4,  double length, double quadraticK ) {
    angles.push_back(AngleInfo(particle1, particle2, particle3, particle4, length, quadraticK ));
    return angles.size()-1;
}

void AmoebaHarmonicInPlaneAngleForce::getAngleParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4,
                                                  double& length, double&  quadraticK ) const {
    particle1       = angles[index].particle1;
    particle2       = angles[index].particle2;
    particle3       = angles[index].particle3;
    particle4       = angles[index].particle4;
    length          = angles[index].length;
    quadraticK      = angles[index].quadraticK;
}

void AmoebaHarmonicInPlaneAngleForce::setAngleParameters(int index, int particle1, int particle2, int particle3, int particle4,
                                                    double length, double quadraticK ) {
    angles[index].particle1  = particle1;
    angles[index].particle2  = particle2;
    angles[index].particle3  = particle3;
    angles[index].particle4  = particle4;
    angles[index].length     = length;
    angles[index].quadraticK = quadraticK;
}

void AmoebaHarmonicInPlaneAngleForce::setAmoebaGlobalHarmonicInPlaneAngleCubic(double cubicK ) {
    _globalCubicK           = cubicK;
}

void AmoebaHarmonicInPlaneAngleForce::setAmoebaGlobalHarmonicInPlaneAngleQuartic(double quarticK ) {
    _globalQuarticK         = quarticK;
}

double AmoebaHarmonicInPlaneAngleForce::getAmoebaGlobalHarmonicInPlaneAngleCubic( void ) const {
    return _globalCubicK;
}

double AmoebaHarmonicInPlaneAngleForce::getAmoebaGlobalHarmonicInPlaneAngleQuartic( void ) const {
    return _globalQuarticK;
}

void AmoebaHarmonicInPlaneAngleForce::setAmoebaGlobalHarmonicInPlaneAnglePentic(double cubicK ) {
    _globalPenticK           = cubicK;
}

void AmoebaHarmonicInPlaneAngleForce::setAmoebaGlobalHarmonicInPlaneAngleSextic(double quarticK ) {
    _globalSexticK         = quarticK;
}

double AmoebaHarmonicInPlaneAngleForce::getAmoebaGlobalHarmonicInPlaneAnglePentic( void ) const {
    return _globalPenticK;
}

double AmoebaHarmonicInPlaneAngleForce::getAmoebaGlobalHarmonicInPlaneAngleSextic( void ) const {
    return _globalSexticK;
}

ForceImpl* AmoebaHarmonicInPlaneAngleForce::createImpl() {
    return new AmoebaHarmonicInPlaneAngleForceImpl(*this);
}

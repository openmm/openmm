/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2021 Stanford University and the Authors.           *
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

#include "ReferencePointFunctions.h"
#include "ReferenceBondIxn.h"
#include "openmm/OpenMMException.h"
#include <cmath>

using namespace OpenMM;

ReferencePointDistanceFunction::ReferencePointDistanceFunction(bool periodic, Vec3** boxVectorHandle) : periodic(periodic), boxVectorHandle(boxVectorHandle) {
}

int ReferencePointDistanceFunction::getNumArguments() const {
    return 6;
}

double ReferencePointDistanceFunction::evaluate(const double* arguments) const {
    Vec3 delta = Vec3(arguments[0], arguments[1], arguments[2])-Vec3(arguments[3], arguments[4], arguments[5]);
    if (periodic) {
        Vec3* boxVectors = *boxVectorHandle;
        delta -= boxVectors[2]*floor(delta[2]/boxVectors[2][2]+0.5);
        delta -= boxVectors[1]*floor(delta[1]/boxVectors[1][1]+0.5);
        delta -= boxVectors[0]*floor(delta[0]/boxVectors[0][0]+0.5);
    }
    return sqrt(delta.dot(delta));
}

double ReferencePointDistanceFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    int argIndex = -1;
    for (int i = 0; i < 6; i++) {
        if (derivOrder[i] > 0) {
            if (derivOrder[i] > 1 || argIndex != -1)
                throw OpenMMException("Unsupported derivative of pointdistance"); // Should be impossible for this to happen.
            argIndex = i;
        }
    }
    Vec3 delta = Vec3(arguments[0], arguments[1], arguments[2])-Vec3(arguments[3], arguments[4], arguments[5]);
    if (periodic) {
        Vec3* boxVectors = *boxVectorHandle;
        delta -= boxVectors[2]*floor(delta[2]/boxVectors[2][2]+0.5);
        delta -= boxVectors[1]*floor(delta[1]/boxVectors[1][1]+0.5);
        delta -= boxVectors[0]*floor(delta[0]/boxVectors[0][0]+0.5);
    }
    double r = sqrt(delta.dot(delta));
    if (r == 0)
        return 0.0;    
    if (argIndex < 3)
        return delta[argIndex]/r;
    return -delta[argIndex-3]/r;
}

Lepton::CustomFunction* ReferencePointDistanceFunction::clone() const {
    return new ReferencePointDistanceFunction(periodic, boxVectorHandle);
}

ReferencePointAngleFunction::ReferencePointAngleFunction(bool periodic, Vec3** boxVectorHandle) : periodic(periodic), boxVectorHandle(boxVectorHandle) {
}

int ReferencePointAngleFunction::getNumArguments() const {
    return 9;
}

double ReferencePointAngleFunction::evaluate(const double* arguments) const {
    Vec3 delta12 = Vec3(arguments[3], arguments[4], arguments[5])-Vec3(arguments[0], arguments[1], arguments[2]);
    Vec3 delta32 = Vec3(arguments[3], arguments[4], arguments[5])-Vec3(arguments[6], arguments[7], arguments[8]);
    if (periodic) {
        Vec3* boxVectors = *boxVectorHandle;
        delta12 -= boxVectors[2]*floor(delta12[2]/boxVectors[2][2]+0.5);
        delta12 -= boxVectors[1]*floor(delta12[1]/boxVectors[1][1]+0.5);
        delta12 -= boxVectors[0]*floor(delta12[0]/boxVectors[0][0]+0.5);
        delta32 -= boxVectors[2]*floor(delta32[2]/boxVectors[2][2]+0.5);
        delta32 -= boxVectors[1]*floor(delta32[1]/boxVectors[1][1]+0.5);
        delta32 -= boxVectors[0]*floor(delta32[0]/boxVectors[0][0]+0.5);
    }
    return ReferenceBondIxn::getAngleBetweenTwoVectors(&delta12[0], &delta32[0]);
}

double ReferencePointAngleFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    int argIndex = -1;
    for (int i = 0; i < 9; i++) {
        if (derivOrder[i] > 0) {
            if (derivOrder[i] > 1 || argIndex != -1)
                throw OpenMMException("Unsupported derivative of pointangle"); // Should be impossible for this to happen.
            argIndex = i;
        }
    }
    Vec3 delta12 = Vec3(arguments[3], arguments[4], arguments[5])-Vec3(arguments[0], arguments[1], arguments[2]);
    Vec3 delta32 = Vec3(arguments[3], arguments[4], arguments[5])-Vec3(arguments[6], arguments[7], arguments[8]);
    if (periodic) {
        Vec3* boxVectors = *boxVectorHandle;
        delta12 -= boxVectors[2]*floor(delta12[2]/boxVectors[2][2]+0.5);
        delta12 -= boxVectors[1]*floor(delta12[1]/boxVectors[1][1]+0.5);
        delta12 -= boxVectors[0]*floor(delta12[0]/boxVectors[0][0]+0.5);
        delta32 -= boxVectors[2]*floor(delta32[2]/boxVectors[2][2]+0.5);
        delta32 -= boxVectors[1]*floor(delta32[1]/boxVectors[1][1]+0.5);
        delta32 -= boxVectors[0]*floor(delta32[0]/boxVectors[0][0]+0.5);
    }
    Vec3 thetaCross = delta12.cross(delta32);
    double lengthThetaCross = sqrt(thetaCross.dot(thetaCross));
    if (lengthThetaCross < 1.0e-6)
        lengthThetaCross = 1.0e-6;
    Vec3 deltaCrossP[3];
    deltaCrossP[0] = delta12.cross(thetaCross)/(delta12.dot(delta12)*lengthThetaCross);
    deltaCrossP[2] = -delta32.cross(thetaCross)/(delta32.dot(delta32)*lengthThetaCross);
    deltaCrossP[1] = -(deltaCrossP[0]+deltaCrossP[2]);
    return -deltaCrossP[argIndex/3][argIndex%3];
}

Lepton::CustomFunction* ReferencePointAngleFunction::clone() const {
    return new ReferencePointAngleFunction(periodic, boxVectorHandle);
}

ReferencePointDihedralFunction::ReferencePointDihedralFunction(bool periodic, Vec3** boxVectorHandle) : periodic(periodic), boxVectorHandle(boxVectorHandle) {
}

int ReferencePointDihedralFunction::getNumArguments() const {
    return 12;
}

double ReferencePointDihedralFunction::evaluate(const double* arguments) const {
    Vec3 delta12 = Vec3(arguments[0], arguments[1], arguments[2])-Vec3(arguments[3], arguments[4], arguments[5]);
    Vec3 delta32 = Vec3(arguments[6], arguments[7], arguments[8])-Vec3(arguments[3], arguments[4], arguments[5]);
    Vec3 delta34 = Vec3(arguments[6], arguments[7], arguments[8])-Vec3(arguments[9], arguments[10], arguments[11]);
    if (periodic) {
        Vec3* boxVectors = *boxVectorHandle;
        delta12 -= boxVectors[2]*floor(delta12[2]/boxVectors[2][2]+0.5);
        delta12 -= boxVectors[1]*floor(delta12[1]/boxVectors[1][1]+0.5);
        delta12 -= boxVectors[0]*floor(delta12[0]/boxVectors[0][0]+0.5);
        delta32 -= boxVectors[2]*floor(delta32[2]/boxVectors[2][2]+0.5);
        delta32 -= boxVectors[1]*floor(delta32[1]/boxVectors[1][1]+0.5);
        delta32 -= boxVectors[0]*floor(delta32[0]/boxVectors[0][0]+0.5);
        delta34 -= boxVectors[2]*floor(delta34[2]/boxVectors[2][2]+0.5);
        delta34 -= boxVectors[1]*floor(delta34[1]/boxVectors[1][1]+0.5);
        delta34 -= boxVectors[0]*floor(delta34[0]/boxVectors[0][0]+0.5);
    }
    return ReferenceBondIxn::getDihedralAngleBetweenThreeVectors(&delta12[0], &delta32[0], &delta34[0], NULL, NULL, &delta12[0]);
}

double ReferencePointDihedralFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    int argIndex = -1;
    for (int i = 0; i < 12; i++) {
        if (derivOrder[i] > 0) {
            if (derivOrder[i] > 1 || argIndex != -1)
                throw OpenMMException("Unsupported derivative of pointdihedral"); // Should be impossible for this to happen.
            argIndex = i;
        }
    }
    Vec3 delta12 = Vec3(arguments[0], arguments[1], arguments[2])-Vec3(arguments[3], arguments[4], arguments[5]);
    Vec3 delta32 = Vec3(arguments[6], arguments[7], arguments[8])-Vec3(arguments[3], arguments[4], arguments[5]);
    Vec3 delta34 = Vec3(arguments[6], arguments[7], arguments[8])-Vec3(arguments[9], arguments[10], arguments[11]);
    if (periodic) {
        Vec3* boxVectors = *boxVectorHandle;
        delta12 -= boxVectors[2]*floor(delta12[2]/boxVectors[2][2]+0.5);
        delta12 -= boxVectors[1]*floor(delta12[1]/boxVectors[1][1]+0.5);
        delta12 -= boxVectors[0]*floor(delta12[0]/boxVectors[0][0]+0.5);
        delta32 -= boxVectors[2]*floor(delta32[2]/boxVectors[2][2]+0.5);
        delta32 -= boxVectors[1]*floor(delta32[1]/boxVectors[1][1]+0.5);
        delta32 -= boxVectors[0]*floor(delta32[0]/boxVectors[0][0]+0.5);
        delta34 -= boxVectors[2]*floor(delta34[2]/boxVectors[2][2]+0.5);
        delta34 -= boxVectors[1]*floor(delta34[1]/boxVectors[1][1]+0.5);
        delta34 -= boxVectors[0]*floor(delta34[0]/boxVectors[0][0]+0.5);
    }
    Vec3 cross1 = delta12.cross(delta32);
    Vec3 cross2 = delta32.cross(delta34);
    double norm2Cross1 = cross1.dot(cross1);
    double norm2Cross2 = cross2.dot(cross2);
    double norm2Delta32 = delta32.dot(delta32);
    double normDelta32 = sqrt(norm2Delta32);
    double forceFactors[4];
    forceFactors[0] = -normDelta32/norm2Cross1;
    forceFactors[3] = normDelta32/norm2Cross2;
    forceFactors[1] = delta12.dot(delta32);
    forceFactors[1] /= norm2Delta32;
    forceFactors[2] = delta34.dot(delta32);
    forceFactors[2] /= norm2Delta32;
    Vec3 internalF[4];
    internalF[0] = forceFactors[0]*cross1;
    internalF[3] = forceFactors[3]*cross2;
    Vec3 s = forceFactors[1]*internalF[0] - forceFactors[2]*internalF[3];
    internalF[1] = internalF[0] - s;
    internalF[2] = internalF[3] + s;
    internalF[0] *= -1;
    internalF[3] *= -1;
    return internalF[argIndex/3][argIndex%3];
}

Lepton::CustomFunction* ReferencePointDihedralFunction::clone() const {
    return new ReferencePointDihedralFunction(periodic, boxVectorHandle);
}

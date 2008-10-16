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
#include "StandardMMForceField.h"
#include "internal/StandardMMForceFieldImpl.h"

using namespace OpenMM;

StandardMMForceField::StandardMMForceField(int numAtoms, int numBonds, int numAngles, int numPeriodicTorsions, int numRBTorsions, int numNonbonded14) :
        atoms(numAtoms), bonds(numBonds), angles(numAngles), periodicTorsions(numPeriodicTorsions), rbTorsions(numRBTorsions), nb14s(numNonbonded14),
        nonbondedMethod(NoCutoff), cutoffDistance(1.0) {
    periodicBoxSize[0] = periodicBoxSize[1] = periodicBoxSize[2] = 2.0;
}

StandardMMForceField::NonbondedMethod StandardMMForceField::getNonbondedMethod() const {
    return nonbondedMethod;
}

void StandardMMForceField::setNonbondedMethod(NonbondedMethod method) {
    nonbondedMethod = method;
}

double StandardMMForceField::getCutoffDistance() const {
    return cutoffDistance;
}

void StandardMMForceField::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

void StandardMMForceField::getPeriodicBoxSize(double& x, double& y, double& z) const {
    x = periodicBoxSize[0];
    y = periodicBoxSize[1];
    z = periodicBoxSize[2];
}

void StandardMMForceField::setPeriodicBoxSize(double x, double y, double z) {
    periodicBoxSize[0] = x;
    periodicBoxSize[1] = y;
    periodicBoxSize[2] = z;
}

void StandardMMForceField::getAtomParameters(int index, double& charge, double& radius, double& depth) const {
    charge = atoms[index].charge;
    radius = atoms[index].radius;
    depth = atoms[index].depth;
}

void StandardMMForceField::setAtomParameters(int index, double charge, double radius, double depth) {
    atoms[index].charge = charge;
    atoms[index].radius = radius;
    atoms[index].depth = depth;
}

void StandardMMForceField::getBondParameters(int index, int& atom1, int& atom2, double& length, double& k) const {
    atom1 = bonds[index].atom1;
    atom2 = bonds[index].atom2;
    length = bonds[index].length;
    k = bonds[index].k;
}

void StandardMMForceField::setBondParameters(int index, int atom1, int atom2, double length, double k) {
    bonds[index].atom1 = atom1;
    bonds[index].atom2 = atom2;
    bonds[index].length = length;
    bonds[index].k = k;
}

void StandardMMForceField::getAngleParameters(int index, int& atom1, int& atom2, int& atom3, double& angle, double& k) const {
    atom1 = angles[index].atom1;
    atom2 = angles[index].atom2;
    atom3 = angles[index].atom3;
    angle = angles[index].angle;
    k = angles[index].k;
}

void StandardMMForceField::setAngleParameters(int index, int atom1, int atom2, int atom3, double angle, double k) {
    angles[index].atom1 = atom1;
    angles[index].atom2 = atom2;
    angles[index].atom3 = atom3;
    angles[index].angle = angle;
    angles[index].k = k;
}

void StandardMMForceField::getPeriodicTorsionParameters(int index, int& atom1, int& atom2, int& atom3, int& atom4, int& periodicity, double& phase, double& k) const {
    atom1 = periodicTorsions[index].atom1;
    atom2 = periodicTorsions[index].atom2;
    atom3 = periodicTorsions[index].atom3;
    atom4 = periodicTorsions[index].atom4;
    periodicity = periodicTorsions[index].periodicity;
    phase = periodicTorsions[index].phase;
    k = periodicTorsions[index].k;
}

void StandardMMForceField::setPeriodicTorsionParameters(int index, int atom1, int atom2, int atom3, int atom4, int periodicity, double phase, double k) {
    periodicTorsions[index].atom1 = atom1;
    periodicTorsions[index].atom2 = atom2;
    periodicTorsions[index].atom3 = atom3;
    periodicTorsions[index].atom4 = atom4;
    periodicTorsions[index].periodicity = periodicity;
    periodicTorsions[index].phase = phase;
    periodicTorsions[index].k = k;
}

void StandardMMForceField::getRBTorsionParameters(int index, int& atom1, int& atom2, int& atom3, int& atom4, double& c0, double& c1, double& c2, double& c3, double& c4, double& c5) const {
    atom1 = rbTorsions[index].atom1;
    atom2 = rbTorsions[index].atom2;
    atom3 = rbTorsions[index].atom3;
    atom4 = rbTorsions[index].atom4;
    c0 = rbTorsions[index].c[0];
    c1 = rbTorsions[index].c[1];
    c2 = rbTorsions[index].c[2];
    c3 = rbTorsions[index].c[3];
    c4 = rbTorsions[index].c[4];
    c5 = rbTorsions[index].c[5];
}

void StandardMMForceField::setRBTorsionParameters(int index, int atom1, int atom2, int atom3, int atom4, double c0, double c1, double c2, double c3, double c4, double c5) {
    rbTorsions[index].atom1 = atom1;
    rbTorsions[index].atom2 = atom2;
    rbTorsions[index].atom3 = atom3;
    rbTorsions[index].atom4 = atom4;
    rbTorsions[index].c[0] = c0;
    rbTorsions[index].c[1] = c1;
    rbTorsions[index].c[2] = c2;
    rbTorsions[index].c[3] = c3;
    rbTorsions[index].c[4] = c4;
    rbTorsions[index].c[5] = c5;
}

void StandardMMForceField::getNonbonded14Parameters(int index, int& atom1, int& atom2, double& charge, double& radius, double& depth) const {
    atom1 = nb14s[index].atom1;
    atom2 = nb14s[index].atom2;
    charge = nb14s[index].charge;
    radius = nb14s[index].radius;
    depth = nb14s[index].depth;
}

void StandardMMForceField::setNonbonded14Parameters(int index, int atom1, int atom2, double charge, double radius, double depth) {
    nb14s[index].atom1 = atom1;
    nb14s[index].atom2 = atom2;
    nb14s[index].charge = charge;
    nb14s[index].radius = radius;
    nb14s[index].depth = depth;
}

ForceImpl* StandardMMForceField::createImpl() {
    return new StandardMMForceFieldImpl(*this);
}

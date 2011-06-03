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

#include "openmm/OpenMMException.h"
#include "openmm/State.h"

using namespace OpenMM;
using namespace std;

double State::getTime() const {
    return time;
}
const vector<Vec3>& State::getPositions() const {
    if ((types&Positions) == 0)
        throw OpenMMException("Invoked getPositions() on a State which does not contain positions.");
    return positions;
}
const vector<Vec3>& State::getVelocities() const {
    if ((types&Velocities) == 0)
        throw OpenMMException("Invoked getVelocities() on a State which does not contain velocities.");
    return velocities;
}
const vector<Vec3>& State::getForces() const {
    if ((types&Forces) == 0)
        throw OpenMMException("Invoked getForces() on a State which does not contain forces.");
    return forces;
}
double State::getKineticEnergy() const {
    if ((types&Energy) == 0)
        throw OpenMMException("Invoked getKineticEnergy() on a State which does not contain energies.");
    return ke;
}
double State::getPotentialEnergy() const {
    if ((types&Energy) == 0)
        throw OpenMMException("Invoked getPotentialEnergy() on a State which does not contain energies.");
    return pe;
}
void State::getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) const {
    a = periodicBoxVectors[0];
    b = periodicBoxVectors[1];
    c = periodicBoxVectors[2];
}
const map<string, double>& State::getParameters() const {
    if ((types&Parameters) == 0)
        throw OpenMMException("Invoked getParameters() on a State which does not contain parameters.");
    return parameters;
}
State::State(double time, int numParticles, int types) : types(types), time(time), ke(0), pe(0),
        positions( (types & Positions) == 0 ? 0 : numParticles), velocities( (types & Velocities) == 0 ? 0 : numParticles),
        forces( (types & Forces) == 0 ? 0 : numParticles) {
}
State::State() : types(0), time(0.0), ke(0), pe(0), positions(0), velocities(0), forces(0) {
}
vector<Vec3>& State::updPositions() {
    return positions;
}
vector<Vec3>& State::updVelocities() {
    return velocities;
}
vector<Vec3>& State::updForces() {
    return forces;
}
map<string, double>& State::updParameters() {
    return parameters;
}
void State::setEnergy(double kinetic, double potential) {
    ke = kinetic;
    pe = potential;
}

void State::setPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) {
    periodicBoxVectors[0] = a;
    periodicBoxVectors[1] = b;
    periodicBoxVectors[2] = c;
}

/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2019 Stanford University and the Authors.      *
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
#include "openmm/HippoNonbondedForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/HippoNonbondedForceImpl.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

HippoNonbondedForce::HippoNonbondedForce() : nonbondedMethod(NoCutoff), cutoffDistance(1.0), switchingDistance(0.9),
        ewaldErrorTol(1e-4), alpha(0.0), dalpha(0.0), nx(0), ny(0), nz(0), dnx(0), dny(0), dnz(0) {
    extrapolationCoefficients = {0.042, 0.635, 0.414};
}

HippoNonbondedForce::NonbondedMethod HippoNonbondedForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void HippoNonbondedForce::setNonbondedMethod(HippoNonbondedForce::NonbondedMethod method) {
    if (method < 0 || method > 1)
        throw OpenMMException("HippoNonbondedForce: Illegal value for nonbonded method");
    nonbondedMethod = method;
}

double HippoNonbondedForce::getCutoffDistance() const {
    return cutoffDistance;
}

void HippoNonbondedForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

double HippoNonbondedForce::getSwitchingDistance() const {
    return switchingDistance;
}

void HippoNonbondedForce::setSwitchingDistance(double distance) {
    switchingDistance = distance;
}

const std::vector<double> & HippoNonbondedForce::getExtrapolationCoefficients() const {
    return extrapolationCoefficients;
}

void HippoNonbondedForce::setExtrapolationCoefficients(const std::vector<double> &coefficients) {
    extrapolationCoefficients = coefficients;
}

void HippoNonbondedForce::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    alpha = this->alpha;
    nx = this->nx;
    ny = this->ny;
    nz = this->nz;
}

void HippoNonbondedForce::getDPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    alpha = this->dalpha;
    nx = this->dnx;
    ny = this->dny;
    nz = this->dnz;
}

void HippoNonbondedForce::setPMEParameters(double alpha, int nx, int ny, int nz) {
    this->alpha = alpha;
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
}

void HippoNonbondedForce::setDPMEParameters(double alpha, int nx, int ny, int nz) {
    this->dalpha = alpha;
    this->dnx = nx;
    this->dny = ny;
    this->dnz = nz;
}

void HippoNonbondedForce::getPMEParametersInContext(const Context& context, double& alpha, int& nx, int& ny, int& nz) const {
    dynamic_cast<const HippoNonbondedForceImpl&>(getImplInContext(context)).getPMEParameters(alpha, nx, ny, nz);
}

void HippoNonbondedForce::getDPMEParametersInContext(const Context& context, double& alpha, int& nx, int& ny, int& nz) const {
    dynamic_cast<const HippoNonbondedForceImpl&>(getImplInContext(context)).getDPMEParameters(alpha, nx, ny, nz);
}

double HippoNonbondedForce::getEwaldErrorTolerance() const {
    return ewaldErrorTol;
}

void HippoNonbondedForce::setEwaldErrorTolerance(double tol) {
    ewaldErrorTol = tol;
}

int HippoNonbondedForce::addParticle(double charge, const std::vector<double>& dipole, const std::vector<double>& quadrupole, double coreCharge,
                                     double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha,
                                     double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY) {
    particles.push_back(ParticleInfo(charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha,
                                     polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY));
    return particles.size()-1;
}

void HippoNonbondedForce::getParticleParameters(int index, double& charge, std::vector<double>& dipole, std::vector<double>& quadrupole, double& coreCharge,
                                                double& alpha, double& epsilon, double& damping, double& c6, double& pauliK, double& pauliQ, double& pauliAlpha,
                                                double& polarizability, int& axisType, int& multipoleAtomZ, int& multipoleAtomX, int& multipoleAtomY) const {
    charge = particles[index].charge;
    dipole = particles[index].dipole;
    quadrupole = particles[index].quadrupole;
    coreCharge = particles[index].coreCharge;
    alpha = particles[index].alpha;
    epsilon = particles[index].epsilon;
    damping = particles[index].damping;
    c6 = particles[index].c6;
    pauliK = particles[index].pauliK;
    pauliQ = particles[index].pauliQ;
    pauliAlpha = particles[index].pauliAlpha;
    polarizability = particles[index].polarizability;
    axisType = particles[index].axisType;
    multipoleAtomZ = particles[index].multipoleAtomZ;
    multipoleAtomX = particles[index].multipoleAtomX;
    multipoleAtomY = particles[index].multipoleAtomY;
}

void HippoNonbondedForce::setParticleParameters(int index, double charge, const std::vector<double>& dipole, const std::vector<double>& quadrupole, double coreCharge,
                                                double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha,
                                                double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY) {
    particles[index].charge = charge;
    particles[index].dipole = dipole;
    particles[index].quadrupole = quadrupole;
    particles[index].coreCharge = coreCharge;
    particles[index].alpha = alpha;
    particles[index].epsilon = epsilon;
    particles[index].damping = damping;
    particles[index].c6 = c6;
    particles[index].pauliK = pauliK;
    particles[index].pauliQ = pauliQ;
    particles[index].pauliAlpha = pauliAlpha;
    particles[index].polarizability = polarizability;
    particles[index].axisType = axisType;
    particles[index].multipoleAtomZ = multipoleAtomZ;
    particles[index].multipoleAtomX = multipoleAtomX;
    particles[index].multipoleAtomY = multipoleAtomY;
}

int HippoNonbondedForce::addException(int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale,
        double dispersionScale, double repulsionScale, double chargeTransferScale, bool replace) {
    map<pair<int, int>, int>::iterator iter = exceptionMap.find(pair<int, int>(particle1, particle2));
    int newIndex;
    if (iter == exceptionMap.end())
        iter = exceptionMap.find(pair<int, int>(particle2, particle1));
    if (iter != exceptionMap.end()) {
        if (!replace) {
            stringstream msg;
            msg << "HippoNonbondedForce: There is already an exception for particles ";
            msg << particle1;
            msg << " and ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        exceptions[iter->second] = ExceptionInfo(particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
        newIndex = iter->second;
        exceptionMap.erase(iter->first);
    }
    else {
        exceptions.push_back(ExceptionInfo(particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale));
        newIndex = exceptions.size()-1;
    }
    exceptionMap[pair<int, int>(particle1, particle2)] = newIndex;
    return newIndex;
}
void HippoNonbondedForce::getExceptionParameters(int index, int& particle1, int& particle2, double& multipoleMultipoleScale, double& dipoleMultipoleScale, double& dipoleDipoleScale,
        double& dispersionScale, double& repulsionScale, double& chargeTransferScale) const {
    ASSERT_VALID_INDEX(index, exceptions);
    particle1 = exceptions[index].particle1;
    particle2 = exceptions[index].particle2;
    multipoleMultipoleScale = exceptions[index].multipoleMultipoleScale;
    dipoleMultipoleScale = exceptions[index].dipoleMultipoleScale;
    dipoleDipoleScale = exceptions[index].dipoleDipoleScale;
    dispersionScale = exceptions[index].dispersionScale;
    repulsionScale = exceptions[index].repulsionScale;
    chargeTransferScale = exceptions[index].chargeTransferScale;
}

void HippoNonbondedForce::setExceptionParameters(int index, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale,
        double dispersionScale, double repulsionScale, double chargeTransferScale) {
    ASSERT_VALID_INDEX(index, exceptions);
    exceptions[index].particle1 = particle1;
    exceptions[index].particle2 = particle2;
    exceptions[index].multipoleMultipoleScale = multipoleMultipoleScale;
    exceptions[index].dipoleMultipoleScale = dipoleMultipoleScale;
    exceptions[index].dipoleDipoleScale = dipoleDipoleScale;
    exceptions[index].dispersionScale = dispersionScale;
    exceptions[index].repulsionScale = repulsionScale;
    exceptions[index].chargeTransferScale = chargeTransferScale;
}

void HippoNonbondedForce::getInducedDipoles(Context& context, vector<Vec3>& dipoles) {
    dynamic_cast<HippoNonbondedForceImpl&>(getImplInContext(context)).getInducedDipoles(getContextImpl(context), dipoles);
}

void HippoNonbondedForce::getLabFramePermanentDipoles(Context& context, vector<Vec3>& dipoles) {
    dynamic_cast<HippoNonbondedForceImpl&>(getImplInContext(context)).getLabFramePermanentDipoles(getContextImpl(context), dipoles);
}

ForceImpl* HippoNonbondedForce::createImpl()  const {
    return new HippoNonbondedForceImpl(*this);
}

void HippoNonbondedForce::updateParametersInContext(Context& context) {
    dynamic_cast<HippoNonbondedForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

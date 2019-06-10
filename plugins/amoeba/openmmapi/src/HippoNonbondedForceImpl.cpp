/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2018 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
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

#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/HippoNonbondedForceImpl.h"
#include "openmm/amoebaKernels.h"

using namespace OpenMM;
using namespace std;

HippoNonbondedForceImpl::HippoNonbondedForceImpl(const HippoNonbondedForce& owner) : owner(owner) {
}

HippoNonbondedForceImpl::~HippoNonbondedForceImpl() {
}

void HippoNonbondedForceImpl::initialize(ContextImpl& context) {
    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    if (owner.getNumParticles() != numParticles)
        throw OpenMMException("HippoNonbondedForce must have exactly as many particles as the System it belongs to.");

    // check cutoff < 0.5*boxSize

    if (owner.getNonbondedMethod() == HippoNonbondedForce::PME) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoffDistance();
        if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
            throw OpenMMException("HippoNonbondedForce: The cutoff distance cannot be greater than half the periodic box size.");
    }

    double quadrupoleValidationTolerance = 1.0e-05;
    for (int i = 0; i < numParticles; i++) {
        int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
        double charge, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability;
        vector<double> dipole, quadrupole;
        owner.getParticleParameters(i, charge, dipole, quadrupole, coreCharge,
                               alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha,
                               polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);

       // check quadrupole is traceless and symmetric

       double trace = fabs(quadrupole[0] + quadrupole[4] + quadrupole[8]);
       if (trace > quadrupoleValidationTolerance) {
             std::stringstream buffer;
             buffer << "HippoNonbondedForce: qudarupole for particle=" << i;
             buffer << " has nonzero trace: " << trace << "; AMOEBA plugin assumes traceless quadrupole.";
             throw OpenMMException(buffer.str());
       }
       if (fabs(quadrupole[1] - quadrupole[3]) > quadrupoleValidationTolerance ) {
             std::stringstream buffer;
             buffer << "HippoNonbondedForce: XY and YX components of quadrupole for particle=" << i;
             buffer << "  are not equal: [" << quadrupole[1] << " " << quadrupole[3] << "];";
             buffer << " AMOEBA plugin assumes symmetric quadrupole tensor.";
             throw OpenMMException(buffer.str());
       }
       if (fabs(quadrupole[2] - quadrupole[6]) > quadrupoleValidationTolerance ) {
             std::stringstream buffer;
             buffer << "HippoNonbondedForce: XZ and ZX components of quadrupole for particle=" << i;
             buffer << "  are not equal: [" << quadrupole[2] << " " << quadrupole[6] << "];";
             buffer << " AMOEBA plugin assumes symmetric quadrupole tensor.";
             throw OpenMMException(buffer.str());
       }
       if (fabs(quadrupole[5] - quadrupole[7]) > quadrupoleValidationTolerance ) {
             std::stringstream buffer;
             buffer << "HippoNonbondedForce: YZ and ZY components of quadrupole for particle=" << i;
             buffer << "  are not equal: [" << quadrupole[5] << " " << quadrupole[7] << "];";
             buffer << " AMOEBA plugin assumes symmetric quadrupole tensor.";
             throw OpenMMException(buffer.str());
       }

       // only 'Z-then-X', 'Bisector', Z-Bisect, ThreeFold  currently handled

        if (axisType != HippoNonbondedForce::ZThenX     && axisType != HippoNonbondedForce::Bisector &&
            axisType != HippoNonbondedForce::ZBisect    && axisType != HippoNonbondedForce::ThreeFold &&
            axisType != HippoNonbondedForce::ZOnly      && axisType != HippoNonbondedForce::NoAxisType) {
             std::stringstream buffer;
             buffer << "HippoNonbondedForce: axis type=" << axisType;
             buffer << " not currently handled - only axisTypes[ ";
             buffer << HippoNonbondedForce::ZThenX   << ", " << HippoNonbondedForce::Bisector  << ", ";
             buffer << HippoNonbondedForce::ZBisect  << ", " << HippoNonbondedForce::ThreeFold << ", ";
             buffer << HippoNonbondedForce::NoAxisType;
             buffer << "] (ZThenX, Bisector, Z-Bisect, ThreeFold, NoAxisType) currently handled .";
             throw OpenMMException(buffer.str());
        }
        if (axisType != HippoNonbondedForce::NoAxisType && (multipoleAtomZ < 0 || multipoleAtomZ >= numParticles)) {
            std::stringstream buffer;
            buffer << "HippoNonbondedForce: invalid z axis particle: " << multipoleAtomZ;
            throw OpenMMException(buffer.str());
        }
        if (axisType != HippoNonbondedForce::NoAxisType && axisType != HippoNonbondedForce::ZOnly &&
                (multipoleAtomX < 0 || multipoleAtomX >= numParticles)) {
            std::stringstream buffer;
            buffer << "HippoNonbondedForce: invalid x axis particle: " << multipoleAtomX;
            throw OpenMMException(buffer.str());
        }
        if ((axisType == HippoNonbondedForce::ZBisect || axisType == HippoNonbondedForce::ThreeFold) &&
                (multipoleAtomY < 0 || multipoleAtomY >= numParticles)) {
            std::stringstream buffer;
            buffer << "HippoNonbondedForce: invalid y axis particle: " << multipoleAtomY;
            throw OpenMMException(buffer.str());
        }
    }
    kernel = context.getPlatform().createKernel(CalcHippoNonbondedForceKernel::Name(), context);
    kernel.getAs<CalcHippoNonbondedForceKernel>().initialize(context.getSystem(), owner);
}

double HippoNonbondedForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcHippoNonbondedForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> HippoNonbondedForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcHippoNonbondedForceKernel::Name());
    return names;
}

void HippoNonbondedForceImpl::getLabFramePermanentDipoles(ContextImpl& context, vector<Vec3>& dipoles) {
    kernel.getAs<CalcHippoNonbondedForceKernel>().getLabFramePermanentDipoles(context, dipoles);
}

void HippoNonbondedForceImpl::getInducedDipoles(ContextImpl& context, vector<Vec3>& dipoles) {
    kernel.getAs<CalcHippoNonbondedForceKernel>().getInducedDipoles(context, dipoles);
}

void HippoNonbondedForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcHippoNonbondedForceKernel>().copyParametersToContext(context, owner);
    context.systemChanged();
}

void HippoNonbondedForceImpl::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    kernel.getAs<CalcHippoNonbondedForceKernel>().getPMEParameters(alpha, nx, ny, nz);
}

void HippoNonbondedForceImpl::getDPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    kernel.getAs<CalcHippoNonbondedForceKernel>().getDPMEParameters(alpha, nx, ny, nz);
}

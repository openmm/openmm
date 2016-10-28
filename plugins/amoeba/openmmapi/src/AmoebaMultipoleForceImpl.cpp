/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/AmoebaMultipoleForceImpl.h"
#include "openmm/amoebaKernels.h"
#include <stdio.h>

using namespace OpenMM;

using std::vector;

bool AmoebaMultipoleForceImpl::initializedCovalentDegrees = false;
int AmoebaMultipoleForceImpl::CovalentDegrees[]           = { 1,2,3,4,0,1,2,3};

AmoebaMultipoleForceImpl::AmoebaMultipoleForceImpl(const AmoebaMultipoleForce& owner) : owner(owner) {
}

AmoebaMultipoleForceImpl::~AmoebaMultipoleForceImpl() {
}

void AmoebaMultipoleForceImpl::initialize(ContextImpl& context) {

    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    if (owner.getNumMultipoles() != numParticles)
        throw OpenMMException("AmoebaMultipoleForce must have exactly as many particles as the System it belongs to.");

    // check cutoff < 0.5*boxSize

    if (owner.getNonbondedMethod() == AmoebaMultipoleForce::PME) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoffDistance();
        if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
            throw OpenMMException("AmoebaMultipoleForce: The cutoff distance cannot be greater than half the periodic box size.");
    }

    double quadrupoleValidationTolerance = 1.0e-05;
    for (int ii = 0; ii < numParticles; ii++) {

        int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
        double charge, thole, dampingFactor, polarity ;
        std::vector<double> molecularDipole;
        std::vector<double> molecularQuadrupole;

        owner.getMultipoleParameters(ii, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY,
                                     thole, dampingFactor, polarity);

       // check quadrupole is traceless and symmetric

       double trace = fabs(molecularQuadrupole[0] + molecularQuadrupole[4] + molecularQuadrupole[8]);
       if (trace > quadrupoleValidationTolerance) {
             std::stringstream buffer;
             buffer << "AmoebaMultipoleForce: qudarupole for particle=" << ii;
             buffer << " has nonzero trace: " << trace << "; AMOEBA plugin assumes traceless quadrupole.";
             throw OpenMMException(buffer.str());
       }
       if (fabs(molecularQuadrupole[1] - molecularQuadrupole[3]) > quadrupoleValidationTolerance ) {
             std::stringstream buffer;
             buffer << "AmoebaMultipoleForce: XY and YX components of quadrupole for particle=" << ii;
             buffer << "  are not equal: [" << molecularQuadrupole[1] << " " << molecularQuadrupole[3] << "];";
             buffer << " AMOEBA plugin assumes symmetric quadrupole tensor.";
             throw OpenMMException(buffer.str());
       }

       if (fabs(molecularQuadrupole[2] - molecularQuadrupole[6]) > quadrupoleValidationTolerance ) {
             std::stringstream buffer;
             buffer << "AmoebaMultipoleForce: XZ and ZX components of quadrupole for particle=" << ii;
             buffer << "  are not equal: [" << molecularQuadrupole[2] << " " << molecularQuadrupole[6] << "];";
             buffer << " AMOEBA plugin assumes symmetric quadrupole tensor.";
             throw OpenMMException(buffer.str());
       }

       if (fabs(molecularQuadrupole[5] - molecularQuadrupole[7]) > quadrupoleValidationTolerance ) {
             std::stringstream buffer;
             buffer << "AmoebaMultipoleForce: YZ and ZY components of quadrupole for particle=" << ii;
             buffer << "  are not equal: [" << molecularQuadrupole[5] << " " << molecularQuadrupole[7] << "];";
             buffer << " AMOEBA plugin assumes symmetric quadrupole tensor.";
             throw OpenMMException(buffer.str());
       }

       // only 'Z-then-X', 'Bisector', Z-Bisect, ThreeFold  currently handled

        if (axisType != AmoebaMultipoleForce::ZThenX     && axisType != AmoebaMultipoleForce::Bisector &&
            axisType != AmoebaMultipoleForce::ZBisect    && axisType != AmoebaMultipoleForce::ThreeFold &&
            axisType != AmoebaMultipoleForce::ZOnly      && axisType != AmoebaMultipoleForce::NoAxisType) {
             std::stringstream buffer;
             buffer << "AmoebaMultipoleForce: axis type=" << axisType;
             buffer << " not currently handled - only axisTypes[ ";
             buffer << AmoebaMultipoleForce::ZThenX   << ", " << AmoebaMultipoleForce::Bisector  << ", ";
             buffer << AmoebaMultipoleForce::ZBisect  << ", " << AmoebaMultipoleForce::ThreeFold << ", ";
             buffer << AmoebaMultipoleForce::NoAxisType;
             buffer << "] (ZThenX, Bisector, Z-Bisect, ThreeFold, NoAxisType) currently handled .";
             throw OpenMMException(buffer.str());
        }
        if (axisType != AmoebaMultipoleForce::NoAxisType && (multipoleAtomZ < 0 || multipoleAtomZ >= numParticles)) {
            std::stringstream buffer;
            buffer << "AmoebaMultipoleForce: invalid z axis particle: " << multipoleAtomZ;
            throw OpenMMException(buffer.str());
        }
        if (axisType != AmoebaMultipoleForce::NoAxisType && axisType != AmoebaMultipoleForce::ZOnly &&
                (multipoleAtomX < 0 || multipoleAtomX >= numParticles)) {
            std::stringstream buffer;
            buffer << "AmoebaMultipoleForce: invalid x axis particle: " << multipoleAtomX;
            throw OpenMMException(buffer.str());
        }
        if ((axisType == AmoebaMultipoleForce::ZBisect || axisType == AmoebaMultipoleForce::ThreeFold) &&
                (multipoleAtomY < 0 || multipoleAtomY >= numParticles)) {
            std::stringstream buffer;
            buffer << "AmoebaMultipoleForce: invalid y axis particle: " << multipoleAtomY;
            throw OpenMMException(buffer.str());
        }
    }
    kernel = context.getPlatform().createKernel(CalcAmoebaMultipoleForceKernel::Name(), context);
    kernel.getAs<CalcAmoebaMultipoleForceKernel>().initialize(context.getSystem(), owner);
}

double AmoebaMultipoleForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcAmoebaMultipoleForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> AmoebaMultipoleForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcAmoebaMultipoleForceKernel::Name());
    return names;
}

const int* AmoebaMultipoleForceImpl::getCovalentDegrees() {
    if (!initializedCovalentDegrees) {
        initializedCovalentDegrees                                      = true;
        CovalentDegrees[AmoebaMultipoleForce::Covalent12]               = 1;
        CovalentDegrees[AmoebaMultipoleForce::Covalent13]               = 2;
        CovalentDegrees[AmoebaMultipoleForce::Covalent14]               = 3;
        CovalentDegrees[AmoebaMultipoleForce::Covalent15]               = 4;
        CovalentDegrees[AmoebaMultipoleForce::PolarizationCovalent11]   = 0;
        CovalentDegrees[AmoebaMultipoleForce::PolarizationCovalent12]   = 1;
        CovalentDegrees[AmoebaMultipoleForce::PolarizationCovalent13]   = 2;
        CovalentDegrees[AmoebaMultipoleForce::PolarizationCovalent14]   = 3;
    }
    return CovalentDegrees;
}

void AmoebaMultipoleForceImpl::getCovalentRange(const AmoebaMultipoleForce& force, int atomIndex, const std::vector<AmoebaMultipoleForce::CovalentType>& lists,
                                                int* minCovalentIndex, int* maxCovalentIndex) {

    *minCovalentIndex =  999999999;
    *maxCovalentIndex = -999999999;
    for (unsigned int kk = 0; kk < lists.size(); kk++) {
        AmoebaMultipoleForce::CovalentType jj = lists[kk];
        std::vector<int> covalentList;
        force.getCovalentMap(atomIndex, jj, covalentList);
        for (unsigned int ii = 0; ii < covalentList.size(); ii++) {
            if (*minCovalentIndex > covalentList[ii]) {
               *minCovalentIndex = covalentList[ii];
            }
            if (*maxCovalentIndex < covalentList[ii]) {
               *maxCovalentIndex = covalentList[ii];
            }
        }
    }
    return;
}

void AmoebaMultipoleForceImpl::getCovalentDegree(const AmoebaMultipoleForce& force, std::vector<int>& covalentDegree) {
    covalentDegree.resize(AmoebaMultipoleForce::CovalentEnd);
    const int* CovalentDegrees = AmoebaMultipoleForceImpl::getCovalentDegrees();
    for (unsigned int kk = 0; kk < AmoebaMultipoleForce::CovalentEnd; kk++) {
        covalentDegree[kk] = CovalentDegrees[kk];
    }
    return;
}

void AmoebaMultipoleForceImpl::getLabFramePermanentDipoles(ContextImpl& context, vector<Vec3>& dipoles) {
    kernel.getAs<CalcAmoebaMultipoleForceKernel>().getLabFramePermanentDipoles(context, dipoles);
}

void AmoebaMultipoleForceImpl::getInducedDipoles(ContextImpl& context, vector<Vec3>& dipoles) {
    kernel.getAs<CalcAmoebaMultipoleForceKernel>().getInducedDipoles(context, dipoles);
}

void AmoebaMultipoleForceImpl::getTotalDipoles(ContextImpl& context, vector<Vec3>& dipoles) {
    kernel.getAs<CalcAmoebaMultipoleForceKernel>().getTotalDipoles(context, dipoles);
}

void AmoebaMultipoleForceImpl::getElectrostaticPotential(ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                                         std::vector< double >& outputElectrostaticPotential) {
    kernel.getAs<CalcAmoebaMultipoleForceKernel>().getElectrostaticPotential(context, inputGrid, outputElectrostaticPotential);
}

void AmoebaMultipoleForceImpl::getSystemMultipoleMoments(ContextImpl& context, std::vector< double >& outputMultipoleMoments) {
    kernel.getAs<CalcAmoebaMultipoleForceKernel>().getSystemMultipoleMoments(context, outputMultipoleMoments);
}

void AmoebaMultipoleForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcAmoebaMultipoleForceKernel>().copyParametersToContext(context, owner);
    context.systemChanged();
}

void AmoebaMultipoleForceImpl::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    kernel.getAs<CalcAmoebaMultipoleForceKernel>().getPMEParameters(alpha, nx, ny, nz);
}

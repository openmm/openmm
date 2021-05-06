/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2020 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "AmoebaCudaKernels.h"
#include "CudaAmoebaKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/AmoebaGeneralizedKirkwoodForceImpl.h"
#include "openmm/internal/AmoebaMultipoleForceImpl.h"
#include "openmm/internal/AmoebaWcaDispersionForceImpl.h"
#include "openmm/internal/AmoebaTorsionTorsionForceImpl.h"
#include "openmm/internal/AmoebaVdwForceImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "CudaBondedUtilities.h"
#include "CudaFFT3D.h"
#include "CudaForceInfo.h"
#include "CudaKernelSources.h"
#include "SimTKOpenMMRealType.h"
#include "jama_lu.h"

#include <algorithm>
#include <cmath>
#ifdef _MSC_VER
#include <windows.h>
#endif

using namespace OpenMM;
using namespace std;

#define CHECK_RESULT(result, prefix) \
    if (result != CUDA_SUCCESS) { \
        std::stringstream m; \
        m<<prefix<<": "<<cu.getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }

static void setPeriodicBoxArgs(ComputeContext& cc, ComputeKernel kernel, int index) {
    Vec3 a, b, c;
    cc.getPeriodicBoxVectors(a, b, c);
    if (cc.getUseDoublePrecision()) {
        kernel->setArg(index++, mm_double4(a[0], b[1], c[2], 0.0));
        kernel->setArg(index++, mm_double4(1.0/a[0], 1.0/b[1], 1.0/c[2], 0.0));
        kernel->setArg(index++, mm_double4(a[0], a[1], a[2], 0.0));
        kernel->setArg(index++, mm_double4(b[0], b[1], b[2], 0.0));
        kernel->setArg(index, mm_double4(c[0], c[1], c[2], 0.0));
    }
    else {
        kernel->setArg(index++, mm_float4((float) a[0], (float) b[1], (float) c[2], 0.0f));
        kernel->setArg(index++, mm_float4(1.0f/(float) a[0], 1.0f/(float) b[1], 1.0f/(float) c[2], 0.0f));
        kernel->setArg(index++, mm_float4((float) a[0], (float) a[1], (float) a[2], 0.0f));
        kernel->setArg(index++, mm_float4((float) b[0], (float) b[1], (float) b[2], 0.0f));
        kernel->setArg(index, mm_float4((float) c[0], (float) c[1], (float) c[2], 0.0f));
    }
}

/* -------------------------------------------------------------------------- *
 *                             AmoebaMultipole                                *
 * -------------------------------------------------------------------------- */

CudaCalcAmoebaMultipoleForceKernel::~CudaCalcAmoebaMultipoleForceKernel() {
    cc.setAsCurrent();
    if (hasInitializedFFT)
        cufftDestroy(fft);
}

void CudaCalcAmoebaMultipoleForceKernel::initialize(const System& system, const AmoebaMultipoleForce& force) {
    CommonCalcAmoebaMultipoleForceKernel::initialize(system, force);
    if (usePME) {
        cufftResult result = cufftPlan3d(&fft, gridSizeX, gridSizeY, gridSizeZ, cc.getUseDoublePrecision() ? CUFFT_Z2Z : CUFFT_C2C);
        if (result != CUFFT_SUCCESS)
            throw OpenMMException("Error initializing FFT: "+cc.intToString(result));
        hasInitializedFFT = true;
    }
}

void CudaCalcAmoebaMultipoleForceKernel::computeFFT(bool forward) {
    CudaArray& grid1 = dynamic_cast<CudaContext&>(cc).unwrap(pmeGrid1);
    CudaArray& grid2 = dynamic_cast<CudaContext&>(cc).unwrap(pmeGrid2);
    if (forward) {
        if (cc.getUseDoublePrecision())
            cufftExecZ2Z(fft, (double2*) grid1.getDevicePointer(), (double2*) grid2.getDevicePointer(), CUFFT_FORWARD);
        else
            cufftExecC2C(fft, (float2*) grid1.getDevicePointer(), (float2*) grid2.getDevicePointer(), CUFFT_FORWARD);
    }
    else {
        if (cc.getUseDoublePrecision())
            cufftExecZ2Z(fft, (double2*) grid2.getDevicePointer(), (double2*) grid1.getDevicePointer(), CUFFT_INVERSE);
        else
            cufftExecC2C(fft, (float2*) grid2.getDevicePointer(), (float2*) grid1.getDevicePointer(), CUFFT_INVERSE);
    }
}

/* -------------------------------------------------------------------------- *
 *                           HippoNonbondedForce                              *
 * -------------------------------------------------------------------------- */

class CudaCalcHippoNonbondedForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const HippoNonbondedForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, coreCharge1, alpha1, epsilon1, damping1, c61, pauliK1, pauliQ1, pauliAlpha1, polarizability1;
        double charge2, coreCharge2, alpha2, epsilon2, damping2, c62, pauliK2, pauliQ2, pauliAlpha2, polarizability2;
        int axisType1, multipoleZ1, multipoleX1, multipoleY1;
        int axisType2, multipoleZ2, multipoleX2, multipoleY2;
        vector<double> dipole1, dipole2, quadrupole1, quadrupole2;
        force.getParticleParameters(particle1, charge1, dipole1, quadrupole1, coreCharge1, alpha1, epsilon1, damping1, c61, pauliK1, pauliQ1, pauliAlpha1,
                                    polarizability1, axisType1, multipoleZ1, multipoleX1, multipoleY1);
        force.getParticleParameters(particle2, charge2, dipole2, quadrupole2, coreCharge2, alpha2, epsilon2, damping2, c62, pauliK2, pauliQ2, pauliAlpha2,
                                    polarizability2, axisType2, multipoleZ2, multipoleX2, multipoleY2);
        if (charge1 != charge2 || coreCharge1 != coreCharge2 || alpha1 != alpha2 || epsilon1 != epsilon1 || damping1 != damping2 || c61 != c62 ||
                pauliK1 != pauliK2 || pauliQ1 != pauliQ2 || pauliAlpha1 != pauliAlpha2 || polarizability1 != polarizability2 || axisType1 != axisType2) {
            return false;
        }
        for (int i = 0; i < dipole1.size(); ++i)
            if (dipole1[i] != dipole2[i])
                return false;
        for (int i = 0; i < quadrupole1.size(); ++i)
            if (quadrupole1[i] != quadrupole2[i])
                return false;
        return true;
    }
    int getNumParticleGroups() {
        return force.getNumExceptions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2;
        double multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale;
        force.getExceptionParameters(index, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        double multipoleMultipoleScale1, dipoleMultipoleScale1, dipoleDipoleScale1, dispersionScale1, repulsionScale1, chargeTransferScale1;
        double multipoleMultipoleScale2, dipoleMultipoleScale2, dipoleDipoleScale2, dispersionScale2, repulsionScale2, chargeTransferScale2;
        force.getExceptionParameters(group1, particle1, particle2, multipoleMultipoleScale1, dipoleMultipoleScale1, dipoleDipoleScale1, dispersionScale1, repulsionScale1, chargeTransferScale1);
        force.getExceptionParameters(group2, particle1, particle2, multipoleMultipoleScale2, dipoleMultipoleScale2, dipoleDipoleScale2, dispersionScale2, repulsionScale2, chargeTransferScale2);
        return (multipoleMultipoleScale1 == multipoleMultipoleScale2 && dipoleMultipoleScale1 == dipoleMultipoleScale2 &&
                dipoleDipoleScale1 == dipoleDipoleScale2 && dispersionScale1 == dispersionScale2 && repulsionScale1 == repulsionScale2 && chargeTransferScale1 == chargeTransferScale2);
    }
private:
    const HippoNonbondedForce& force;
};

class CudaCalcHippoNonbondedForceKernel::TorquePostComputation : public CudaContext::ForcePostComputation {
public:
    TorquePostComputation(CudaCalcHippoNonbondedForceKernel& owner) : owner(owner) {
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        owner.addTorquesToForces();
        return 0.0;
    }
private:
    CudaCalcHippoNonbondedForceKernel& owner;
};

CudaCalcHippoNonbondedForceKernel::CudaCalcHippoNonbondedForceKernel(const std::string& name, const Platform& platform, CudaContext& cu, const System& system) :
        CalcHippoNonbondedForceKernel(name, platform), cu(cu), system(system), sort(NULL), hasInitializedKernels(false), hasInitializedFFT(false), multipolesAreValid(false) {
}

CudaCalcHippoNonbondedForceKernel::~CudaCalcHippoNonbondedForceKernel() {
    cu.setAsCurrent();
    if (sort != NULL)
        delete sort;
    if (hasInitializedFFT) {
        cufftDestroy(fftForward);
        cufftDestroy(fftBackward);
        cufftDestroy(dfftForward);
        cufftDestroy(dfftBackward);
    }
}

void CudaCalcHippoNonbondedForceKernel::initialize(const System& system, const HippoNonbondedForce& force) {
    cu.setAsCurrent();
    extrapolationCoefficients = force.getExtrapolationCoefficients();
    usePME = (force.getNonbondedMethod() == HippoNonbondedForce::PME);

    // Initialize particle parameters.

    numParticles = force.getNumParticles();
    vector<double> coreChargeVec, valenceChargeVec, alphaVec, epsilonVec, dampingVec, c6Vec, pauliKVec, pauliQVec, pauliAlphaVec, polarizabilityVec;
    vector<double> localDipolesVec, localQuadrupolesVec;
    vector<int4> multipoleParticlesVec;
    vector<vector<int> > exclusions(numParticles);
    for (int i = 0; i < numParticles; i++) {
        double charge, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability;
        int axisType, atomX, atomY, atomZ;
        vector<double> dipole, quadrupole;
        force.getParticleParameters(i, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha,
                                    polarizability, axisType, atomZ, atomX, atomY);
        coreChargeVec.push_back(coreCharge);
        valenceChargeVec.push_back(charge-coreCharge);
        alphaVec.push_back(alpha);
        epsilonVec.push_back(epsilon);
        dampingVec.push_back(damping);
        c6Vec.push_back(c6);
        pauliKVec.push_back(pauliK);
        pauliQVec.push_back(pauliQ);
        pauliAlphaVec.push_back(pauliAlpha);
        polarizabilityVec.push_back(polarizability);
        multipoleParticlesVec.push_back(make_int4(atomX, atomY, atomZ, axisType));
        for (int j = 0; j < 3; j++)
            localDipolesVec.push_back(dipole[j]);
        localQuadrupolesVec.push_back(quadrupole[0]);
        localQuadrupolesVec.push_back(quadrupole[1]);
        localQuadrupolesVec.push_back(quadrupole[2]);
        localQuadrupolesVec.push_back(quadrupole[4]);
        localQuadrupolesVec.push_back(quadrupole[5]);
        exclusions[i].push_back(i);
    }
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    for (int i = numParticles; i < paddedNumAtoms; i++) {
        coreChargeVec.push_back(0);
        valenceChargeVec.push_back(0);
        alphaVec.push_back(0);
        epsilonVec.push_back(0);
        dampingVec.push_back(0);
        c6Vec.push_back(0);
        pauliKVec.push_back(0);
        pauliQVec.push_back(0);
        pauliAlphaVec.push_back(0);
        polarizabilityVec.push_back(0);
        multipoleParticlesVec.push_back(make_int4(0, 0, 0, 0));
        for (int j = 0; j < 3; j++)
            localDipolesVec.push_back(0);
        for (int j = 0; j < 5; j++)
            localQuadrupolesVec.push_back(0);
    }
    int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    coreCharge.initialize(cu, paddedNumAtoms, elementSize, "coreCharge");
    valenceCharge.initialize(cu, paddedNumAtoms, elementSize, "valenceCharge");
    alpha.initialize(cu, paddedNumAtoms, elementSize, "alpha");
    epsilon.initialize(cu, paddedNumAtoms, elementSize, "epsilon");
    damping.initialize(cu, paddedNumAtoms, elementSize, "damping");
    c6.initialize(cu, paddedNumAtoms, elementSize, "c6");
    pauliK.initialize(cu, paddedNumAtoms, elementSize, "pauliK");
    pauliQ.initialize(cu, paddedNumAtoms, elementSize, "pauliQ");
    pauliAlpha.initialize(cu, paddedNumAtoms, elementSize, "pauliAlpha");
    polarizability.initialize(cu, paddedNumAtoms, elementSize, "polarizability");
    multipoleParticles.initialize<int4>(cu, paddedNumAtoms, "multipoleParticles");
    localDipoles.initialize(cu, 3*paddedNumAtoms, elementSize, "localDipoles");
    localQuadrupoles.initialize(cu, 5*paddedNumAtoms, elementSize, "localQuadrupoles");
    lastPositions.initialize(cu, cu.getPosq().getSize(), cu.getPosq().getElementSize(), "lastPositions");
    coreCharge.upload(coreChargeVec, true);
    valenceCharge.upload(valenceChargeVec, true);
    alpha.upload(alphaVec, true);
    epsilon.upload(epsilonVec, true);
    damping.upload(dampingVec, true);
    c6.upload(c6Vec, true);
    pauliK.upload(pauliKVec, true);
    pauliQ.upload(pauliQVec, true);
    pauliAlpha.upload(pauliAlphaVec, true);
    polarizability.upload(polarizabilityVec, true);
    multipoleParticles.upload(multipoleParticlesVec);
    localDipoles.upload(localDipolesVec, true);
    localQuadrupoles.upload(localQuadrupolesVec, true);
    
    // Create workspace arrays.
    
    labDipoles.initialize(cu, paddedNumAtoms, 3*elementSize, "dipole");
    labQuadrupoles[0].initialize(cu, paddedNumAtoms, elementSize, "qXX");
    labQuadrupoles[1].initialize(cu, paddedNumAtoms, elementSize, "qXY");
    labQuadrupoles[2].initialize(cu, paddedNumAtoms, elementSize, "qXZ");
    labQuadrupoles[3].initialize(cu, paddedNumAtoms, elementSize, "qYY");
    labQuadrupoles[4].initialize(cu, paddedNumAtoms, elementSize, "qYZ");
    fracDipoles.initialize(cu, paddedNumAtoms, 3*elementSize, "fracDipoles");
    fracQuadrupoles.initialize(cu, 6*paddedNumAtoms, elementSize, "fracQuadrupoles");
    field.initialize(cu, 3*paddedNumAtoms, sizeof(long long), "field");
    inducedField.initialize(cu, 3*paddedNumAtoms, sizeof(long long), "inducedField");
    torque.initialize(cu, 3*paddedNumAtoms, sizeof(long long), "torque");
    inducedDipole.initialize(cu, paddedNumAtoms, 3*elementSize, "inducedDipole");
    int numOrders = extrapolationCoefficients.size();
    extrapolatedDipole.initialize(cu, 3*numParticles*numOrders, elementSize, "extrapolatedDipole");
    extrapolatedPhi.initialize(cu, 10*numParticles*numOrders, elementSize, "extrapolatedPhi");
    cu.addAutoclearBuffer(field);
    cu.addAutoclearBuffer(torque);
    
    // Record exceptions and exclusions.
    
    vector<double> exceptionScaleVec[6];
    vector<int2> exceptionAtomsVec;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale;
        force.getExceptionParameters(i, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
        exclusions[particle1].push_back(particle2);
        exclusions[particle2].push_back(particle1);
        if (usePME || multipoleMultipoleScale != 0 || dipoleMultipoleScale != 0 || dipoleDipoleScale != 0 || dispersionScale != 0 || repulsionScale != 0 || chargeTransferScale != 0) {
            exceptionAtomsVec.push_back(make_int2(particle1, particle2));
            exceptionScaleVec[0].push_back(multipoleMultipoleScale);
            exceptionScaleVec[1].push_back(dipoleMultipoleScale);
            exceptionScaleVec[2].push_back(dipoleDipoleScale);
            exceptionScaleVec[3].push_back(dispersionScale);
            exceptionScaleVec[4].push_back(repulsionScale);
            exceptionScaleVec[5].push_back(chargeTransferScale);
        }
    }
    if (exceptionAtomsVec.size() > 0) {
        exceptionAtoms.initialize<int2>(cu, exceptionAtomsVec.size(), "exceptionAtoms");
        exceptionAtoms.upload(exceptionAtomsVec);
        for (int i = 0; i < 6; i++) {
            exceptionScales[i].initialize(cu, exceptionAtomsVec.size(), elementSize, "exceptionScales");
            exceptionScales[i].upload(exceptionScaleVec[i], true);
        }
    }
    
    // Create the kernels.

    bool useShuffle = (cu.getComputeCapability() >= 3.0 && !cu.getUseDoublePrecision());
    map<string, string> defines;
    defines["HIPPO"] = "1";
    defines["NUM_ATOMS"] = cu.intToString(numParticles);
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    defines["NUM_BLOCKS"] = cu.intToString(cu.getNumAtomBlocks());
    defines["ENERGY_SCALE_FACTOR"] = cu.doubleToString(ONE_4PI_EPS0);
    if (useShuffle)
        defines["USE_SHUFFLE"] = "";
    maxExtrapolationOrder = extrapolationCoefficients.size();
    defines["MAX_EXTRAPOLATION_ORDER"] = cu.intToString(maxExtrapolationOrder);
    stringstream coefficients;
    for (int i = 0; i < maxExtrapolationOrder; i++) {
        if (i > 0)
            coefficients << ",";
        double sum = 0;
        for (int j = i; j < maxExtrapolationOrder; j++)
            sum += extrapolationCoefficients[j];
        coefficients << cu.doubleToString(sum);
    }
    defines["EXTRAPOLATION_COEFFICIENTS_SUM"] = coefficients.str();
    cutoff = force.getCutoffDistance();
    if (usePME) {
        int nx, ny, nz;
        force.getPMEParameters(pmeAlpha, nx, ny, nz);
        if (nx == 0 || pmeAlpha == 0) {
            NonbondedForce nb;
            nb.setEwaldErrorTolerance(force.getEwaldErrorTolerance());
            nb.setCutoffDistance(force.getCutoffDistance());
            NonbondedForceImpl::calcPMEParameters(system, nb, pmeAlpha, gridSizeX, gridSizeY, gridSizeZ, false);
            gridSizeX = CudaFFT3D::findLegalDimension(gridSizeX);
            gridSizeY = CudaFFT3D::findLegalDimension(gridSizeY);
            gridSizeZ = CudaFFT3D::findLegalDimension(gridSizeZ);
        } else {
            gridSizeX = CudaFFT3D::findLegalDimension(nx);
            gridSizeY = CudaFFT3D::findLegalDimension(ny);
            gridSizeZ = CudaFFT3D::findLegalDimension(nz);
        }
        force.getDPMEParameters(dpmeAlpha, nx, ny, nz);
        if (nx == 0 || dpmeAlpha == 0) {
            NonbondedForce nb;
            nb.setEwaldErrorTolerance(force.getEwaldErrorTolerance());
            nb.setCutoffDistance(force.getCutoffDistance());
            NonbondedForceImpl::calcPMEParameters(system, nb, dpmeAlpha, dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ, true);
            dispersionGridSizeX = CudaFFT3D::findLegalDimension(dispersionGridSizeX);
            dispersionGridSizeY = CudaFFT3D::findLegalDimension(dispersionGridSizeY);
            dispersionGridSizeZ = CudaFFT3D::findLegalDimension(dispersionGridSizeZ);
        } else {
            dispersionGridSizeX = CudaFFT3D::findLegalDimension(nx);
            dispersionGridSizeY = CudaFFT3D::findLegalDimension(ny);
            dispersionGridSizeZ = CudaFFT3D::findLegalDimension(nz);
        }
        defines["EWALD_ALPHA"] = cu.doubleToString(pmeAlpha);
        defines["SQRT_PI"] = cu.doubleToString(sqrt(M_PI));
        defines["USE_EWALD"] = "";
        defines["USE_CUTOFF"] = "";
        defines["USE_PERIODIC"] = "";
        defines["CUTOFF_SQUARED"] = cu.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
    }
    ComputeProgram program = cu.compileProgram(CudaAmoebaKernelSources::hippoMultipoles, defines);
    computeMomentsKernel = program->createKernel("computeLabFrameMoments");
    computeMomentsKernel->addArg(cu.getPosq());
    computeMomentsKernel->addArg(multipoleParticles);
    computeMomentsKernel->addArg(localDipoles);
    computeMomentsKernel->addArg(localQuadrupoles);
    computeMomentsKernel->addArg(labDipoles);
    computeMomentsKernel->addArg(labQuadrupoles[0]);
    computeMomentsKernel->addArg(labQuadrupoles[1]);
    computeMomentsKernel->addArg(labQuadrupoles[2]);
    computeMomentsKernel->addArg(labQuadrupoles[3]);
    computeMomentsKernel->addArg(labQuadrupoles[4]);
    recordInducedDipolesKernel = program->createKernel("recordInducedDipoles");
    recordInducedDipolesKernel->addArg(field);
    recordInducedDipolesKernel->addArg(inducedDipole);
    recordInducedDipolesKernel->addArg(polarizability);
    mapTorqueKernel = program->createKernel("mapTorqueToForce");
    mapTorqueKernel->addArg(cu.getLongForceBuffer());
    mapTorqueKernel->addArg(torque);
    mapTorqueKernel->addArg(cu.getPosq());
    mapTorqueKernel->addArg(multipoleParticles);
    program = cu.compileProgram(CudaAmoebaKernelSources::multipoleInducedField, defines);
    initExtrapolatedKernel = program->createKernel("initExtrapolatedDipoles");
    initExtrapolatedKernel->addArg(inducedDipole);
    initExtrapolatedKernel->addArg(extrapolatedDipole);
    iterateExtrapolatedKernel = program->createKernel("iterateExtrapolatedDipoles");
    iterateExtrapolatedKernel->addArg();
    iterateExtrapolatedKernel->addArg(inducedDipole);
    iterateExtrapolatedKernel->addArg(extrapolatedDipole);
    iterateExtrapolatedKernel->addArg(inducedField);
    iterateExtrapolatedKernel->addArg(polarizability);
    computeExtrapolatedKernel = program->createKernel("computeExtrapolatedDipoles");
    computeExtrapolatedKernel->addArg(inducedDipole);
    computeExtrapolatedKernel->addArg(extrapolatedDipole);
    polarizationEnergyKernel = program->createKernel("computePolarizationEnergy");
    polarizationEnergyKernel->addArg(cu.getEnergyBuffer());
    polarizationEnergyKernel->addArg(inducedDipole);
    polarizationEnergyKernel->addArg(extrapolatedDipole);
    polarizationEnergyKernel->addArg(polarizability);

    // Set up PME.
    
    if (usePME) {
        // Create the PME kernels.

        map<string, string> pmeDefines;
        pmeDefines["HIPPO"] = "1";
        pmeDefines["EWALD_ALPHA"] = cu.doubleToString(pmeAlpha);
        pmeDefines["DISPERSION_EWALD_ALPHA"] = cu.doubleToString(dpmeAlpha);
        pmeDefines["PME_ORDER"] = cu.intToString(PmeOrder);
        pmeDefines["NUM_ATOMS"] = cu.intToString(numParticles);
        pmeDefines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
        pmeDefines["EPSILON_FACTOR"] = cu.doubleToString(ONE_4PI_EPS0);
        pmeDefines["GRID_SIZE_X"] = cu.intToString(gridSizeX);
        pmeDefines["GRID_SIZE_Y"] = cu.intToString(gridSizeY);
        pmeDefines["GRID_SIZE_Z"] = cu.intToString(gridSizeZ);
        pmeDefines["M_PI"] = cu.doubleToString(M_PI);
        pmeDefines["SQRT_PI"] = cu.doubleToString(sqrt(M_PI));
        pmeDefines["EXTRAPOLATION_COEFFICIENTS_SUM"] = coefficients.str();
        pmeDefines["MAX_EXTRAPOLATION_ORDER"] = cu.intToString(maxExtrapolationOrder);
        program = cu.compileProgram(CudaAmoebaKernelSources::multipolePme, pmeDefines);
        pmeTransformMultipolesKernel = program->createKernel("transformMultipolesToFractionalCoordinates");
        pmeTransformMultipolesKernel->addArg(labDipoles);
        pmeTransformMultipolesKernel->addArg(labQuadrupoles[0]);
        pmeTransformMultipolesKernel->addArg(labQuadrupoles[1]);
        pmeTransformMultipolesKernel->addArg(labQuadrupoles[2]);
        pmeTransformMultipolesKernel->addArg(labQuadrupoles[3]);
        pmeTransformMultipolesKernel->addArg(labQuadrupoles[4]);
        pmeTransformMultipolesKernel->addArg(fracDipoles);
        pmeTransformMultipolesKernel->addArg(fracQuadrupoles);
        for (int i = 0; i < 3; i++)
            pmeTransformMultipolesKernel->addArg();
        pmeTransformPotentialKernel = program->createKernel("transformPotentialToCartesianCoordinates");
        pmeTransformPotentialKernel->addArg();
        pmeTransformPotentialKernel->addArg(pmeCphi);
        for (int i = 0; i < 3; i++)
            pmeTransformPotentialKernel->addArg();
        pmeSpreadFixedMultipolesKernel = program->createKernel("gridSpreadFixedMultipoles");
        pmeSpreadFixedMultipolesKernel->addArg(cu.getPosq());
        pmeSpreadFixedMultipolesKernel->addArg(fracDipoles);
        pmeSpreadFixedMultipolesKernel->addArg(fracQuadrupoles);
        pmeSpreadFixedMultipolesKernel->addArg(pmeGrid1);
        pmeSpreadFixedMultipolesKernel->addArg(coreCharge);
        pmeSpreadFixedMultipolesKernel->addArg(valenceCharge);
        for (int i = 0; i < 6; i++)
            pmeSpreadFixedMultipolesKernel->addArg();
        pmeSpreadInducedDipolesKernel = program->createKernel("gridSpreadInducedDipoles");
        pmeSpreadInducedDipolesKernel->addArg(cu.getPosq());
        pmeSpreadInducedDipolesKernel->addArg(inducedDipole);
        pmeSpreadInducedDipolesKernel->addArg(pmeGrid1);
        for (int i = 0; i < 6; i++)
            pmeSpreadInducedDipolesKernel->addArg();
        pmeFinishSpreadChargeKernel = program->createKernel("finishSpreadCharge");
        pmeFinishSpreadChargeKernel->addArg(pmeGrid1);
        pmeConvolutionKernel = program->createKernel("reciprocalConvolution");
        pmeConvolutionKernel->addArg(pmeGrid2);
        pmeConvolutionKernel->addArg(pmeBsplineModuliX);
        pmeConvolutionKernel->addArg(pmeBsplineModuliY);
        pmeConvolutionKernel->addArg(pmeBsplineModuliX);
        for (int i = 0; i < 4; i++)
            pmeConvolutionKernel->addArg();
        pmeFixedPotentialKernel = program->createKernel("computeFixedPotentialFromGrid");
        pmeFixedPotentialKernel->addArg(pmeGrid1);
        pmeFixedPotentialKernel->addArg(pmePhi);
        pmeFixedPotentialKernel->addArg(field);
        pmeFixedPotentialKernel->addArg(cu.getPosq());
        pmeFixedPotentialKernel->addArg(labDipoles);
        for (int i = 0; i < 6; i++)
            pmeFixedPotentialKernel->addArg();
        pmeInducedPotentialKernel = program->createKernel("computeInducedPotentialFromGrid");
        pmeInducedPotentialKernel->addArg(pmeGrid1);
        pmeInducedPotentialKernel->addArg(extrapolatedPhi);
        pmeInducedPotentialKernel->addArg();
        pmeInducedPotentialKernel->addArg(pmePhidp);
        pmeInducedPotentialKernel->addArg(cu.getPosq());
        for (int i = 0; i < 6; i++)
            pmeInducedPotentialKernel->addArg();
        pmeFixedForceKernel = program->createKernel("computeFixedMultipoleForceAndEnergy");
        pmeFixedForceKernel->addArg(cu.getPosq());
        pmeFixedForceKernel->addArg(cu.getLongForceBuffer());
        pmeFixedForceKernel->addArg(torque);
        pmeFixedForceKernel->addArg(cu.getEnergyBuffer());
        pmeFixedForceKernel->addArg(labDipoles);
        pmeFixedForceKernel->addArg(coreCharge);
        pmeFixedForceKernel->addArg(valenceCharge);
        pmeFixedForceKernel->addArg(labQuadrupoles[0]);
        pmeFixedForceKernel->addArg(labQuadrupoles[1]);
        pmeFixedForceKernel->addArg(labQuadrupoles[2]);
        pmeFixedForceKernel->addArg(labQuadrupoles[3]);
        pmeFixedForceKernel->addArg(labQuadrupoles[4]);
        pmeFixedForceKernel->addArg(fracDipoles);
        pmeFixedForceKernel->addArg(fracQuadrupoles);
        pmeFixedForceKernel->addArg(pmePhi);
        pmeFixedForceKernel->addArg(pmeCphi);
        for (int i = 0; i < 3; i++)
            pmeFixedForceKernel->addArg();
        pmeInducedForceKernel = program->createKernel("computeInducedDipoleForceAndEnergy");
        pmeInducedForceKernel->addArg(cu.getPosq());
        pmeInducedForceKernel->addArg(cu.getLongForceBuffer());
        pmeInducedForceKernel->addArg(torque);
        pmeInducedForceKernel->addArg(cu.getEnergyBuffer());
        pmeInducedForceKernel->addArg(labDipoles);
        pmeInducedForceKernel->addArg(coreCharge);
        pmeInducedForceKernel->addArg(valenceCharge);
        pmeInducedForceKernel->addArg(extrapolatedDipole);
        pmeInducedForceKernel->addArg(extrapolatedPhi);
        pmeInducedForceKernel->addArg(labQuadrupoles[0]);
        pmeInducedForceKernel->addArg(labQuadrupoles[1]);
        pmeInducedForceKernel->addArg(labQuadrupoles[2]);
        pmeInducedForceKernel->addArg(labQuadrupoles[3]);
        pmeInducedForceKernel->addArg(labQuadrupoles[4]);
        pmeInducedForceKernel->addArg(fracDipoles);
        pmeInducedForceKernel->addArg(fracQuadrupoles);
        pmeInducedForceKernel->addArg(inducedDipole);
        pmeInducedForceKernel->addArg(pmePhi);
        pmeInducedForceKernel->addArg(pmePhidp);
        pmeInducedForceKernel->addArg(pmeCphi);
        for (int i = 0; i < 3; i++)
            pmeInducedForceKernel->addArg();
        pmeRecordInducedFieldDipolesKernel = program->createKernel("recordInducedFieldDipoles");
        pmeRecordInducedFieldDipolesKernel->addArg(pmePhidp);
        pmeRecordInducedFieldDipolesKernel->addArg(inducedField);
        pmeRecordInducedFieldDipolesKernel->addArg(inducedDipole);
        for (int i = 0; i < 3; i++)
            pmeRecordInducedFieldDipolesKernel->addArg();
        pmeSelfEnergyKernel = program->createKernel("calculateSelfEnergyAndTorque");
        pmeSelfEnergyKernel->addArg(torque);
        pmeSelfEnergyKernel->addArg(cu.getEnergyBuffer());
        pmeSelfEnergyKernel->addArg(labDipoles);
        pmeSelfEnergyKernel->addArg(coreCharge);
        pmeSelfEnergyKernel->addArg(valenceCharge);
        pmeSelfEnergyKernel->addArg(c6);
        pmeSelfEnergyKernel->addArg(inducedDipole);
        pmeSelfEnergyKernel->addArg(labQuadrupoles[0]);
        pmeSelfEnergyKernel->addArg(labQuadrupoles[1]);
        pmeSelfEnergyKernel->addArg(labQuadrupoles[2]);
        pmeSelfEnergyKernel->addArg(labQuadrupoles[3]);
        pmeSelfEnergyKernel->addArg(labQuadrupoles[4]);

        // Create the dispersion PME kernels.

        pmeDefines["EWALD_ALPHA"] = cu.doubleToString(dpmeAlpha);
        pmeDefines["EPSILON_FACTOR"] = "1";
        pmeDefines["GRID_SIZE_X"] = cu.intToString(dispersionGridSizeX);
        pmeDefines["GRID_SIZE_Y"] = cu.intToString(dispersionGridSizeY);
        pmeDefines["GRID_SIZE_Z"] = cu.intToString(dispersionGridSizeZ);
        pmeDefines["RECIP_EXP_FACTOR"] = cu.doubleToString(M_PI*M_PI/(dpmeAlpha*dpmeAlpha));
        pmeDefines["CHARGE"] = "charges[atom]";
        pmeDefines["USE_LJPME"] = "1";
        program = cu.compileProgram(CudaKernelSources::pme, pmeDefines);
        dpmeFinishSpreadChargeKernel = program->createKernel("finishSpreadCharge");
        dpmeFinishSpreadChargeKernel->addArg(pmeGrid2);
        dpmeFinishSpreadChargeKernel->addArg(pmeGrid1);
        dpmeGridIndexKernel = program->createKernel("findAtomGridIndex");
        dpmeGridIndexKernel->addArg(cu.getPosq());
        dpmeGridIndexKernel->addArg(pmeAtomGridIndex);
        for (int i = 0; i < 8; i++)
            dpmeGridIndexKernel->addArg();
        dpmeSpreadChargeKernel = program->createKernel("gridSpreadCharge");
        dpmeSpreadChargeKernel->addArg(cu.getPosq());
        dpmeSpreadChargeKernel->addArg(pmeGrid2);
        for (int i = 0; i < 8; i++)
            dpmeSpreadChargeKernel->addArg();
        dpmeSpreadChargeKernel->addArg(pmeAtomGridIndex);
        dpmeSpreadChargeKernel->addArg(c6);
        dpmeConvolutionKernel = program->createKernel("reciprocalConvolution");
        dpmeConvolutionKernel->addArg(pmeGrid2);
        dpmeConvolutionKernel->addArg(cu.getEnergyBuffer());
        dpmeConvolutionKernel->addArg(dpmeBsplineModuliX);
        dpmeConvolutionKernel->addArg(dpmeBsplineModuliY);
        dpmeConvolutionKernel->addArg(dpmeBsplineModuliZ);
        for (int i = 0; i < 4; i++)
            dpmeConvolutionKernel->addArg();
        dpmeEvalEnergyKernel = program->createKernel("gridEvaluateEnergy");
        dpmeEvalEnergyKernel->addArg(pmeGrid2);
        dpmeEvalEnergyKernel->addArg(cu.getEnergyBuffer());
        dpmeEvalEnergyKernel->addArg(dpmeBsplineModuliX);
        dpmeEvalEnergyKernel->addArg(dpmeBsplineModuliY);
        dpmeEvalEnergyKernel->addArg(dpmeBsplineModuliZ);
        for (int i = 0; i < 4; i++)
            dpmeEvalEnergyKernel->addArg();
        dpmeInterpolateForceKernel = program->createKernel("gridInterpolateForce");
        dpmeInterpolateForceKernel->addArg(cu.getPosq());
        dpmeInterpolateForceKernel->addArg(cu.getLongForceBuffer());
        dpmeInterpolateForceKernel->addArg(pmeGrid1);
        for (int i = 0; i < 8; i++)
            dpmeInterpolateForceKernel->addArg();
        dpmeInterpolateForceKernel->addArg(pmeAtomGridIndex);
        dpmeInterpolateForceKernel->addArg(c6);

        // Create required data structures.

        int roundedZSize = PmeOrder*(int) ceil(gridSizeZ/(double) PmeOrder);
        int gridElements = gridSizeX*gridSizeY*roundedZSize;
        roundedZSize = PmeOrder*(int) ceil(dispersionGridSizeZ/(double) PmeOrder);
        gridElements = max(gridElements, dispersionGridSizeX*dispersionGridSizeY*roundedZSize);
        pmeGrid1.initialize(cu, gridElements, elementSize, "pmeGrid1");
        pmeGrid2.initialize(cu, gridElements, 2*elementSize, "pmeGrid2");
        cu.addAutoclearBuffer(pmeGrid1);
        pmeBsplineModuliX.initialize(cu, gridSizeX, elementSize, "pmeBsplineModuliX");
        pmeBsplineModuliY.initialize(cu, gridSizeY, elementSize, "pmeBsplineModuliY");
        pmeBsplineModuliZ.initialize(cu, gridSizeZ, elementSize, "pmeBsplineModuliZ");
        dpmeBsplineModuliX.initialize(cu, dispersionGridSizeX, elementSize, "dpmeBsplineModuliX");
        dpmeBsplineModuliY.initialize(cu, dispersionGridSizeY, elementSize, "dpmeBsplineModuliY");
        dpmeBsplineModuliZ.initialize(cu, dispersionGridSizeZ, elementSize, "dpmeBsplineModuliZ");
        pmePhi.initialize(cu, 20*numParticles, elementSize, "pmePhi");
        pmePhidp.initialize(cu, 20*numParticles, elementSize, "pmePhidp");
        pmeCphi.initialize(cu, 10*numParticles, elementSize, "pmeCphi");
        pmeAtomGridIndex.initialize<int2>(cu, numParticles, "pmeAtomGridIndex");
        sort = new CudaSort(cu, new SortTrait(), cu.getNumAtoms());
        cufftResult result = cufftPlan3d(&fftForward, gridSizeX, gridSizeY, gridSizeZ, cu.getUseDoublePrecision() ? CUFFT_D2Z : CUFFT_R2C);
        if (result != CUFFT_SUCCESS)
            throw OpenMMException("Error initializing FFT: "+cu.intToString(result));
        result = cufftPlan3d(&fftBackward, gridSizeX, gridSizeY, gridSizeZ, cu.getUseDoublePrecision() ? CUFFT_Z2D : CUFFT_C2R);
        if (result != CUFFT_SUCCESS)
            throw OpenMMException("Error initializing FFT: "+cu.intToString(result));
        result = cufftPlan3d(&dfftForward, dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ, cu.getUseDoublePrecision() ? CUFFT_D2Z : CUFFT_R2C);
        if (result != CUFFT_SUCCESS)
            throw OpenMMException("Error initializing FFT: "+cu.intToString(result));
        result = cufftPlan3d(&dfftBackward, dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ, cu.getUseDoublePrecision() ? CUFFT_Z2D : CUFFT_C2R);
        if (result != CUFFT_SUCCESS)
            throw OpenMMException("Error initializing FFT: "+cu.intToString(result));
        hasInitializedFFT = true;

        // Initialize the B-spline moduli.

        double data[PmeOrder];
        double x = 0.0;
        data[0] = 1.0 - x;
        data[1] = x;
        for (int i = 2; i < PmeOrder; i++) {
            double denom = 1.0/i;
            data[i] = x*data[i-1]*denom;
            for (int j = 1; j < i; j++)
                data[i-j] = ((x+j)*data[i-j-1] + ((i-j+1)-x)*data[i-j])*denom;
            data[0] = (1.0-x)*data[0]*denom;
        }
        int maxSize = max(max(gridSizeX, gridSizeY), gridSizeZ);
        vector<double> bsplines_data(maxSize+1, 0.0);
        for (int i = 2; i <= PmeOrder+1; i++)
            bsplines_data[i] = data[i-2];
        for (int dim = 0; dim < 3; dim++) {
            int ndata = (dim == 0 ? gridSizeX : dim == 1 ? gridSizeY : gridSizeZ);
            vector<double> moduli(ndata);

            // get the modulus of the discrete Fourier transform

            double factor = 2.0*M_PI/ndata;
            for (int i = 0; i < ndata; i++) {
                double sc = 0.0;
                double ss = 0.0;
                for (int j = 1; j <= ndata; j++) {
                    double arg = factor*i*(j-1);
                    sc += bsplines_data[j]*cos(arg);
                    ss += bsplines_data[j]*sin(arg);
                }
                moduli[i] = sc*sc+ss*ss;
            }

            // Fix for exponential Euler spline interpolation failure.

            double eps = 1.0e-7;
            if (moduli[0] < eps)
                moduli[0] = 0.9*moduli[1];
            for (int i = 1; i < ndata-1; i++)
                if (moduli[i] < eps)
                    moduli[i] = 0.9*(moduli[i-1]+moduli[i+1]);
            if (moduli[ndata-1] < eps)
                moduli[ndata-1] = 0.9*moduli[ndata-2];

            // Compute and apply the optimal zeta coefficient.

            int jcut = 50;
            for (int i = 1; i <= ndata; i++) {
                int k = i - 1;
                if (i > ndata/2)
                    k = k - ndata;
                double zeta;
                if (k == 0)
                    zeta = 1.0;
                else {
                    double sum1 = 1.0;
                    double sum2 = 1.0;
                    factor = M_PI*k/ndata;
                    for (int j = 1; j <= jcut; j++) {
                        double arg = factor/(factor+M_PI*j);
                        sum1 += pow(arg, PmeOrder);
                        sum2 += pow(arg, 2*PmeOrder);
                    }
                    for (int j = 1; j <= jcut; j++) {
                        double arg = factor/(factor-M_PI*j);
                        sum1 += pow(arg, PmeOrder);
                        sum2 += pow(arg, 2*PmeOrder);
                    }
                    zeta = sum2/sum1;
                }
                moduli[i-1] = moduli[i-1]*zeta*zeta;
            }
            if (cu.getUseDoublePrecision()) {
                if (dim == 0)
                    pmeBsplineModuliX.upload(moduli);
                else if (dim == 1)
                    pmeBsplineModuliY.upload(moduli);
                else
                    pmeBsplineModuliZ.upload(moduli);
            }
            else {
                vector<float> modulif(ndata);
                for (int i = 0; i < ndata; i++)
                    modulif[i] = (float) moduli[i];
                if (dim == 0)
                    pmeBsplineModuliX.upload(modulif);
                else if (dim == 1)
                    pmeBsplineModuliY.upload(modulif);
                else
                    pmeBsplineModuliZ.upload(modulif);
            }
        }

        // Initialize the b-spline moduli for dispersion PME.

        maxSize = max(max(dispersionGridSizeX, dispersionGridSizeY), dispersionGridSizeZ);
        vector<double> ddata(PmeOrder);
        bsplines_data.resize(maxSize);
        data[PmeOrder-1] = 0.0;
        data[1] = 0.0;
        data[0] = 1.0;
        for (int i = 3; i < PmeOrder; i++) {
            double div = 1.0/(i-1.0);
            data[i-1] = 0.0;
            for (int j = 1; j < (i-1); j++)
                data[i-j-1] = div*(j*data[i-j-2]+(i-j)*data[i-j-1]);
            data[0] = div*data[0];
        }

        // Differentiate.

        ddata[0] = -data[0];
        for (int i = 1; i < PmeOrder; i++)
            ddata[i] = data[i-1]-data[i];
        double div = 1.0/(PmeOrder-1);
        data[PmeOrder-1] = 0.0;
        for (int i = 1; i < (PmeOrder-1); i++)
            data[PmeOrder-i-1] = div*(i*data[PmeOrder-i-2]+(PmeOrder-i)*data[PmeOrder-i-1]);
        data[0] = div*data[0];
        for (int i = 0; i < maxSize; i++)
            bsplines_data[i] = 0.0;
        for (int i = 1; i <= PmeOrder; i++)
            bsplines_data[i] = data[i-1];

        // Evaluate the actual bspline moduli for X/Y/Z.

        for(int dim = 0; dim < 3; dim++) {
            int ndata = (dim == 0 ? dispersionGridSizeX : dim == 1 ? dispersionGridSizeY : dispersionGridSizeZ);
            vector<double> moduli(ndata);
            for (int i = 0; i < ndata; i++) {
                double sc = 0.0;
                double ss = 0.0;
                for (int j = 0; j < ndata; j++) {
                    double arg = (2.0*M_PI*i*j)/ndata;
                    sc += bsplines_data[j]*cos(arg);
                    ss += bsplines_data[j]*sin(arg);
                }
                moduli[i] = sc*sc+ss*ss;
            }
            for (int i = 0; i < ndata; i++)
                if (moduli[i] < 1.0e-7)
                    moduli[i] = (moduli[i-1]+moduli[i+1])*0.5;
            if (dim == 0)
                dpmeBsplineModuliX.upload(moduli, true);
            else if (dim == 1)
                dpmeBsplineModuliY.upload(moduli, true);
            else
                dpmeBsplineModuliZ.upload(moduli, true);
        }
    }

    // Add the interaction to the default nonbonded kernel.
    
    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
    nb.setKernelSource(CudaAmoebaKernelSources::hippoInteractionHeader+CudaAmoebaKernelSources::hippoNonbonded);
    nb.addArgument(CudaNonbondedUtilities::ParameterInfo("torqueBuffers", "unsigned long long", 1, torque.getElementSize(), torque.getDevicePointer(), false));
    nb.addArgument(CudaNonbondedUtilities::ParameterInfo("extrapolatedDipole", "real3", 1, extrapolatedDipole.getElementSize(), extrapolatedDipole.getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("coreCharge", "real", 1, coreCharge.getElementSize(), coreCharge.getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("valenceCharge", "real", 1, valenceCharge.getElementSize(), valenceCharge.getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("alpha", "real", 1, alpha.getElementSize(), alpha.getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("epsilon", "real", 1, epsilon.getElementSize(), epsilon.getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("damping", "real", 1, damping.getElementSize(), damping.getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("c6", "real", 1, c6.getElementSize(), c6.getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("pauliK", "real", 1, pauliK.getElementSize(), pauliK.getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("pauliQ", "real", 1, pauliQ.getElementSize(), pauliQ.getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("pauliAlpha", "real", 1, pauliAlpha.getElementSize(), pauliAlpha.getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("dipole", "real", 3, labDipoles.getElementSize(), labDipoles.getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("inducedDipole", "real", 3, inducedDipole.getElementSize(), inducedDipole.getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("qXX", "real", 1, labQuadrupoles[0].getElementSize(), labQuadrupoles[0].getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("qXY", "real", 1, labQuadrupoles[1].getElementSize(), labQuadrupoles[1].getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("qXZ", "real", 1, labQuadrupoles[2].getElementSize(), labQuadrupoles[2].getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("qYY", "real", 1, labQuadrupoles[3].getElementSize(), labQuadrupoles[3].getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("qYZ", "real", 1, labQuadrupoles[4].getElementSize(), labQuadrupoles[4].getDevicePointer()));
    map<string, string> replacements;
    replacements["ENERGY_SCALE_FACTOR"] = cu.doubleToString(ONE_4PI_EPS0);
    replacements["SWITCH_CUTOFF"] = cu.doubleToString(force.getSwitchingDistance());
    replacements["SWITCH_C3"] = cu.doubleToString(10/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 3.0));
    replacements["SWITCH_C4"] = cu.doubleToString(15/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 4.0));
    replacements["SWITCH_C5"] = cu.doubleToString(6/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 5.0));
    replacements["MAX_EXTRAPOLATION_ORDER"] = cu.intToString(maxExtrapolationOrder);
    replacements["EXTRAPOLATION_COEFFICIENTS_SUM"] = coefficients.str();
    replacements["USE_EWALD"] = (usePME ? "1" : "0");
    replacements["PME_ALPHA"] = (usePME ? cu.doubleToString(pmeAlpha) : "0");
    replacements["DPME_ALPHA"] = (usePME ? cu.doubleToString(dpmeAlpha) : "0");
    replacements["SQRT_PI"] = cu.doubleToString(sqrt(M_PI));
    string interactionSource = cu.replaceStrings(CudaAmoebaKernelSources::hippoInteraction, replacements);
    nb.addInteraction(usePME, usePME, true, force.getCutoffDistance(), exclusions, interactionSource, force.getForceGroup());
    nb.setUsePadding(false);
    
    // Create the kernel for computing exceptions.
    
    if (exceptionAtoms.isInitialized()) {
        replacements["COMPUTE_INTERACTION"] = interactionSource;
        string exceptionsSrc = CudaAmoebaKernelSources::hippoInteractionHeader+CudaAmoebaKernelSources::hippoNonbondedExceptions;
        exceptionsSrc = cu.replaceStrings(exceptionsSrc, replacements);
        defines["NUM_EXCEPTIONS"] = cu.intToString(exceptionAtoms.getSize());
        program = cu.compileProgram(exceptionsSrc, defines);
        computeExceptionsKernel = program->createKernel("computeNonbondedExceptions");
        computeExceptionsKernel->addArg(cu.getForce());
        computeExceptionsKernel->addArg(cu.getEnergyBuffer());
        computeExceptionsKernel->addArg(torque);
        computeExceptionsKernel->addArg(cu.getPosq());
        computeExceptionsKernel->addArg(extrapolatedDipole);
        computeExceptionsKernel->addArg(exceptionAtoms);
        computeExceptionsKernel->addArg(exceptionScales[0]);
        computeExceptionsKernel->addArg(exceptionScales[1]);
        computeExceptionsKernel->addArg(exceptionScales[2]);
        computeExceptionsKernel->addArg(exceptionScales[3]);
        computeExceptionsKernel->addArg(exceptionScales[4]);
        computeExceptionsKernel->addArg(exceptionScales[5]);
        computeExceptionsKernel->addArg(coreCharge);
        computeExceptionsKernel->addArg(valenceCharge);
        computeExceptionsKernel->addArg(alpha);
        computeExceptionsKernel->addArg(epsilon);
        computeExceptionsKernel->addArg(damping);
        computeExceptionsKernel->addArg(c6);
        computeExceptionsKernel->addArg(pauliK);
        computeExceptionsKernel->addArg(pauliQ);
        computeExceptionsKernel->addArg(pauliAlpha);
        computeExceptionsKernel->addArg(labDipoles);
        computeExceptionsKernel->addArg(inducedDipole);
        computeExceptionsKernel->addArg(labQuadrupoles[0]);
        computeExceptionsKernel->addArg(labQuadrupoles[1]);
        computeExceptionsKernel->addArg(labQuadrupoles[2]);
        computeExceptionsKernel->addArg(labQuadrupoles[3]);
        computeExceptionsKernel->addArg(labQuadrupoles[4]);
        computeExceptionsKernel->addArg(extrapolatedDipole);
        if (nb.getUseCutoff())
            for (int i = 0; i < 5; i++)
                computeExceptionsKernel->addArg();
    }
    cu.addForce(new ForceInfo(force));
    cu.addPostComputation(new TorquePostComputation(*this));
}

void CudaCalcHippoNonbondedForceKernel::createFieldKernel(const string& interactionSrc, vector<CudaArray*> params,
            CudaArray& fieldBuffer, ComputeKernel& kernel, ComputeKernel& exceptionKernel, CudaArray& exceptionScale) {
    // Create the kernel source.

    map<string, string> replacements;
    replacements["COMPUTE_FIELD"] = interactionSrc;
    stringstream extraArgs, atomParams, loadLocal1, loadLocal2, load1, load2, load3;
    for (auto param : params) {
        string name = param->getName();
        string type = (param->getElementSize() == 4 || param->getElementSize() == 8 ? "real" : "real3");
        extraArgs << ", const " << type << "* __restrict__ " << name;
        atomParams << type << " " << name << ";\n";
        loadLocal1 << "localData[localAtomIndex]." << name << " = " << name << "1;\n";
        loadLocal2 << "localData[localAtomIndex]." << name << " = " << name << "[j];\n";
        load1 << type << " " << name << "1 = " << name << "[atom1];\n";
        load2 << type << " " << name << "2 = localData[atom2]." << name << ";\n";
        load3 << type << " " << name << "2 = " << name << "[atom2];\n";
    }
    replacements["PARAMETER_ARGUMENTS"] = extraArgs.str();
    replacements["ATOM_PARAMETER_DATA"] = atomParams.str();
    replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
    replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
    replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
    replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
    replacements["LOAD_ATOM2_PARAMETERS_FROM_GLOBAL"] = load3.str();
    string src = cu.replaceStrings(CudaAmoebaKernelSources::hippoComputeField, replacements);

    // Set defines and create the kernel.

    map<string, string> defines;
    if (usePME) {
        defines["USE_CUTOFF"] = "1";
        defines["USE_PERIODIC"] = "1";
        defines["USE_EWALD"] = "1";
        defines["PME_ALPHA"] = cu.doubleToString(pmeAlpha);
        defines["SQRT_PI"] = cu.doubleToString(sqrt(M_PI));
    }
    defines["WARPS_PER_GROUP"] = cu.intToString(cu.getNonbondedUtilities().getForceThreadBlockSize()/CudaContext::TileSize);
    defines["THREAD_BLOCK_SIZE"] = cu.intToString(cu.getNonbondedUtilities().getForceThreadBlockSize());
    defines["CUTOFF"] = cu.doubleToString(cutoff);
    defines["CUTOFF_SQUARED"] = cu.doubleToString(cutoff*cutoff);
    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    defines["NUM_BLOCKS"] = cu.intToString(cu.getNumAtomBlocks());
    defines["TILE_SIZE"] = cu.intToString(CudaContext::TileSize);
    defines["NUM_TILES_WITH_EXCLUSIONS"] = cu.intToString(cu.getNonbondedUtilities().getExclusionTiles().getSize());
    defines["NUM_EXCEPTIONS"] = cu.intToString(exceptionAtoms.isInitialized() ? exceptionAtoms.getSize() : 0);
    ComputeProgram program = cu.compileProgram(src, defines);
    kernel = program->createKernel("computeField");

    // Build the list of arguments.

    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
    kernel->addArg(cu.getPosq());
    kernel->addArg(cu.getNonbondedUtilities().getExclusions());
    kernel->addArg(cu.getNonbondedUtilities().getExclusionTiles());
    kernel->addArg(fieldBuffer);
    if (nb.getUseCutoff()) {
        kernel->addArg(nb.getInteractingTiles());
        kernel->addArg(nb.getInteractionCount());
        for (int i = 0; i < 5; i++)
            kernel->addArg();
        kernel->addArg(maxTiles);
        kernel->addArg(nb.getBlockCenters());
        kernel->addArg(nb.getBlockBoundingBoxes());
        kernel->addArg(nb.getInteractingAtoms());
    }
    else
        kernel->addArg(maxTiles);
    for (auto param : params)
        kernel->addArg(*param);
    
    // If there are any exceptions, build the kernel and arguments to compute them.
    
    if (exceptionAtoms.isInitialized()) {
        exceptionKernel = program->createKernel("computeFieldExceptions");
        exceptionKernel->addArg(cu.getPosq());
        exceptionKernel->addArg(fieldBuffer);
        exceptionKernel->addArg(exceptionAtoms);
        exceptionKernel->addArg(exceptionScale);
        if (nb.getUseCutoff())
            for (int i = 0; i < 5; i++)
                exceptionKernel->addArg();
        for (auto param : params)
            exceptionKernel->addArg(*param);
    }
}

double CudaCalcHippoNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        
        // These kernels can't be compiled in initialize(), because the nonbonded utilities object
        // has not yet been initialized then.

        maxTiles = (nb.getUseCutoff() ? nb.getInteractingTiles().getSize() : cu.getNumAtomBlocks()*(cu.getNumAtomBlocks()+1)/2);
        createFieldKernel(CudaAmoebaKernelSources::hippoFixedField, {&coreCharge, &valenceCharge, &alpha, &labDipoles, &labQuadrupoles[0],
                &labQuadrupoles[1], &labQuadrupoles[2], &labQuadrupoles[3], &labQuadrupoles[4]}, field, fixedFieldKernel,
                fixedFieldExceptionKernel, exceptionScales[1]);
        createFieldKernel(CudaAmoebaKernelSources::hippoMutualField, {&alpha, &inducedDipole}, inducedField, mutualFieldKernel,
                mutualFieldExceptionKernel, exceptionScales[2]);
    }

    // Compute the lab frame moments.

    computeMomentsKernel->execute(cu.getNumAtoms());

    void* recipBoxVectorPointer[3];
    if (usePME) {
        setPeriodicBoxArgs(cu, dpmeGridIndexKernel, 2);
        setPeriodicBoxArgs(cu, dpmeSpreadChargeKernel, 2);
        setPeriodicBoxArgs(cu, dpmeInterpolateForceKernel, 3);
        
        // Compute reciprocal box vectors.
        
        Vec3 a, b, c;
        cu.getPeriodicBoxVectors(a, b, c);
        double determinant = a[0]*b[1]*c[2];
        double scale = 1.0/determinant;
        double4 recipBoxVectors[3];
        recipBoxVectors[0] = make_double4(b[1]*c[2]*scale, 0, 0, 0);
        recipBoxVectors[1] = make_double4(-b[0]*c[2]*scale, a[0]*c[2]*scale, 0, 0);
        recipBoxVectors[2] = make_double4((b[0]*c[1]-b[1]*c[0])*scale, -a[0]*c[1]*scale, a[0]*b[1]*scale, 0);
        if (cu.getUseDoublePrecision()) {
            mm_double4 boxVectors[] = {mm_double4(a[0], a[1], a[2], 0), mm_double4(b[0], b[1], b[2], 0), mm_double4(c[0], c[1], c[2], 0)};
            pmeConvolutionKernel->setArg(4, mm_double4(a[0], b[1], c[2], 0));
            dpmeConvolutionKernel->setArg(5, mm_double4(a[0], b[1], c[2], 0));
            dpmeEvalEnergyKernel->setArg(5, mm_double4(a[0], b[1], c[2], 0));
            for (int i = 0; i < 3; i++) {
                pmeTransformMultipolesKernel->setArg(8+i, recipBoxVectors[i]);
                pmeTransformPotentialKernel->setArg(2+i, recipBoxVectors[i]);
                pmeSpreadFixedMultipolesKernel->setArg(6+i, boxVectors[i]);
                pmeSpreadFixedMultipolesKernel->setArg(9+i, recipBoxVectors[i]);
                pmeSpreadInducedDipolesKernel->setArg(3+i, boxVectors[i]);
                pmeSpreadInducedDipolesKernel->setArg(6+i, recipBoxVectors[i]);
                pmeConvolutionKernel->setArg(5+i, recipBoxVectors[i]);
                pmeFixedPotentialKernel->setArg(5+i, boxVectors[i]);
                pmeFixedPotentialKernel->setArg(8+i, recipBoxVectors[i]);
                pmeInducedPotentialKernel->setArg(5+i, boxVectors[i]);
                pmeInducedPotentialKernel->setArg(8+i, recipBoxVectors[i]);
                pmeFixedForceKernel->setArg(16+i, recipBoxVectors[i]);
                pmeInducedForceKernel->setArg(20+i, recipBoxVectors[i]);
                pmeRecordInducedFieldDipolesKernel->setArg(3+i, recipBoxVectors[i]);
                dpmeGridIndexKernel->setArg(7+i, recipBoxVectors[i]);
                dpmeSpreadChargeKernel->setArg(7+i, recipBoxVectors[i]);
                dpmeConvolutionKernel->setArg(6+i, recipBoxVectors[i]);
                dpmeEvalEnergyKernel->setArg(6+i, recipBoxVectors[i]);
                dpmeInterpolateForceKernel->setArg(8+i, recipBoxVectors[i]);
            }
        }
        else {
            mm_float4 boxVectors[] = {mm_float4(a[0], a[1], a[2], 0), mm_float4(b[0], b[1], b[2], 0), mm_float4(c[0], c[1], c[2], 0)};
            pmeConvolutionKernel->setArg(4, mm_float4(a[0], b[1], c[2], 0));
            dpmeConvolutionKernel->setArg(5, mm_float4(a[0], b[1], c[2], 0));
            dpmeEvalEnergyKernel->setArg(5, mm_float4(a[0], b[1], c[2], 0));
            float4 recipBoxVectorsFloat[3];
            recipBoxVectorsFloat[0] = make_float4((float) recipBoxVectors[0].x, 0, 0, 0);
            recipBoxVectorsFloat[1] = make_float4((float) recipBoxVectors[1].x, (float) recipBoxVectors[1].y, 0, 0);
            recipBoxVectorsFloat[2] = make_float4((float) recipBoxVectors[2].x, (float) recipBoxVectors[2].y, (float) recipBoxVectors[2].z, 0);
            for (int i = 0; i < 3; i++) {
                pmeTransformMultipolesKernel->setArg(8+i, recipBoxVectorsFloat[i]);
                pmeTransformPotentialKernel->setArg(2+i, recipBoxVectorsFloat[i]);
                pmeSpreadFixedMultipolesKernel->setArg(6+i, boxVectors[i]);
                pmeSpreadFixedMultipolesKernel->setArg(9+i, recipBoxVectorsFloat[i]);
                pmeSpreadInducedDipolesKernel->setArg(3+i, boxVectors[i]);
                pmeSpreadInducedDipolesKernel->setArg(6+i, recipBoxVectorsFloat[i]);
                pmeConvolutionKernel->setArg(5+i, recipBoxVectorsFloat[i]);
                pmeFixedPotentialKernel->setArg(5+i, boxVectors[i]);
                pmeFixedPotentialKernel->setArg(8+i, recipBoxVectorsFloat[i]);
                pmeInducedPotentialKernel->setArg(5+i, boxVectors[i]);
                pmeInducedPotentialKernel->setArg(8+i, recipBoxVectorsFloat[i]);
                pmeFixedForceKernel->setArg(16+i, recipBoxVectorsFloat[i]);
                pmeInducedForceKernel->setArg(20+i, recipBoxVectorsFloat[i]);
                pmeRecordInducedFieldDipolesKernel->setArg(3+i, recipBoxVectorsFloat[i]);
                dpmeGridIndexKernel->setArg(7+i, recipBoxVectorsFloat[i]);
                dpmeSpreadChargeKernel->setArg(7+i, recipBoxVectorsFloat[i]);
                dpmeConvolutionKernel->setArg(6+i, recipBoxVectorsFloat[i]);
                dpmeEvalEnergyKernel->setArg(6+i, recipBoxVectorsFloat[i]);
                dpmeInterpolateForceKernel->setArg(8+i, recipBoxVectorsFloat[i]);
            }
        }



        float3 recipBoxVectorsFloat[3];
        if (cu.getUseDoublePrecision()) {
            recipBoxVectorPointer[0] = &recipBoxVectors[0];
            recipBoxVectorPointer[1] = &recipBoxVectors[1];
            recipBoxVectorPointer[2] = &recipBoxVectors[2];
        }
        else {
            recipBoxVectorsFloat[0] = make_float3((float) recipBoxVectors[0].x, 0, 0);
            recipBoxVectorsFloat[1] = make_float3((float) recipBoxVectors[1].x, (float) recipBoxVectors[1].y, 0);
            recipBoxVectorsFloat[2] = make_float3((float) recipBoxVectors[2].x, (float) recipBoxVectors[2].y, (float) recipBoxVectors[2].z);
            recipBoxVectorPointer[0] = &recipBoxVectorsFloat[0];
            recipBoxVectorPointer[1] = &recipBoxVectorsFloat[1];
            recipBoxVectorPointer[2] = &recipBoxVectorsFloat[2];
        }

        // Reciprocal space calculation for electrostatics.
        
        pmeTransformMultipolesKernel->execute(cu.getNumAtoms());
        pmeSpreadFixedMultipolesKernel->execute(cu.getNumAtoms());
        if (cu.getUseDoublePrecision()) {
            pmeFinishSpreadChargeKernel->execute(pmeGrid1.getSize());
            cufftExecD2Z(fftForward, (double*) pmeGrid1.getDevicePointer(), (double2*) pmeGrid2.getDevicePointer());
        }
        else
            cufftExecR2C(fftForward, (float*) pmeGrid1.getDevicePointer(), (float2*) pmeGrid2.getDevicePointer());
        pmeConvolutionKernel->execute(gridSizeX*gridSizeY*gridSizeZ, 256);
        if (cu.getUseDoublePrecision())
            cufftExecZ2D(fftBackward, (double2*) pmeGrid2.getDevicePointer(), (double*) pmeGrid1.getDevicePointer());
        else
            cufftExecC2R(fftBackward, (float2*) pmeGrid2.getDevicePointer(), (float*) pmeGrid1.getDevicePointer());
        pmeFixedPotentialKernel->execute(cu.getNumAtoms());
        pmeTransformPotentialKernel->setArg(0, pmePhi);
        pmeTransformPotentialKernel->execute(cu.getNumAtoms());
        pmeFixedForceKernel->execute(cu.getNumAtoms());

        // Reciprocal space calculation for dispersion.

        dpmeGridIndexKernel->execute(cu.getNumAtoms());
        sort->sort(pmeAtomGridIndex);
        cu.clearBuffer(pmeGrid2);
        dpmeSpreadChargeKernel->execute(cu.getNumAtoms(), 128);
        dpmeFinishSpreadChargeKernel->execute(dispersionGridSizeX*dispersionGridSizeY*dispersionGridSizeZ, 256);
        if (cu.getUseDoublePrecision())
            cufftExecD2Z(dfftForward, (double*) pmeGrid1.getDevicePointer(), (double2*) pmeGrid2.getDevicePointer());
        else
            cufftExecR2C(dfftForward, (float*) pmeGrid1.getDevicePointer(), (float2*) pmeGrid2.getDevicePointer());
        if (includeEnergy)
            dpmeEvalEnergyKernel->execute(dispersionGridSizeX*dispersionGridSizeY*dispersionGridSizeZ);
        dpmeConvolutionKernel->execute(dispersionGridSizeX*dispersionGridSizeY*dispersionGridSizeZ, 256);
        if (cu.getUseDoublePrecision())
            cufftExecZ2D(dfftBackward, (double2*) pmeGrid2.getDevicePointer(), (double*) pmeGrid1.getDevicePointer());
        else
            cufftExecC2R(dfftBackward, (float2*) pmeGrid2.getDevicePointer(), (float*)  pmeGrid1.getDevicePointer());
        dpmeInterpolateForceKernel->execute(cu.getNumAtoms(), 128);
    }

    // Compute the field from fixed multipoles.

    if (nb.getUseCutoff())
        setPeriodicBoxArgs(cu, fixedFieldKernel, 6);
    fixedFieldKernel->execute(nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    if (exceptionAtoms.isInitialized()) {
        if (nb.getUseCutoff())
            setPeriodicBoxArgs(cu, fixedFieldExceptionKernel, 4);
        fixedFieldExceptionKernel->execute(exceptionAtoms.getSize());
    }

    // Iterate the induced dipoles.

    computeExtrapolatedDipoles(recipBoxVectorPointer);

    // Add the polarization energy.

    if (includeEnergy)
        polarizationEnergyKernel->execute(cu.getNumAtoms());

    // Compute the forces due to the reciprocal space PME calculation for induced dipoles.

    if (usePME) {
        pmeTransformPotentialKernel->setArg(0, pmePhidp);
        pmeTransformPotentialKernel->execute(cu.getNumAtoms());
        pmeInducedForceKernel->execute(cu.getNumAtoms());
        pmeSelfEnergyKernel->execute(cu.getNumAtoms());
    }

    // Compute nonbonded exceptions.

    if (exceptionAtoms.isInitialized()) {
        if (nb.getUseCutoff())
            setPeriodicBoxArgs(cu, computeExceptionsKernel, 29);
        computeExceptionsKernel->execute(exceptionAtoms.getSize());
    }

    // Record the current atom positions so we can tell later if they have changed.
    
    cu.getPosq().copyTo(lastPositions);
    multipolesAreValid = true;
    return 0.0;
}

void CudaCalcHippoNonbondedForceKernel::computeInducedField(void** recipBoxVectorPointer, int optOrder) {
    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
    cu.clearBuffer(inducedField);
    if (nb.getUseCutoff())
        setPeriodicBoxArgs(cu, mutualFieldKernel, 6);
    mutualFieldKernel->execute(nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    if (exceptionAtoms.isInitialized()) {
        if (nb.getUseCutoff())
            setPeriodicBoxArgs(cu, mutualFieldExceptionKernel, 4);
        mutualFieldExceptionKernel->execute(exceptionAtoms.getSize());
    }
    if (usePME) {
        cu.clearBuffer(pmeGrid1);
        pmeSpreadInducedDipolesKernel->execute(cu.getNumAtoms());
        if (cu.getUseDoublePrecision()) {
            pmeFinishSpreadChargeKernel->execute(pmeGrid1.getSize());
            cufftExecD2Z(fftForward, (double*) pmeGrid1.getDevicePointer(), (double2*) pmeGrid2.getDevicePointer());
        }
        else
            cufftExecR2C(fftForward, (float*) pmeGrid1.getDevicePointer(), (float2*) pmeGrid2.getDevicePointer());
        pmeConvolutionKernel->execute(gridSizeX*gridSizeY*gridSizeZ, 256);
        if (cu.getUseDoublePrecision())
            cufftExecZ2D(fftBackward, (double2*) pmeGrid2.getDevicePointer(), (double*) pmeGrid1.getDevicePointer());
        else
            cufftExecC2R(fftBackward, (float2*) pmeGrid2.getDevicePointer(), (float*) pmeGrid1.getDevicePointer());
        pmeInducedPotentialKernel->setArg(2, optOrder);
        pmeInducedPotentialKernel->execute(cu.getNumAtoms());
        pmeRecordInducedFieldDipolesKernel->execute(cu.getNumAtoms());
    }
}

void CudaCalcHippoNonbondedForceKernel::computeExtrapolatedDipoles(void** recipBoxVectorPointer) {
    // Start by storing the direct dipoles as PT0

    recordInducedDipolesKernel->execute(cu.getNumAtoms());
    initExtrapolatedKernel->execute(extrapolatedDipole.getSize());

    // Recursively apply alpha.Tau to the µ_(n) components to generate µ_(n+1), and store the result

    for (int order = 1; order < maxExtrapolationOrder; ++order) {
        computeInducedField(recipBoxVectorPointer, order-1);
        iterateExtrapolatedKernel->setArg(0, order);
        iterateExtrapolatedKernel->execute(extrapolatedDipole.getSize());
    }
    
    // Take a linear combination of the µ_(n) components to form the total dipole

    computeExtrapolatedKernel->execute(extrapolatedDipole.getSize());
    computeInducedField(recipBoxVectorPointer, maxExtrapolationOrder-1);
}

void CudaCalcHippoNonbondedForceKernel::addTorquesToForces() {
    mapTorqueKernel->execute(cu.getNumAtoms());
}

void CudaCalcHippoNonbondedForceKernel::getInducedDipoles(ContextImpl& context, vector<Vec3>& dipoles) {
    ensureMultipolesValid(context);
    int numParticles = cu.getNumAtoms();
    dipoles.resize(numParticles);
    const vector<int>& order = cu.getAtomIndex();
    if (cu.getUseDoublePrecision()) {
        vector<double3> d;
        inducedDipole.download(d);
        for (int i = 0; i < numParticles; i++)
            dipoles[order[i]] = Vec3(d[i].x, d[i].y, d[i].z);
    }
    else {
        vector<float3> d;
        inducedDipole.download(d);
        for (int i = 0; i < numParticles; i++)
            dipoles[order[i]] = Vec3(d[i].x, d[i].y, d[i].z);
    }
}

void CudaCalcHippoNonbondedForceKernel::ensureMultipolesValid(ContextImpl& context) {
    if (multipolesAreValid) {
        int numParticles = cu.getNumAtoms();
        if (cu.getUseDoublePrecision()) {
            vector<double4> pos1, pos2;
            cu.getPosq().download(pos1);
            lastPositions.download(pos2);
            for (int i = 0; i < numParticles; i++)
                if (pos1[i].x != pos2[i].x || pos1[i].y != pos2[i].y || pos1[i].z != pos2[i].z) {
                    multipolesAreValid = false;
                    break;
                }
        }
        else {
            vector<float4> pos1, pos2;
            cu.getPosq().download(pos1);
            lastPositions.download(pos2);
            for (int i = 0; i < numParticles; i++)
                if (pos1[i].x != pos2[i].x || pos1[i].y != pos2[i].y || pos1[i].z != pos2[i].z) {
                    multipolesAreValid = false;
                    break;
                }
        }
    }
    if (!multipolesAreValid)
        context.calcForcesAndEnergy(false, false, context.getIntegrator().getIntegrationForceGroups());
}

void CudaCalcHippoNonbondedForceKernel::getLabFramePermanentDipoles(ContextImpl& context, vector<Vec3>& dipoles) {
    ensureMultipolesValid(context);
    int numParticles = cu.getNumAtoms();
    dipoles.resize(numParticles);
    const vector<int>& order = cu.getAtomIndex();
    if (cu.getUseDoublePrecision()) {
        vector<double3> labDipoleVec;
        labDipoles.download(labDipoleVec);
        for (int i = 0; i < numParticles; i++)
            dipoles[order[i]] = Vec3(labDipoleVec[i].x, labDipoleVec[i].y, labDipoleVec[i].z);
    }
    else {
        vector<float3> labDipoleVec;
        labDipoles.download(labDipoleVec);
        for (int i = 0; i < numParticles; i++)
            dipoles[order[i]] = Vec3(labDipoleVec[i].x, labDipoleVec[i].y, labDipoleVec[i].z);
    }
}

void CudaCalcHippoNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const HippoNonbondedForce& force) {
    // Make sure the new parameters are acceptable.
    
    cu.setAsCurrent();
    if (force.getNumParticles() != cu.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    
    // Record the per-particle parameters.
    
    vector<double> coreChargeVec, valenceChargeVec, alphaVec, epsilonVec, dampingVec, c6Vec, pauliKVec, pauliQVec, pauliAlphaVec, polarizabilityVec;
    vector<double> localDipolesVec, localQuadrupolesVec;
    vector<int4> multipoleParticlesVec;
    for (int i = 0; i < numParticles; i++) {
        double charge, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability;
        int axisType, atomX, atomY, atomZ;
        vector<double> dipole, quadrupole;
        force.getParticleParameters(i, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha,
                                    polarizability, axisType, atomZ, atomX, atomY);
        coreChargeVec.push_back(coreCharge);
        valenceChargeVec.push_back(charge-coreCharge);
        alphaVec.push_back(alpha);
        epsilonVec.push_back(epsilon);
        dampingVec.push_back(damping);
        c6Vec.push_back(c6);
        pauliKVec.push_back(pauliK);
        pauliQVec.push_back(pauliQ);
        pauliAlphaVec.push_back(pauliAlpha);
        polarizabilityVec.push_back(polarizability);
        multipoleParticlesVec.push_back(make_int4(atomX, atomY, atomZ, axisType));
        for (int j = 0; j < 3; j++)
            localDipolesVec.push_back(dipole[j]);
        localQuadrupolesVec.push_back(quadrupole[0]);
        localQuadrupolesVec.push_back(quadrupole[1]);
        localQuadrupolesVec.push_back(quadrupole[2]);
        localQuadrupolesVec.push_back(quadrupole[4]);
        localQuadrupolesVec.push_back(quadrupole[5]);
    }
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    for (int i = numParticles; i < paddedNumAtoms; i++) {
        coreChargeVec.push_back(0);
        valenceChargeVec.push_back(0);
        alphaVec.push_back(0);
        epsilonVec.push_back(0);
        dampingVec.push_back(0);
        c6Vec.push_back(0);
        pauliKVec.push_back(0);
        pauliQVec.push_back(0);
        pauliAlphaVec.push_back(0);
        polarizabilityVec.push_back(0);
        multipoleParticlesVec.push_back(make_int4(0, 0, 0, 0));
        for (int j = 0; j < 3; j++)
            localDipolesVec.push_back(0);
        for (int j = 0; j < 5; j++)
            localQuadrupolesVec.push_back(0);
    }
    coreCharge.upload(coreChargeVec, true);
    valenceCharge.upload(valenceChargeVec, true);
    alpha.upload(alphaVec, true);
    epsilon.upload(epsilonVec, true);
    damping.upload(dampingVec, true);
    c6.upload(c6Vec, true);
    pauliK.upload(pauliKVec, true);
    pauliQ.upload(pauliQVec, true);
    pauliAlpha.upload(pauliAlphaVec, true);
    polarizability.upload(polarizabilityVec, true);
    multipoleParticles.upload(multipoleParticlesVec);
    localDipoles.upload(localDipolesVec, true);
    localQuadrupoles.upload(localQuadrupolesVec, true);
    
    // Record the per-exception parameters.

    vector<double> exceptionScaleVec[6];
    vector<int2> exceptionAtomsVec;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale;
        force.getExceptionParameters(i, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
        if (usePME || multipoleMultipoleScale != 0 || dipoleMultipoleScale != 0 || dipoleDipoleScale != 0 || dispersionScale != 0 || repulsionScale != 0 || chargeTransferScale != 0) {
            exceptionAtomsVec.push_back(make_int2(particle1, particle2));
            exceptionScaleVec[0].push_back(multipoleMultipoleScale);
            exceptionScaleVec[1].push_back(dipoleMultipoleScale);
            exceptionScaleVec[2].push_back(dipoleDipoleScale);
            exceptionScaleVec[3].push_back(dispersionScale);
            exceptionScaleVec[4].push_back(repulsionScale);
            exceptionScaleVec[5].push_back(chargeTransferScale);
        }
    }
    if (exceptionAtomsVec.size() > 0) {
        if (!exceptionAtoms.isInitialized() || exceptionAtoms.getSize() != exceptionAtomsVec.size())
            throw OpenMMException("updateParametersInContext: The number of exceptions has changed");
        exceptionAtoms.upload(exceptionAtomsVec);
        for (int i = 0; i < 6; i++)
            exceptionScales[i].upload(exceptionScaleVec[i], true);
    }
    else if (exceptionAtoms.isInitialized())
        throw OpenMMException("updateParametersInContext: The number of exceptions has changed");
    cu.invalidateMolecules();
    multipolesAreValid = false;
}

void CudaCalcHippoNonbondedForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    alpha = pmeAlpha;
    nx = gridSizeX;
    ny = gridSizeY;
    nz = gridSizeZ;
}

void CudaCalcHippoNonbondedForceKernel::getDPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    alpha = dpmeAlpha;
    nx = dispersionGridSizeX;
    ny = dispersionGridSizeY;
    nz = dispersionGridSizeZ;
}

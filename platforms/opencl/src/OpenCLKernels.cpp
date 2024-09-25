/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

#include "OpenCLKernels.h"
#include "OpenCLForceInfo.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "CommonKernelSources.h"
#include "OpenCLBondedUtilities.h"
#include "OpenCLExpressionUtilities.h"
#include "OpenCLIntegrationUtilities.h"
#include "OpenCLNonbondedUtilities.h"
#include "OpenCLKernelSources.h"
#include "SimTKOpenMMRealType.h"
#include "SimTKOpenMMUtilities.h"
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iterator>
#include <set>

using namespace OpenMM;
using namespace std;

static void setPeriodicBoxSizeArg(OpenCLContext& cl, cl::Kernel& kernel, int index) {
    if (cl.getUseDoublePrecision())
        kernel.setArg<mm_double4>(index, cl.getPeriodicBoxSizeDouble());
    else
        kernel.setArg<mm_float4>(index, cl.getPeriodicBoxSize());
}

static void setPeriodicBoxArgs(OpenCLContext& cl, cl::Kernel& kernel, int index) {
    if (cl.getUseDoublePrecision()) {
        kernel.setArg<mm_double4>(index++, cl.getPeriodicBoxSizeDouble());
        kernel.setArg<mm_double4>(index++, cl.getInvPeriodicBoxSizeDouble());
        kernel.setArg<mm_double4>(index++, cl.getPeriodicBoxVecXDouble());
        kernel.setArg<mm_double4>(index++, cl.getPeriodicBoxVecYDouble());
        kernel.setArg<mm_double4>(index, cl.getPeriodicBoxVecZDouble());
    }
    else {
        kernel.setArg<mm_float4>(index++, cl.getPeriodicBoxSize());
        kernel.setArg<mm_float4>(index++, cl.getInvPeriodicBoxSize());
        kernel.setArg<mm_float4>(index++, cl.getPeriodicBoxVecX());
        kernel.setArg<mm_float4>(index++, cl.getPeriodicBoxVecY());
        kernel.setArg<mm_float4>(index, cl.getPeriodicBoxVecZ());
    }
}

void OpenCLCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void OpenCLCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    cl.setForcesValid(true);
    cl.clearAutoclearBuffers();
    for (auto computation : cl.getPreComputations())
        computation->computeForceAndEnergy(includeForces, includeEnergy, groups);
    OpenCLNonbondedUtilities& nb = cl.getNonbondedUtilities();
    cl.setComputeForceCount(cl.getComputeForceCount()+1);
    nb.prepareInteractions(groups);
    map<string, double>& derivs = cl.getEnergyParamDerivWorkspace();
    for (auto& param : context.getParameters())
        derivs[param.first] = 0;
}

double OpenCLCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups, bool& valid) {
    cl.getBondedUtilities().computeInteractions(groups);
    cl.getNonbondedUtilities().computeInteractions(groups, includeForces, includeEnergy);
    double sum = 0.0;
    for (auto computation : cl.getPostComputations())
        sum += computation->computeForceAndEnergy(includeForces, includeEnergy, groups);
    cl.reduceForces();
    cl.getIntegrationUtilities().distributeForcesFromVirtualSites();
    if (includeEnergy)
        sum += cl.reduceEnergy();
    if (!cl.getForcesValid())
        valid = false;
    return sum;
}

class OpenCLCalcNonbondedForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(int requiredBuffers, const NonbondedForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, charge2, sigma1, sigma2, epsilon1, epsilon2;
        force.getParticleParameters(particle1, charge1, sigma1, epsilon1);
        force.getParticleParameters(particle2, charge2, sigma2, epsilon2);
        return (charge1 == charge2 && sigma1 == sigma2 && epsilon1 == epsilon2);
    }
    int getNumParticleGroups() {
        return force.getNumExceptions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(index, particle1, particle2, chargeProd, sigma, epsilon);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        double chargeProd1, chargeProd2, sigma1, sigma2, epsilon1, epsilon2;
        force.getExceptionParameters(group1, particle1, particle2, chargeProd1, sigma1, epsilon1);
        force.getExceptionParameters(group2, particle1, particle2, chargeProd2, sigma2, epsilon2);
        return (chargeProd1 == chargeProd2 && sigma1 == sigma2 && epsilon1 == epsilon2);
    }
private:
    const NonbondedForce& force;
};

class OpenCLCalcNonbondedForceKernel::PmeIO : public CalcPmeReciprocalForceKernel::IO {
public:
    PmeIO(OpenCLContext& cl, cl::Kernel addForcesKernel) : cl(cl), addForcesKernel(addForcesKernel) {
        forceTemp.initialize<mm_float4>(cl, cl.getNumAtoms(), "PmeForce");
        addForcesKernel.setArg<cl::Buffer>(0, forceTemp.getDeviceBuffer());
    }
    float* getPosq() {
        cl.getPosq().download(posq);
        return (float*) &posq[0];
    }
    void setForce(float* force) {
        forceTemp.upload(force);
        addForcesKernel.setArg<cl::Buffer>(1, cl.getLongForceBuffer().getDeviceBuffer());
        cl.executeKernel(addForcesKernel, cl.getNumAtoms());
    }
private:
    OpenCLContext& cl;
    vector<mm_float4> posq;
    OpenCLArray forceTemp;
    cl::Kernel addForcesKernel;
};

class OpenCLCalcNonbondedForceKernel::PmePreComputation : public OpenCLContext::ForcePreComputation {
public:
    PmePreComputation(OpenCLContext& cl, Kernel& pme, CalcPmeReciprocalForceKernel::IO& io) : cl(cl), pme(pme), io(io) {
    }
    void computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        Vec3 boxVectors[3] = {Vec3(cl.getPeriodicBoxSize().x, 0, 0), Vec3(0, cl.getPeriodicBoxSize().y, 0), Vec3(0, 0, cl.getPeriodicBoxSize().z)};
        pme.getAs<CalcPmeReciprocalForceKernel>().beginComputation(io, boxVectors, includeEnergy);
    }
private:
    OpenCLContext& cl;
    Kernel pme;
    CalcPmeReciprocalForceKernel::IO& io;
};

class OpenCLCalcNonbondedForceKernel::PmePostComputation : public OpenCLContext::ForcePostComputation {
public:
    PmePostComputation(Kernel& pme, CalcPmeReciprocalForceKernel::IO& io) : pme(pme), io(io) {
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        return pme.getAs<CalcPmeReciprocalForceKernel>().finishComputation(io);
    }
private:
    Kernel pme;
    CalcPmeReciprocalForceKernel::IO& io;
};

class OpenCLCalcNonbondedForceKernel::SyncQueuePreComputation : public OpenCLContext::ForcePreComputation {
public:
    SyncQueuePreComputation(OpenCLContext& cl, cl::CommandQueue queue, int forceGroup) : cl(cl), queue(queue), forceGroup(forceGroup) {
    }
    void computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        if ((groups&(1<<forceGroup)) != 0) {
            vector<cl::Event> events(1);
            cl.getQueue().enqueueMarkerWithWaitList(NULL, &events[0]);
            queue.enqueueBarrierWithWaitList(&events);
        }
    }
private:
    OpenCLContext& cl;
    cl::CommandQueue queue;
    int forceGroup;
};

class OpenCLCalcNonbondedForceKernel::SyncQueuePostComputation : public OpenCLContext::ForcePostComputation {
public:
    SyncQueuePostComputation(OpenCLContext& cl, cl::Event& event, OpenCLArray& pmeEnergyBuffer, int forceGroup) : cl(cl), event(event),
            pmeEnergyBuffer(pmeEnergyBuffer), forceGroup(forceGroup) {
    }
    void setKernel(cl::Kernel kernel) {
        addEnergyKernel = kernel;
        addEnergyKernel.setArg<cl::Buffer>(0, pmeEnergyBuffer.getDeviceBuffer());
        addEnergyKernel.setArg<cl::Buffer>(1, cl.getEnergyBuffer().getDeviceBuffer());
        addEnergyKernel.setArg<cl_int>(2, pmeEnergyBuffer.getSize());
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        if ((groups&(1<<forceGroup)) != 0) {
            vector<cl::Event> events(1);
            events[0] = event;
            event = cl::Event();
            cl.getQueue().enqueueBarrierWithWaitList(&events);
            if (includeEnergy)
                cl.executeKernel(addEnergyKernel, pmeEnergyBuffer.getSize());
        }
        return 0.0;
    }
private:
    OpenCLContext& cl;
    cl::Event& event;
    cl::Kernel addEnergyKernel;
    OpenCLArray& pmeEnergyBuffer;
    int forceGroup;
};

OpenCLCalcNonbondedForceKernel::~OpenCLCalcNonbondedForceKernel() {
    if (sort != NULL)
        delete sort;
    if (fft != NULL)
        delete fft;
    if (dispersionFft != NULL)
        delete dispersionFft;
    if (pmeio != NULL)
        delete pmeio;
}

void OpenCLCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
        ;
    string prefix = "nonbonded"+cl.intToString(forceIndex)+"_";

    // Identify which exceptions are 1-4 interactions.

    set<int> exceptionsWithOffsets;
    for (int i = 0; i < force.getNumExceptionParameterOffsets(); i++) {
        string param;
        int exception;
        double charge, sigma, epsilon;
        force.getExceptionParameterOffset(i, param, exception, charge, sigma, epsilon);
        exceptionsWithOffsets.insert(exception);
    }
    vector<pair<int, int> > exclusions;
    vector<int> exceptions;
    map<int, int> exceptionIndex;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        exclusions.push_back(pair<int, int>(particle1, particle2));
        if (chargeProd != 0.0 || epsilon != 0.0 || exceptionsWithOffsets.find(i) != exceptionsWithOffsets.end()) {
            exceptionIndex[i] = exceptions.size();
            exceptions.push_back(i);
        }
    }

    // Initialize nonbonded interactions.

    int numParticles = force.getNumParticles();
    vector<mm_float4> baseParticleParamVec(cl.getPaddedNumAtoms(), mm_float4(0, 0, 0, 0));
    vector<vector<int> > exclusionList(numParticles);
    hasCoulomb = false;
    hasLJ = false;
    for (int i = 0; i < numParticles; i++) {
        double charge, sigma, epsilon;
        force.getParticleParameters(i, charge, sigma, epsilon);
        baseParticleParamVec[i] = mm_float4(charge, sigma, epsilon, 0);
        exclusionList[i].push_back(i);
        if (charge != 0.0)
            hasCoulomb = true;
        if (epsilon != 0.0)
            hasLJ = true;
    }
    for (int i = 0; i < force.getNumParticleParameterOffsets(); i++) {
        string param;
        int particle;
        double charge, sigma, epsilon;
        force.getParticleParameterOffset(i, param, particle, charge, sigma, epsilon);
        if (charge != 0.0)
            hasCoulomb = true;
        if (epsilon != 0.0)
            hasLJ = true;
    }
    for (auto exclusion : exclusions) {
        exclusionList[exclusion.first].push_back(exclusion.second);
        exclusionList[exclusion.second].push_back(exclusion.first);
    }
    nonbondedMethod = CalcNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
    bool useCutoff = (nonbondedMethod != NoCutoff);
    bool usePeriodic = (nonbondedMethod != NoCutoff && nonbondedMethod != CutoffNonPeriodic);
    doLJPME = (nonbondedMethod == LJPME && hasLJ);
    usePosqCharges = hasCoulomb ? cl.requestPosqCharges() : false;
    map<string, string> defines;
    defines["HAS_COULOMB"] = (hasCoulomb ? "1" : "0");
    defines["HAS_LENNARD_JONES"] = (hasLJ ? "1" : "0");
    defines["USE_LJ_SWITCH"] = (useCutoff && force.getUseSwitchingFunction() ? "1" : "0");
    if (useCutoff) {
        // Compute the reaction field constants.

        double reactionFieldK = pow(force.getCutoffDistance(), -3.0)*(force.getReactionFieldDielectric()-1.0)/(2.0*force.getReactionFieldDielectric()+1.0);
        double reactionFieldC = (1.0 / force.getCutoffDistance())*(3.0*force.getReactionFieldDielectric())/(2.0*force.getReactionFieldDielectric()+1.0);
        defines["REACTION_FIELD_K"] = cl.doubleToString(reactionFieldK);
        defines["REACTION_FIELD_C"] = cl.doubleToString(reactionFieldC);
        
        // Compute the switching coefficients.
        
        if (force.getUseSwitchingFunction()) {
            defines["LJ_SWITCH_CUTOFF"] = cl.doubleToString(force.getSwitchingDistance());
            defines["LJ_SWITCH_C3"] = cl.doubleToString(10/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 3.0));
            defines["LJ_SWITCH_C4"] = cl.doubleToString(15/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 4.0));
            defines["LJ_SWITCH_C5"] = cl.doubleToString(6/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 5.0));
        }
    }
    if (force.getUseDispersionCorrection() && cl.getContextIndex() == 0 && !doLJPME)
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(system, force);
    else
        dispersionCoefficient = 0.0;
    alpha = 0;
    ewaldSelfEnergy = 0.0;
    map<string, string> paramsDefines;
    paramsDefines["ONE_4PI_EPS0"] = cl.doubleToString(ONE_4PI_EPS0);
    hasOffsets = (force.getNumParticleParameterOffsets() > 0 || force.getNumExceptionParameterOffsets() > 0);
    if (hasOffsets)
        paramsDefines["HAS_OFFSETS"] = "1";
    if (force.getNumParticleParameterOffsets() > 0)
        paramsDefines["HAS_PARTICLE_OFFSETS"] = "1";
    if (force.getNumExceptionParameterOffsets() > 0)
        paramsDefines["HAS_EXCEPTION_OFFSETS"] = "1";
    if (usePosqCharges)
        paramsDefines["USE_POSQ_CHARGES"] = "1";
    if (doLJPME)
        paramsDefines["INCLUDE_LJPME_EXCEPTIONS"] = "1";
    if (nonbondedMethod == Ewald) {
        // Compute the Ewald parameters.

        int kmaxx, kmaxy, kmaxz;
        NonbondedForceImpl::calcEwaldParameters(system, force, alpha, kmaxx, kmaxy, kmaxz);
        defines["EWALD_ALPHA"] = cl.doubleToString(alpha);
        defines["TWO_OVER_SQRT_PI"] = cl.doubleToString(2.0/sqrt(M_PI));
        defines["USE_EWALD"] = "1";
        if (cl.getContextIndex() == 0) {
            paramsDefines["INCLUDE_EWALD"] = "1";
            paramsDefines["EWALD_SELF_ENERGY_SCALE"] = cl.doubleToString(ONE_4PI_EPS0*alpha/sqrt(M_PI));
            for (int i = 0; i < numParticles; i++)
                ewaldSelfEnergy -= baseParticleParamVec[i].x*baseParticleParamVec[i].x*ONE_4PI_EPS0*alpha/sqrt(M_PI);

            // Create the reciprocal space kernels.

            map<string, string> replacements;
            replacements["NUM_ATOMS"] = cl.intToString(numParticles);
            replacements["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
            replacements["KMAX_X"] = cl.intToString(kmaxx);
            replacements["KMAX_Y"] = cl.intToString(kmaxy);
            replacements["KMAX_Z"] = cl.intToString(kmaxz);
            replacements["EXP_COEFFICIENT"] = cl.doubleToString(-1.0/(4.0*alpha*alpha));
            replacements["ONE_4PI_EPS0"] = cl.doubleToString(ONE_4PI_EPS0);
            replacements["M_PI"] = cl.doubleToString(M_PI);
            cl::Program program = cl.createProgram(CommonKernelSources::ewald, replacements);
            ewaldSumsKernel = cl::Kernel(program, "calculateEwaldCosSinSums");
            ewaldForcesKernel = cl::Kernel(program, "calculateEwaldForces");
            int elementSize = (cl.getUseDoublePrecision() ? sizeof(mm_double2) : sizeof(mm_float2));
            cosSinSums.initialize(cl, (2*kmaxx-1)*(2*kmaxy-1)*(2*kmaxz-1), elementSize, "cosSinSums");
        }
    }
    else if (((nonbondedMethod == PME || nonbondedMethod == LJPME) && hasCoulomb) || doLJPME) {
        // Compute the PME parameters.

        NonbondedForceImpl::calcPMEParameters(system, force, alpha, gridSizeX, gridSizeY, gridSizeZ, false);
        gridSizeX = OpenCLFFT3D::findLegalDimension(gridSizeX);
        gridSizeY = OpenCLFFT3D::findLegalDimension(gridSizeY);
        gridSizeZ = OpenCLFFT3D::findLegalDimension(gridSizeZ);
        if (doLJPME) {
            NonbondedForceImpl::calcPMEParameters(system, force, dispersionAlpha, dispersionGridSizeX,
                                                  dispersionGridSizeY, dispersionGridSizeZ, true);
            dispersionGridSizeX = OpenCLFFT3D::findLegalDimension(dispersionGridSizeX);
            dispersionGridSizeY = OpenCLFFT3D::findLegalDimension(dispersionGridSizeY);
            dispersionGridSizeZ = OpenCLFFT3D::findLegalDimension(dispersionGridSizeZ);
        }
        defines["EWALD_ALPHA"] = cl.doubleToString(alpha);
        defines["TWO_OVER_SQRT_PI"] = cl.doubleToString(2.0/sqrt(M_PI));
        defines["USE_EWALD"] = "1";
        defines["DO_LJPME"] = doLJPME ? "1" : "0";
        if (doLJPME) {
            defines["EWALD_DISPERSION_ALPHA"] = cl.doubleToString(dispersionAlpha);
            double invRCut6 = pow(force.getCutoffDistance(), -6);
            double dalphaR = dispersionAlpha * force.getCutoffDistance();
            double dar2 = dalphaR*dalphaR;
            double dar4 = dar2*dar2;
            double multShift6 = -invRCut6*(1.0 - exp(-dar2) * (1.0 + dar2 + 0.5*dar4));
            defines["INVCUT6"] = cl.doubleToString(invRCut6);
            defines["MULTSHIFT6"] = cl.doubleToString(multShift6);
        }
        if (cl.getContextIndex() == 0) {
            paramsDefines["INCLUDE_EWALD"] = "1";
            paramsDefines["EWALD_SELF_ENERGY_SCALE"] = cl.doubleToString(ONE_4PI_EPS0*alpha/sqrt(M_PI));
            for (int i = 0; i < numParticles; i++)
                ewaldSelfEnergy -= baseParticleParamVec[i].x*baseParticleParamVec[i].x*ONE_4PI_EPS0*alpha/sqrt(M_PI);
            if (doLJPME) {
                paramsDefines["INCLUDE_LJPME"] = "1";
                paramsDefines["LJPME_SELF_ENERGY_SCALE"] = cl.doubleToString(pow(dispersionAlpha, 6)/3.0);
                for (int i = 0; i < numParticles; i++)
                    ewaldSelfEnergy += baseParticleParamVec[i].z*pow(baseParticleParamVec[i].y*dispersionAlpha, 6)/3.0;
            }
            pmeDefines["PME_ORDER"] = cl.intToString(PmeOrder);
            pmeDefines["NUM_ATOMS"] = cl.intToString(numParticles);
            pmeDefines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
            pmeDefines["RECIP_EXP_FACTOR"] = cl.doubleToString(M_PI*M_PI/(alpha*alpha));
            pmeDefines["GRID_SIZE_X"] = cl.intToString(gridSizeX);
            pmeDefines["GRID_SIZE_Y"] = cl.intToString(gridSizeY);
            pmeDefines["GRID_SIZE_Z"] = cl.intToString(gridSizeZ);
            pmeDefines["EPSILON_FACTOR"] = cl.doubleToString(sqrt(ONE_4PI_EPS0));
            pmeDefines["M_PI"] = cl.doubleToString(M_PI);
            pmeDefines["USE_FIXED_POINT_CHARGE_SPREADING"] = "1";
            bool deviceIsCpu = (cl.getDevice().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
            if (deviceIsCpu)
                pmeDefines["DEVICE_IS_CPU"] = "1";
            if (cl.getPlatformData().useCpuPme && !doLJPME && usePosqCharges) {
                // Create the CPU PME kernel.

                try {
                    cpuPme = getPlatform().createKernel(CalcPmeReciprocalForceKernel::Name(), *cl.getPlatformData().context);
                    cpuPme.getAs<CalcPmeReciprocalForceKernel>().initialize(gridSizeX, gridSizeY, gridSizeZ, numParticles, alpha, false);
                    cl::Program program = cl.createProgram(CommonKernelSources::pme, pmeDefines);
                    cl::Kernel addForcesKernel = cl::Kernel(program, "addForces");
                    pmeio = new PmeIO(cl, addForcesKernel);
                    cl.addPreComputation(new PmePreComputation(cl, cpuPme, *pmeio));
                    cl.addPostComputation(new PmePostComputation(cpuPme, *pmeio));
                }
                catch (OpenMMException& ex) {
                    // The CPU PME plugin isn't available.
                }
            }
            if (pmeio == NULL) {
                // Create required data structures.

                int elementSize = (cl.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
                int gridElements = gridSizeX*gridSizeY*gridSizeZ;
                if (doLJPME) {
                    gridElements = max(gridElements, dispersionGridSizeX*dispersionGridSizeY*dispersionGridSizeZ);
                }
                pmeGrid1.initialize(cl, gridElements, 2*elementSize, "pmeGrid1");
                pmeGrid2.initialize(cl, gridElements, 2*elementSize, "pmeGrid2");
                cl.addAutoclearBuffer(pmeGrid2);
                pmeBsplineModuliX.initialize(cl, gridSizeX, elementSize, "pmeBsplineModuliX");
                pmeBsplineModuliY.initialize(cl, gridSizeY, elementSize, "pmeBsplineModuliY");
                pmeBsplineModuliZ.initialize(cl, gridSizeZ, elementSize, "pmeBsplineModuliZ");
                if (doLJPME) {
                    pmeDispersionBsplineModuliX.initialize(cl, dispersionGridSizeX, elementSize, "pmeDispersionBsplineModuliX");
                    pmeDispersionBsplineModuliY.initialize(cl, dispersionGridSizeY, elementSize, "pmeDispersionBsplineModuliY");
                    pmeDispersionBsplineModuliZ.initialize(cl, dispersionGridSizeZ, elementSize, "pmeDispersionBsplineModuliZ");
                }
                pmeBsplineTheta.initialize(cl, PmeOrder*numParticles, 4*elementSize, "pmeBsplineTheta");
                pmeAtomRange.initialize<cl_int>(cl, gridSizeX*gridSizeY*gridSizeZ+1, "pmeAtomRange");
                pmeAtomGridIndex.initialize<mm_int2>(cl, numParticles, "pmeAtomGridIndex");
                int energyElementSize = (cl.getUseDoublePrecision() || cl.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
                pmeEnergyBuffer.initialize(cl, cl.getNumThreadBlocks()*OpenCLContext::ThreadBlockSize, energyElementSize, "pmeEnergyBuffer");
                cl.clearBuffer(pmeEnergyBuffer);
                sort = new OpenCLSort(cl, new SortTrait(), cl.getNumAtoms());
                fft = new OpenCLFFT3D(cl, gridSizeX, gridSizeY, gridSizeZ, true);
                if (doLJPME)
                    dispersionFft = new OpenCLFFT3D(cl, dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ, true);
                string vendor = cl.getDevice().getInfo<CL_DEVICE_VENDOR>();
                bool isNvidia = (vendor.size() >= 6 && vendor.substr(0, 6) == "NVIDIA");
                usePmeQueue = (!cl.getPlatformData().disablePmeStream && !cl.getPlatformData().useCpuPme && isNvidia);
                if (usePmeQueue) {
                    pmeDefines["USE_PME_STREAM"] = "1";
                    pmeQueue = cl::CommandQueue(cl.getContext(), cl.getDevice());
                    int recipForceGroup = force.getReciprocalSpaceForceGroup();
                    if (recipForceGroup < 0)
                        recipForceGroup = force.getForceGroup();
                    cl.addPreComputation(new SyncQueuePreComputation(cl, pmeQueue, recipForceGroup));
                    cl.addPostComputation(syncQueue = new SyncQueuePostComputation(cl, pmeSyncEvent, pmeEnergyBuffer, recipForceGroup));
                }

                // Initialize the b-spline moduli.

                for (int grid = 0; grid < 2; grid++) {
                    int xsize, ysize, zsize;
                    OpenCLArray *xmoduli, *ymoduli, *zmoduli;
                    if (grid == 0) {
                        xsize = gridSizeX;
                        ysize = gridSizeY;
                        zsize = gridSizeZ;
                        xmoduli = &pmeBsplineModuliX;
                        ymoduli = &pmeBsplineModuliY;
                        zmoduli = &pmeBsplineModuliZ;
                    }
                    else {
                        if (!doLJPME)
                            continue;
                        xsize = dispersionGridSizeX;
                        ysize = dispersionGridSizeY;
                        zsize = dispersionGridSizeZ;
                        xmoduli = &pmeDispersionBsplineModuliX;
                        ymoduli = &pmeDispersionBsplineModuliY;
                        zmoduli = &pmeDispersionBsplineModuliZ;
                    }
                    int maxSize = max(max(xsize, ysize), zsize);
                    vector<double> data(PmeOrder);
                    vector<double> ddata(PmeOrder);
                    vector<double> bsplines_data(maxSize);
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

                    for (int dim = 0; dim < 3; dim++) {
                        int ndata = (dim == 0 ? xsize : dim == 1 ? ysize : zsize);
                        vector<cl_double> moduli(ndata);
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
                        {
                            if (moduli[i] < 1.0e-7)
                                moduli[i] = (moduli[(i-1+ndata)%ndata]+moduli[(i+1)%ndata])*0.5;
                        }
                        if (dim == 0)
                            xmoduli->upload(moduli, true);
                        else if (dim == 1)
                            ymoduli->upload(moduli, true);
                        else
                            zmoduli->upload(moduli, true);
                    }
                }
            }
        }
    }

    // Add code to subtract off the reciprocal part of excluded interactions.

    if ((nonbondedMethod == Ewald || nonbondedMethod == PME || nonbondedMethod == LJPME) && pmeio == NULL) {
        int numContexts = cl.getPlatformData().contexts.size();
        int startIndex = cl.getContextIndex()*force.getNumExceptions()/numContexts;
        int endIndex = (cl.getContextIndex()+1)*force.getNumExceptions()/numContexts;
        int numExclusions = endIndex-startIndex;
        if (numExclusions > 0) {
            paramsDefines["HAS_EXCLUSIONS"] = "1";
            vector<vector<int> > atoms(numExclusions, vector<int>(2));
            exclusionAtoms.initialize<mm_int2>(cl, numExclusions, "exclusionAtoms");
            exclusionParams.initialize<mm_float4>(cl, numExclusions, "exclusionParams");
            vector<mm_int2> exclusionAtomsVec(numExclusions);
            for (int i = 0; i < numExclusions; i++) {
                int j = i+startIndex;
                exclusionAtomsVec[i] = mm_int2(exclusions[j].first, exclusions[j].second);
                atoms[i][0] = exclusions[j].first;
                atoms[i][1] = exclusions[j].second;
            }
            exclusionAtoms.upload(exclusionAtomsVec);
            map<string, string> replacements;
            replacements["PARAMS"] = cl.getBondedUtilities().addArgument(exclusionParams, "float4");
            replacements["EWALD_ALPHA"] = cl.doubleToString(alpha);
            replacements["TWO_OVER_SQRT_PI"] = cl.doubleToString(2.0/sqrt(M_PI));
            replacements["DO_LJPME"] = doLJPME ? "1" : "0";
            replacements["USE_PERIODIC"] = force.getExceptionsUsePeriodicBoundaryConditions() ? "1" : "0";
            if (doLJPME)
                replacements["EWALD_DISPERSION_ALPHA"] = cl.doubleToString(dispersionAlpha);
            if (force.getIncludeDirectSpace())
                cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(CommonKernelSources::pmeExclusions, replacements), force.getForceGroup());
        }
    }

    // Add the interaction to the default nonbonded kernel.
    
    string source = cl.replaceStrings(CommonKernelSources::coulombLennardJones, defines);
    charges.initialize(cl, cl.getPaddedNumAtoms(), cl.getUseDoublePrecision() ? sizeof(double) : sizeof(float), "charges");
    baseParticleParams.initialize<mm_float4>(cl, cl.getPaddedNumAtoms(), "baseParticleParams");
    baseParticleParams.upload(baseParticleParamVec);
    map<string, string> replacements;
    replacements["ONE_4PI_EPS0"] = cl.doubleToString(ONE_4PI_EPS0);
    if (usePosqCharges) {
        replacements["CHARGE1"] = "posq1.w";
        replacements["CHARGE2"] = "posq2.w";
    }
    else {
        replacements["CHARGE1"] = prefix+"charge1";
        replacements["CHARGE2"] = prefix+"charge2";
    }
    if (hasCoulomb && !usePosqCharges)
        cl.getNonbondedUtilities().addParameter(OpenCLNonbondedUtilities::ParameterInfo(prefix+"charge", "real", 1, charges.getElementSize(), charges.getDeviceBuffer()));
    sigmaEpsilon.initialize<mm_float2>(cl, cl.getPaddedNumAtoms(), "sigmaEpsilon");
    if (hasLJ) {
        replacements["SIGMA_EPSILON1"] = prefix+"sigmaEpsilon1";
        replacements["SIGMA_EPSILON2"] = prefix+"sigmaEpsilon2";
        cl.getNonbondedUtilities().addParameter(OpenCLNonbondedUtilities::ParameterInfo(prefix+"sigmaEpsilon", "float", 2, sizeof(cl_float2), sigmaEpsilon.getDeviceBuffer()));
    }
    source = cl.replaceStrings(source, replacements);
    if (force.getIncludeDirectSpace())
        cl.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, source, force.getForceGroup(), numParticles > 3000);

    // Initialize the exceptions.

    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*exceptions.size()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*exceptions.size()/numContexts;
    int numExceptions = endIndex-startIndex;
    if (numExceptions > 0) {
        paramsDefines["HAS_EXCEPTIONS"] = "1";
        exceptionAtoms.resize(numExceptions);
        vector<vector<int> > atoms(numExceptions, vector<int>(2));
        exceptionParams.initialize<mm_float4>(cl, numExceptions, "exceptionParams");
        baseExceptionParams.initialize<mm_float4>(cl, numExceptions, "baseExceptionParams");
        vector<mm_float4> baseExceptionParamsVec(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            double chargeProd, sigma, epsilon;
            force.getExceptionParameters(exceptions[startIndex+i], atoms[i][0], atoms[i][1], chargeProd, sigma, epsilon);
            baseExceptionParamsVec[i] = mm_float4(chargeProd, sigma, epsilon, 0);
            exceptionAtoms[i] = make_pair(atoms[i][0], atoms[i][1]);
        }
        baseExceptionParams.upload(baseExceptionParamsVec);
        map<string, string> replacements;
        replacements["APPLY_PERIODIC"] = (usePeriodic && force.getExceptionsUsePeriodicBoundaryConditions() ? "1" : "0");
        replacements["PARAMS"] = cl.getBondedUtilities().addArgument(exceptionParams, "float4");
        if (force.getIncludeDirectSpace())
            cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(CommonKernelSources::nonbondedExceptions, replacements), force.getForceGroup());
    }
    
    // Initialize parameter offsets.

    vector<vector<mm_float4> > particleOffsetVec(force.getNumParticles());
    vector<vector<mm_float4> > exceptionOffsetVec(numExceptions);
    for (int i = 0; i < force.getNumParticleParameterOffsets(); i++) {
        string param;
        int particle;
        double charge, sigma, epsilon;
        force.getParticleParameterOffset(i, param, particle, charge, sigma, epsilon);
        auto paramPos = find(paramNames.begin(), paramNames.end(), param);
        int paramIndex;
        if (paramPos == paramNames.end()) {
            paramIndex = paramNames.size();
            paramNames.push_back(param);
        }
        else
            paramIndex = paramPos-paramNames.begin();
        particleOffsetVec[particle].push_back(mm_float4(charge, sigma, epsilon, paramIndex));
    }
    for (int i = 0; i < force.getNumExceptionParameterOffsets(); i++) {
        string param;
        int exception;
        double charge, sigma, epsilon;
        force.getExceptionParameterOffset(i, param, exception, charge, sigma, epsilon);
        int index = exceptionIndex[exception];
        if (index < startIndex || index >= endIndex)
            continue;
        auto paramPos = find(paramNames.begin(), paramNames.end(), param);
        int paramIndex;
        if (paramPos == paramNames.end()) {
            paramIndex = paramNames.size();
            paramNames.push_back(param);
        }
        else
            paramIndex = paramPos-paramNames.begin();
        exceptionOffsetVec[index-startIndex].push_back(mm_float4(charge, sigma, epsilon, paramIndex));
    }
    paramValues.resize(paramNames.size(), 0.0);
    particleParamOffsets.initialize<mm_float4>(cl, max(force.getNumParticleParameterOffsets(), 1), "particleParamOffsets");
    particleOffsetIndices.initialize<cl_int>(cl, cl.getPaddedNumAtoms()+1, "particleOffsetIndices");
    vector<cl_int> particleOffsetIndicesVec, exceptionOffsetIndicesVec;
    vector<mm_float4> p, e;
    for (int i = 0; i < particleOffsetVec.size(); i++) {
        particleOffsetIndicesVec.push_back(p.size());
        for (int j = 0; j < particleOffsetVec[i].size(); j++)
            p.push_back(particleOffsetVec[i][j]);
    }
    while (particleOffsetIndicesVec.size() < particleOffsetIndices.getSize())
        particleOffsetIndicesVec.push_back(p.size());
    for (int i = 0; i < exceptionOffsetVec.size(); i++) {
        exceptionOffsetIndicesVec.push_back(e.size());
        for (int j = 0; j < exceptionOffsetVec[i].size(); j++)
            e.push_back(exceptionOffsetVec[i][j]);
    }
    exceptionOffsetIndicesVec.push_back(e.size());
    if (force.getNumParticleParameterOffsets() > 0) {
        particleParamOffsets.upload(p);
        particleOffsetIndices.upload(particleOffsetIndicesVec);
    }
    exceptionParamOffsets.initialize<mm_float4>(cl, max((int) e.size(), 1), "exceptionParamOffsets");
    exceptionOffsetIndices.initialize<cl_int>(cl, exceptionOffsetIndicesVec.size(), "exceptionOffsetIndices");
    if (e.size() > 0) {
        exceptionParamOffsets.upload(e);
        exceptionOffsetIndices.upload(exceptionOffsetIndicesVec);
    }
    globalParams.initialize(cl, max((int) paramValues.size(), 1), cl.getUseDoublePrecision() ? sizeof(double) : sizeof(float), "globalParams");
    if (paramValues.size() > 0)
        globalParams.upload(paramValues, true);
    recomputeParams = true;
    
    // Initialize the kernel for updating parameters.
    
    cl::Program program = cl.createProgram(CommonKernelSources::nonbondedParameters, paramsDefines);
    computeParamsKernel = cl::Kernel(program, "computeParameters");
    computeExclusionParamsKernel = cl::Kernel(program, "computeExclusionParameters");
    info = new ForceInfo(0, force);
    cl.addForce(info);
}

double OpenCLCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) {
    bool deviceIsCpu = (cl.getDevice().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        int index = 0;
        computeParamsKernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        index++;
        computeParamsKernel.setArg<cl::Buffer>(index++, globalParams.getDeviceBuffer());
        computeParamsKernel.setArg<cl_int>(index++, cl.getPaddedNumAtoms());
        computeParamsKernel.setArg<cl::Buffer>(index++, baseParticleParams.getDeviceBuffer());
        computeParamsKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        computeParamsKernel.setArg<cl::Buffer>(index++, charges.getDeviceBuffer());
        computeParamsKernel.setArg<cl::Buffer>(index++, sigmaEpsilon.getDeviceBuffer());
        computeParamsKernel.setArg<cl::Buffer>(index++, particleParamOffsets.getDeviceBuffer());
        computeParamsKernel.setArg<cl::Buffer>(index++, particleOffsetIndices.getDeviceBuffer());
        if (exceptionParams.isInitialized()) {
            computeParamsKernel.setArg<cl_int>(index++, exceptionParams.getSize());
            computeParamsKernel.setArg<cl::Buffer>(index++, baseExceptionParams.getDeviceBuffer());
            computeParamsKernel.setArg<cl::Buffer>(index++, exceptionParams.getDeviceBuffer());
            computeParamsKernel.setArg<cl::Buffer>(index++, exceptionParamOffsets.getDeviceBuffer());
            computeParamsKernel.setArg<cl::Buffer>(index++, exceptionOffsetIndices.getDeviceBuffer());
        }
        if (exclusionParams.isInitialized()) {
            computeExclusionParamsKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
            computeExclusionParamsKernel.setArg<cl::Buffer>(1, charges.getDeviceBuffer());
            computeExclusionParamsKernel.setArg<cl::Buffer>(2, sigmaEpsilon.getDeviceBuffer());
            computeExclusionParamsKernel.setArg<cl_int>(3, exclusionParams.getSize());
            computeExclusionParamsKernel.setArg<cl::Buffer>(4, exclusionAtoms.getDeviceBuffer());
            computeExclusionParamsKernel.setArg<cl::Buffer>(5, exclusionParams.getDeviceBuffer());
        }
        if (cosSinSums.isInitialized()) {
            ewaldSumsKernel.setArg<cl::Buffer>(0, cl.getEnergyBuffer().getDeviceBuffer());
            ewaldSumsKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
            ewaldSumsKernel.setArg<cl::Buffer>(2, cosSinSums.getDeviceBuffer());
            ewaldForcesKernel.setArg<cl::Buffer>(0, cl.getLongForceBuffer().getDeviceBuffer());
            ewaldForcesKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
            ewaldForcesKernel.setArg<cl::Buffer>(2, cosSinSums.getDeviceBuffer());
        }
        if (pmeGrid1.isInitialized()) {
            // Create kernels for Coulomb PME.
            
            map<string, string> replacements;
            replacements["CHARGE"] = (usePosqCharges ? "pos.w" : "charges[atom]");
            cl::Program program = cl.createProgram(cl.replaceStrings(CommonKernelSources::pme, replacements), pmeDefines);
            pmeGridIndexKernel = cl::Kernel(program, "findAtomGridIndex");
            pmeSpreadChargeKernel = cl::Kernel(program, "gridSpreadCharge");
            pmeConvolutionKernel = cl::Kernel(program, "reciprocalConvolution");
            pmeEvalEnergyKernel = cl::Kernel(program, "gridEvaluateEnergy");
            pmeInterpolateForceKernel = cl::Kernel(program, "gridInterpolateForce");
            int elementSize = (cl.getUseDoublePrecision() ? sizeof(mm_double4) : sizeof(mm_float4));
            pmeGridIndexKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
            pmeGridIndexKernel.setArg<cl::Buffer>(1, pmeAtomGridIndex.getDeviceBuffer());
            pmeSpreadChargeKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
            pmeSpreadChargeKernel.setArg<cl::Buffer>(1, pmeGrid2.getDeviceBuffer());
            pmeSpreadChargeKernel.setArg<cl::Buffer>(10, pmeAtomGridIndex.getDeviceBuffer());
            pmeSpreadChargeKernel.setArg<cl::Buffer>(11, charges.getDeviceBuffer());
            pmeConvolutionKernel.setArg<cl::Buffer>(0, pmeGrid2.getDeviceBuffer());
            pmeConvolutionKernel.setArg<cl::Buffer>(1, pmeBsplineModuliX.getDeviceBuffer());
            pmeConvolutionKernel.setArg<cl::Buffer>(2, pmeBsplineModuliY.getDeviceBuffer());
            pmeConvolutionKernel.setArg<cl::Buffer>(3, pmeBsplineModuliZ.getDeviceBuffer());
            pmeEvalEnergyKernel.setArg<cl::Buffer>(0, pmeGrid2.getDeviceBuffer());
            pmeEvalEnergyKernel.setArg<cl::Buffer>(1, usePmeQueue ? pmeEnergyBuffer.getDeviceBuffer() : cl.getEnergyBuffer().getDeviceBuffer());
            pmeEvalEnergyKernel.setArg<cl::Buffer>(2, pmeBsplineModuliX.getDeviceBuffer());
            pmeEvalEnergyKernel.setArg<cl::Buffer>(3, pmeBsplineModuliY.getDeviceBuffer());
            pmeEvalEnergyKernel.setArg<cl::Buffer>(4, pmeBsplineModuliZ.getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(1, cl.getLongForceBuffer().getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(2, pmeGrid1.getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(11, pmeAtomGridIndex.getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(12, charges.getDeviceBuffer());
            pmeFinishSpreadChargeKernel = cl::Kernel(program, "finishSpreadCharge");
            pmeFinishSpreadChargeKernel.setArg<cl::Buffer>(0, pmeGrid2.getDeviceBuffer());
            pmeFinishSpreadChargeKernel.setArg<cl::Buffer>(1, pmeGrid1.getDeviceBuffer());
            if (usePmeQueue)
                syncQueue->setKernel(cl::Kernel(program, "addEnergy"));

            if (doLJPME) {
                // Create kernels for LJ PME.

                pmeDefines["EWALD_ALPHA"] = cl.doubleToString(dispersionAlpha);
                pmeDefines["GRID_SIZE_X"] = cl.intToString(dispersionGridSizeX);
                pmeDefines["GRID_SIZE_Y"] = cl.intToString(dispersionGridSizeY);
                pmeDefines["GRID_SIZE_Z"] = cl.intToString(dispersionGridSizeZ);
                pmeDefines["EPSILON_FACTOR"] = "1";
                pmeDefines["RECIP_EXP_FACTOR"] = cl.doubleToString(M_PI*M_PI/(dispersionAlpha*dispersionAlpha));
                pmeDefines["USE_LJPME"] = "1";
                pmeDefines["CHARGE_FROM_SIGEPS"] = "1";
                program = cl.createProgram(CommonKernelSources::pme, pmeDefines);
                pmeDispersionGridIndexKernel = cl::Kernel(program, "findAtomGridIndex");
                pmeDispersionSpreadChargeKernel = cl::Kernel(program, "gridSpreadCharge");
                pmeDispersionConvolutionKernel = cl::Kernel(program, "reciprocalConvolution");
                pmeDispersionEvalEnergyKernel = cl::Kernel(program, "gridEvaluateEnergy");
                pmeDispersionInterpolateForceKernel = cl::Kernel(program, "gridInterpolateForce");
                pmeDispersionGridIndexKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
                pmeDispersionGridIndexKernel.setArg<cl::Buffer>(1, pmeAtomGridIndex.getDeviceBuffer());
                pmeDispersionSpreadChargeKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
                pmeDispersionSpreadChargeKernel.setArg<cl::Buffer>(1, pmeGrid2.getDeviceBuffer());
                pmeDispersionSpreadChargeKernel.setArg<cl::Buffer>(10, pmeAtomGridIndex.getDeviceBuffer());
                pmeDispersionSpreadChargeKernel.setArg<cl::Buffer>(11, sigmaEpsilon.getDeviceBuffer());
                pmeDispersionConvolutionKernel.setArg<cl::Buffer>(0, pmeGrid2.getDeviceBuffer());
                pmeDispersionConvolutionKernel.setArg<cl::Buffer>(1, pmeDispersionBsplineModuliX.getDeviceBuffer());
                pmeDispersionConvolutionKernel.setArg<cl::Buffer>(2, pmeDispersionBsplineModuliY.getDeviceBuffer());
                pmeDispersionConvolutionKernel.setArg<cl::Buffer>(3, pmeDispersionBsplineModuliZ.getDeviceBuffer());
                pmeDispersionEvalEnergyKernel.setArg<cl::Buffer>(0, pmeGrid2.getDeviceBuffer());
                pmeDispersionEvalEnergyKernel.setArg<cl::Buffer>(1, usePmeQueue ? pmeEnergyBuffer.getDeviceBuffer() : cl.getEnergyBuffer().getDeviceBuffer());
                pmeDispersionEvalEnergyKernel.setArg<cl::Buffer>(2, pmeDispersionBsplineModuliX.getDeviceBuffer());
                pmeDispersionEvalEnergyKernel.setArg<cl::Buffer>(3, pmeDispersionBsplineModuliY.getDeviceBuffer());
                pmeDispersionEvalEnergyKernel.setArg<cl::Buffer>(4, pmeDispersionBsplineModuliZ.getDeviceBuffer());
                pmeDispersionInterpolateForceKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
                pmeDispersionInterpolateForceKernel.setArg<cl::Buffer>(1, cl.getLongForceBuffer().getDeviceBuffer());
                pmeDispersionInterpolateForceKernel.setArg<cl::Buffer>(2, pmeGrid1.getDeviceBuffer());
                pmeDispersionInterpolateForceKernel.setArg<cl::Buffer>(11, pmeAtomGridIndex.getDeviceBuffer());
                pmeDispersionInterpolateForceKernel.setArg<cl::Buffer>(12, sigmaEpsilon.getDeviceBuffer());
                pmeDispersionFinishSpreadChargeKernel = cl::Kernel(program, "finishSpreadCharge");
                pmeDispersionFinishSpreadChargeKernel.setArg<cl::Buffer>(0, pmeGrid2.getDeviceBuffer());
                pmeDispersionFinishSpreadChargeKernel.setArg<cl::Buffer>(1, pmeGrid1.getDeviceBuffer());
            }
       }
    }
    
    // Update particle and exception parameters.

    bool paramChanged = false;
    for (int i = 0; i < paramNames.size(); i++) {
        double value = context.getParameter(paramNames[i]);
        if (value != paramValues[i]) {
            paramValues[i] = value;;
            paramChanged = true;
        }
    }
    if (paramChanged) {
        recomputeParams = true;
        globalParams.upload(paramValues, true);
    }
    double energy = (includeReciprocal ? ewaldSelfEnergy : 0.0);
    if (recomputeParams || hasOffsets) {
        computeParamsKernel.setArg<cl_int>(1, includeEnergy && includeReciprocal);
        cl.executeKernel(computeParamsKernel, cl.getPaddedNumAtoms());
        if (exclusionParams.isInitialized())
            cl.executeKernel(computeExclusionParamsKernel, exclusionParams.getSize());
        if (usePmeQueue) {
            vector<cl::Event> events(1);
            cl.getQueue().enqueueMarkerWithWaitList(NULL, &events[0]);
            pmeQueue.enqueueBarrierWithWaitList(&events);
        }
        if (hasOffsets)
            energy = 0.0; // The Ewald self energy was computed in the kernel.
        recomputeParams = false;
    }
    
    // Do reciprocal space calculations.
    
    if (cosSinSums.isInitialized() && includeReciprocal) {
        mm_double4 boxSize = cl.getPeriodicBoxSizeDouble();
        if (cl.getUseDoublePrecision()) {
            ewaldSumsKernel.setArg<mm_double4>(3, boxSize);
            ewaldForcesKernel.setArg<mm_double4>(3, boxSize);
        }
        else {
            ewaldSumsKernel.setArg<mm_float4>(3, mm_float4((float) boxSize.x, (float) boxSize.y, (float) boxSize.z, 0));
            ewaldForcesKernel.setArg<mm_float4>(3, mm_float4((float) boxSize.x, (float) boxSize.y, (float) boxSize.z, 0));
        }
        cl.executeKernel(ewaldSumsKernel, cosSinSums.getSize());
        cl.executeKernel(ewaldForcesKernel, cl.getNumAtoms());
    }
    if (pmeGrid1.isInitialized() && includeReciprocal) {
        if (usePmeQueue && !includeEnergy)
            cl.setQueue(pmeQueue);
        
        // Invert the periodic box vectors.
        
        Vec3 boxVectors[3];
        cl.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double determinant = boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2];
        double scale = 1.0/determinant;
        mm_double4 recipBoxVectors[3];
        recipBoxVectors[0] = mm_double4(boxVectors[1][1]*boxVectors[2][2]*scale, 0, 0, 0);
        recipBoxVectors[1] = mm_double4(-boxVectors[1][0]*boxVectors[2][2]*scale, boxVectors[0][0]*boxVectors[2][2]*scale, 0, 0);
        recipBoxVectors[2] = mm_double4((boxVectors[1][0]*boxVectors[2][1]-boxVectors[1][1]*boxVectors[2][0])*scale, -boxVectors[0][0]*boxVectors[2][1]*scale, boxVectors[0][0]*boxVectors[1][1]*scale, 0);
        mm_float4 recipBoxVectorsFloat[3];
        for (int i = 0; i < 3; i++)
            recipBoxVectorsFloat[i] = mm_float4((float) recipBoxVectors[i].x, (float) recipBoxVectors[i].y, (float) recipBoxVectors[i].z, 0);
        
        // Execute the reciprocal space kernels.

        if (hasCoulomb) {
            setPeriodicBoxArgs(cl, pmeGridIndexKernel, 2);
            if (cl.getUseDoublePrecision()) {
                pmeGridIndexKernel.setArg<mm_double4>(7, recipBoxVectors[0]);
                pmeGridIndexKernel.setArg<mm_double4>(8, recipBoxVectors[1]);
                pmeGridIndexKernel.setArg<mm_double4>(9, recipBoxVectors[2]);
            }
            else {
                pmeGridIndexKernel.setArg<mm_float4>(7, recipBoxVectorsFloat[0]);
                pmeGridIndexKernel.setArg<mm_float4>(8, recipBoxVectorsFloat[1]);
                pmeGridIndexKernel.setArg<mm_float4>(9, recipBoxVectorsFloat[2]);
            }
            cl.executeKernel(pmeGridIndexKernel, cl.getNumAtoms());
            sort->sort(pmeAtomGridIndex);
            setPeriodicBoxArgs(cl, pmeSpreadChargeKernel, 2);
            if (cl.getUseDoublePrecision()) {
                pmeSpreadChargeKernel.setArg<mm_double4>(7, recipBoxVectors[0]);
                pmeSpreadChargeKernel.setArg<mm_double4>(8, recipBoxVectors[1]);
                pmeSpreadChargeKernel.setArg<mm_double4>(9, recipBoxVectors[2]);
            }
            else {
                pmeSpreadChargeKernel.setArg<mm_float4>(7, recipBoxVectorsFloat[0]);
                pmeSpreadChargeKernel.setArg<mm_float4>(8, recipBoxVectorsFloat[1]);
                pmeSpreadChargeKernel.setArg<mm_float4>(9, recipBoxVectorsFloat[2]);
            }
            cl.executeKernel(pmeSpreadChargeKernel, cl.getNumAtoms());
            cl.executeKernel(pmeFinishSpreadChargeKernel, gridSizeX*gridSizeY*gridSizeZ);
            fft->execFFT(pmeGrid1, pmeGrid2, true);
            mm_double4 boxSize = cl.getPeriodicBoxSizeDouble();
            if (cl.getUseDoublePrecision()) {
                pmeConvolutionKernel.setArg<mm_double4>(4, recipBoxVectors[0]);
                pmeConvolutionKernel.setArg<mm_double4>(5, recipBoxVectors[1]);
                pmeConvolutionKernel.setArg<mm_double4>(6, recipBoxVectors[2]);
                pmeEvalEnergyKernel.setArg<mm_double4>(5, recipBoxVectors[0]);
                pmeEvalEnergyKernel.setArg<mm_double4>(6, recipBoxVectors[1]);
                pmeEvalEnergyKernel.setArg<mm_double4>(7, recipBoxVectors[2]);
            }
            else {
                pmeConvolutionKernel.setArg<mm_float4>(4, recipBoxVectorsFloat[0]);
                pmeConvolutionKernel.setArg<mm_float4>(5, recipBoxVectorsFloat[1]);
                pmeConvolutionKernel.setArg<mm_float4>(6, recipBoxVectorsFloat[2]);
                pmeEvalEnergyKernel.setArg<mm_float4>(5, recipBoxVectorsFloat[0]);
                pmeEvalEnergyKernel.setArg<mm_float4>(6, recipBoxVectorsFloat[1]);
                pmeEvalEnergyKernel.setArg<mm_float4>(7, recipBoxVectorsFloat[2]);
            }
            if (includeEnergy)
                cl.executeKernel(pmeEvalEnergyKernel, gridSizeX*gridSizeY*gridSizeZ);
            cl.executeKernel(pmeConvolutionKernel, gridSizeX*gridSizeY*gridSizeZ);
            fft->execFFT(pmeGrid2, pmeGrid1, false);
            setPeriodicBoxArgs(cl, pmeInterpolateForceKernel, 3);
            if (cl.getUseDoublePrecision()) {
                pmeInterpolateForceKernel.setArg<mm_double4>(8, recipBoxVectors[0]);
                pmeInterpolateForceKernel.setArg<mm_double4>(9, recipBoxVectors[1]);
                pmeInterpolateForceKernel.setArg<mm_double4>(10, recipBoxVectors[2]);
            }
            else {
                pmeInterpolateForceKernel.setArg<mm_float4>(8, recipBoxVectorsFloat[0]);
                pmeInterpolateForceKernel.setArg<mm_float4>(9, recipBoxVectorsFloat[1]);
                pmeInterpolateForceKernel.setArg<mm_float4>(10, recipBoxVectorsFloat[2]);
            }
            if (deviceIsCpu)
                cl.executeKernel(pmeInterpolateForceKernel, 2*cl.getDevice().getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>(), 1);
            else
                cl.executeKernel(pmeInterpolateForceKernel, cl.getNumAtoms());
        }
        
        if (doLJPME && hasLJ) {
            setPeriodicBoxArgs(cl, pmeDispersionGridIndexKernel, 2);
            if (cl.getUseDoublePrecision()) {
                pmeDispersionGridIndexKernel.setArg<mm_double4>(7, recipBoxVectors[0]);
                pmeDispersionGridIndexKernel.setArg<mm_double4>(8, recipBoxVectors[1]);
                pmeDispersionGridIndexKernel.setArg<mm_double4>(9, recipBoxVectors[2]);
            }
            else {
                pmeDispersionGridIndexKernel.setArg<mm_float4>(7, recipBoxVectorsFloat[0]);
                pmeDispersionGridIndexKernel.setArg<mm_float4>(8, recipBoxVectorsFloat[1]);
                pmeDispersionGridIndexKernel.setArg<mm_float4>(9, recipBoxVectorsFloat[2]);
            }
            cl.executeKernel(pmeDispersionGridIndexKernel, cl.getNumAtoms());
            if (!hasCoulomb)
                sort->sort(pmeAtomGridIndex);
            cl.clearBuffer(pmeGrid2);
            setPeriodicBoxArgs(cl, pmeDispersionSpreadChargeKernel, 2);
            if (cl.getUseDoublePrecision()) {
                pmeDispersionSpreadChargeKernel.setArg<mm_double4>(7, recipBoxVectors[0]);
                pmeDispersionSpreadChargeKernel.setArg<mm_double4>(8, recipBoxVectors[1]);
                pmeDispersionSpreadChargeKernel.setArg<mm_double4>(9, recipBoxVectors[2]);
            }
            else {
                pmeDispersionSpreadChargeKernel.setArg<mm_float4>(7, recipBoxVectorsFloat[0]);
                pmeDispersionSpreadChargeKernel.setArg<mm_float4>(8, recipBoxVectorsFloat[1]);
                pmeDispersionSpreadChargeKernel.setArg<mm_float4>(9, recipBoxVectorsFloat[2]);
            }
            cl.executeKernel(pmeDispersionSpreadChargeKernel, cl.getNumAtoms());
            cl.executeKernel(pmeDispersionFinishSpreadChargeKernel, gridSizeX*gridSizeY*gridSizeZ);
            dispersionFft->execFFT(pmeGrid1, pmeGrid2, true);
            if (cl.getUseDoublePrecision()) {
                pmeDispersionConvolutionKernel.setArg<mm_double4>(4, recipBoxVectors[0]);
                pmeDispersionConvolutionKernel.setArg<mm_double4>(5, recipBoxVectors[1]);
                pmeDispersionConvolutionKernel.setArg<mm_double4>(6, recipBoxVectors[2]);
                pmeDispersionEvalEnergyKernel.setArg<mm_double4>(5, recipBoxVectors[0]);
                pmeDispersionEvalEnergyKernel.setArg<mm_double4>(6, recipBoxVectors[1]);
                pmeDispersionEvalEnergyKernel.setArg<mm_double4>(7, recipBoxVectors[2]);
            }
            else {
                pmeDispersionConvolutionKernel.setArg<mm_float4>(4, recipBoxVectorsFloat[0]);
                pmeDispersionConvolutionKernel.setArg<mm_float4>(5, recipBoxVectorsFloat[1]);
                pmeDispersionConvolutionKernel.setArg<mm_float4>(6, recipBoxVectorsFloat[2]);
                pmeDispersionEvalEnergyKernel.setArg<mm_float4>(5, recipBoxVectorsFloat[0]);
                pmeDispersionEvalEnergyKernel.setArg<mm_float4>(6, recipBoxVectorsFloat[1]);
                pmeDispersionEvalEnergyKernel.setArg<mm_float4>(7, recipBoxVectorsFloat[2]);
            }
            if (!hasCoulomb) cl.clearBuffer(pmeEnergyBuffer);
            if (includeEnergy)
                cl.executeKernel(pmeDispersionEvalEnergyKernel, gridSizeX*gridSizeY*gridSizeZ);
            cl.executeKernel(pmeDispersionConvolutionKernel, gridSizeX*gridSizeY*gridSizeZ);
            dispersionFft->execFFT(pmeGrid2, pmeGrid1, false);
            setPeriodicBoxArgs(cl, pmeDispersionInterpolateForceKernel, 3);
            if (cl.getUseDoublePrecision()) {
                pmeDispersionInterpolateForceKernel.setArg<mm_double4>(8, recipBoxVectors[0]);
                pmeDispersionInterpolateForceKernel.setArg<mm_double4>(9, recipBoxVectors[1]);
                pmeDispersionInterpolateForceKernel.setArg<mm_double4>(10, recipBoxVectors[2]);
            }
            else {
                pmeDispersionInterpolateForceKernel.setArg<mm_float4>(8, recipBoxVectorsFloat[0]);
                pmeDispersionInterpolateForceKernel.setArg<mm_float4>(9, recipBoxVectorsFloat[1]);
                pmeDispersionInterpolateForceKernel.setArg<mm_float4>(10, recipBoxVectorsFloat[2]);
            }
            if (deviceIsCpu)
                cl.executeKernel(pmeDispersionInterpolateForceKernel, 2*cl.getDevice().getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>(), 1);
            else
                cl.executeKernel(pmeDispersionInterpolateForceKernel, cl.getNumAtoms());
        }
        if (usePmeQueue) {
            pmeQueue.enqueueMarkerWithWaitList(NULL, &pmeSyncEvent);
            cl.restoreDefaultQueue();
        }
    }
    if (dispersionCoefficient != 0.0 && includeDirect) {
        mm_double4 boxSize = cl.getPeriodicBoxSizeDouble();
        energy += dispersionCoefficient/(boxSize.x*boxSize.y*boxSize.z);
    }
    return energy;
}

void OpenCLCalcNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const NonbondedForce& force, int firstParticle, int lastParticle, int firstException, int lastException) {
    // Make sure the new parameters are acceptable.

    if (force.getNumParticles() != cl.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    if (!hasCoulomb || !hasLJ) {
        for (int i = 0; i < force.getNumParticles(); i++) {
            double charge, sigma, epsilon;
            force.getParticleParameters(i, charge, sigma, epsilon);
            if (!hasCoulomb && charge != 0.0)
                throw OpenMMException("updateParametersInContext: The nonbonded force kernel does not include Coulomb interactions, because all charges were originally 0");
            if (!hasLJ && epsilon != 0.0)
                throw OpenMMException("updateParametersInContext: The nonbonded force kernel does not include Lennard-Jones interactions, because all epsilons were originally 0");
        }
    }
    set<int> exceptionsWithOffsets;
    for (int i = 0; i < force.getNumExceptionParameterOffsets(); i++) {
        string param;
        int exception;
        double charge, sigma, epsilon;
        force.getExceptionParameterOffset(i, param, exception, charge, sigma, epsilon);
        exceptionsWithOffsets.insert(exception);
    }
    vector<int> exceptions;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        if (chargeProd != 0.0 || epsilon != 0.0 || exceptionsWithOffsets.find(i) != exceptionsWithOffsets.end())
            exceptions.push_back(i);
    }
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*exceptions.size()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*exceptions.size()/numContexts;
    int numExceptions = endIndex-startIndex;
    if (numExceptions != exceptionAtoms.size())
        throw OpenMMException("updateParametersInContext: The set of non-excluded exceptions has changed");

    // Record the per-particle parameters.

    if (firstParticle <= lastParticle) {
        vector<mm_float4> baseParticleParamVec(cl.getPaddedNumAtoms(), mm_float4(0, 0, 0, 0));
        for (int i = 0; i < force.getNumParticles(); i++) {
            double charge, sigma, epsilon;
            force.getParticleParameters(i, charge, sigma, epsilon);
            baseParticleParamVec[i] = mm_float4(charge, sigma, epsilon, 0);
        }
        baseParticleParams.uploadSubArray(&baseParticleParamVec[firstParticle], firstParticle, lastParticle-firstParticle+1);

        // Compute the self energy.

        ewaldSelfEnergy = 0.0;
        if (nonbondedMethod == Ewald || nonbondedMethod == PME || nonbondedMethod == LJPME) {
            if (cl.getContextIndex() == 0) {
                for (int i = 0; i < force.getNumParticles(); i++) {
                    ewaldSelfEnergy -= baseParticleParamVec[i].x*baseParticleParamVec[i].x*ONE_4PI_EPS0*alpha/sqrt(M_PI);
                    if (doLJPME)
                        ewaldSelfEnergy += baseParticleParamVec[i].z*pow(baseParticleParamVec[i].y*dispersionAlpha, 6)/3.0;
                }
            }
        }
    }
    
    // Record the exceptions.

    if (firstException <= lastException) {
        vector<mm_float4> baseExceptionParamsVec(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            int particle1, particle2;
            double chargeProd, sigma, epsilon;
            force.getExceptionParameters(exceptions[startIndex+i], particle1, particle2, chargeProd, sigma, epsilon);
            if (make_pair(particle1, particle2) != exceptionAtoms[i])
                throw OpenMMException("updateParametersInContext: The set of non-excluded exceptions has changed");
            baseExceptionParamsVec[i] = mm_float4(chargeProd, sigma, epsilon, 0);
        }
        baseExceptionParams.upload(baseExceptionParamsVec);
    }
    
    // Compute other values.
    
    if (force.getUseDispersionCorrection() && cl.getContextIndex() == 0 && (nonbondedMethod == CutoffPeriodic || nonbondedMethod == Ewald || nonbondedMethod == PME))
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(context.getSystem(), force);
    cl.invalidateMolecules(info, firstParticle <= lastParticle, firstException <= lastException);
    recomputeParams = true;
}

void OpenCLCalcNonbondedForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (nonbondedMethod != PME)
        throw OpenMMException("getPMEParametersInContext: This Context is not using PME");
    if (cl.getPlatformData().useCpuPme)
        cpuPme.getAs<CalcPmeReciprocalForceKernel>().getPMEParameters(alpha, nx, ny, nz);
    else {
        alpha = this->alpha;
        nx = gridSizeX;
        ny = gridSizeY;
        nz = gridSizeZ;
    }
}

void OpenCLCalcNonbondedForceKernel::getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (nonbondedMethod != LJPME)
        throw OpenMMException("getPMEParametersInContext: This Context is not using PME");
    if (cl.getPlatformData().useCpuPme)
        //cpuPme.getAs<CalcPmeReciprocalForceKernel>().getLJPMEParameters(alpha, nx, ny, nz);
        throw OpenMMException("getPMEParametersInContext: CPUPME has not been implemented for LJPME yet.");
    else {
        alpha = this->dispersionAlpha;
        nx = dispersionGridSizeX;
        ny = dispersionGridSizeY;
        nz = dispersionGridSizeZ;
    }
}

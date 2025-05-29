/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2024 Stanford University and the Authors.      *
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

#include "CommonDrudeKernels.h"
#include "CommonDrudeKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/common/BondedUtilities.h"
#include "openmm/common/ComputeForceInfo.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/common/IntegrationUtilities.h"
#include "CommonKernelSources.h"
#include "SimTKOpenMMRealType.h"
#include <set>

using namespace OpenMM;
using namespace std;

class CommonDrudeForceInfo : public ComputeForceInfo {
public:
    CommonDrudeForceInfo(const DrudeForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumParticles()+force.getNumScreenedPairs();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        particles.clear();
        if (index < force.getNumParticles()) {
            int p, p1, p2, p3, p4;
            double charge, polarizability, aniso12, aniso34;
            force.getParticleParameters(index, p, p1, p2, p3, p4, charge, polarizability, aniso12, aniso34);
            particles.push_back(p);
            particles.push_back(p1);
            if (p2 != -1)
                particles.push_back(p2);
            if (p3 != -1)
                particles.push_back(p3);
            if (p4 != -1)
                particles.push_back(p4);
        }
        else {
            int drude1, drude2;
            double thole;
            force.getScreenedPairParameters(index-force.getNumParticles(), drude1, drude2, thole);
            int p, p1, p2, p3, p4;
            double charge, polarizability, aniso12, aniso34;
            force.getParticleParameters(drude1, p, p1, p2, p3, p4, charge, polarizability, aniso12, aniso34);
            particles.push_back(p);
            particles.push_back(p1);
            force.getParticleParameters(drude2, p, p1, p2, p3, p4, charge, polarizability, aniso12, aniso34);
            particles.push_back(p);
            particles.push_back(p1);
        }
    }
    bool areGroupsIdentical(int group1, int group2) {
        if (group1 < force.getNumParticles() && group2 < force.getNumParticles()) {
            int p, p1, p2, p3, p4;
            double charge1, polarizability1, aniso12_1, aniso34_1;
            double charge2, polarizability2, aniso12_2, aniso34_2;
            force.getParticleParameters(group1, p, p1, p2, p3, p4, charge1, polarizability1, aniso12_1, aniso34_1);
            force.getParticleParameters(group2, p, p1, p2, p3, p4, charge2, polarizability2, aniso12_2, aniso34_2);
            return (charge1 == charge2 && polarizability1 == polarizability2 && aniso12_1 == aniso12_2 && aniso34_1 == aniso34_2);
        }
        if (group1 >= force.getNumParticles() && group2 >= force.getNumParticles()) {
            int drude1, drude2;
            double thole1, thole2;
            force.getScreenedPairParameters(group1-force.getNumParticles(), drude1, drude2, thole1);
            force.getScreenedPairParameters(group1-force.getNumParticles(), drude1, drude2, thole2);
            return (thole1 == thole2);
        }
        return false;
    }
private:
    const DrudeForce& force;
};

void CommonCalcDrudeForceKernel::initialize(const System& system, const DrudeForce& force) {
    if (cc.getContextIndex() != 0)
        return; // This is run entirely on one device
    ContextSelector selector(cc);
    int numParticles = force.getNumParticles();
    if (numParticles > 0) {
        // Create the harmonic interaction .
        
        vector<vector<int> > atoms(numParticles, vector<int>(5));
        particleParams.initialize<mm_float4>(cc, numParticles, "drudeParticleParams");
        vector<mm_float4> paramVector(numParticles);
        for (int i = 0; i < numParticles; i++) {
            double charge, polarizability, aniso12, aniso34;
            force.getParticleParameters(i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], atoms[i][4], charge, polarizability, aniso12, aniso34);
            double a1 = (atoms[i][2] == -1 ? 1 : aniso12);
            double a2 = (atoms[i][3] == -1 || atoms[i][4] == -1 ? 1 : aniso34);
            double a3 = 3-a1-a2;
            double k3 = ONE_4PI_EPS0*charge*charge/(polarizability*a3);
            double k1 = ONE_4PI_EPS0*charge*charge/(polarizability*a1) - k3;
            double k2 = ONE_4PI_EPS0*charge*charge/(polarizability*a2) - k3;
            if (atoms[i][2] == -1) {
                atoms[i][2] = atoms[i][0];
                k1 = 0;
            }
            if (atoms[i][3] == -1 || atoms[i][4] == -1) {
                atoms[i][3] = atoms[i][0];
                atoms[i][4] = atoms[i][0];
                k2 = 0;
            }
            paramVector[i] = mm_float4((float) k1, (float) k2, (float) k3, 0.0f);
        }
        particleParams.upload(paramVector);
        map<string, string> replacements;
        replacements["PARAMS"] = cc.getBondedUtilities().addArgument(particleParams, "float4");
        cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonDrudeKernelSources::drudeParticleForce, replacements), force.getForceGroup());
    }
    int numPairs = force.getNumScreenedPairs();
    if (numPairs > 0) {
        // Create the screened interaction between dipole pairs.
        
        vector<vector<int> > atoms(numPairs, vector<int>(4));
        pairParams.initialize<mm_float2>(cc, numPairs, "drudePairParams");
        vector<mm_float2> paramVector(numPairs);
        for (int i = 0; i < numPairs; i++) {
            int drude1, drude2;
            double thole;
            force.getScreenedPairParameters(i, drude1, drude2, thole);
            int p2, p3, p4;
            double charge1, charge2, polarizability1, polarizability2, aniso12, aniso34;
            force.getParticleParameters(drude1, atoms[i][0], atoms[i][1], p2, p3, p4, charge1, polarizability1, aniso12, aniso34);
            force.getParticleParameters(drude2, atoms[i][2], atoms[i][3], p2, p3, p4, charge2, polarizability2, aniso12, aniso34);
            double screeningScale = thole/pow(polarizability1*polarizability2, 1.0/6.0);
            double energyScale = ONE_4PI_EPS0*charge1*charge2;
            paramVector[i] = mm_float2((float) screeningScale, (float) energyScale);
        }
        pairParams.upload(paramVector);
        map<string, string> replacements;
        replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
        replacements["PARAMS"] = cc.getBondedUtilities().addArgument(pairParams, "float2");
        cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonDrudeKernelSources::drudePairForce, replacements), force.getForceGroup());
    }
    cc.addForce(new CommonDrudeForceInfo(force));
}

double CommonCalcDrudeForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CommonCalcDrudeForceKernel::copyParametersToContext(ContextImpl& context, const DrudeForce& force) {
    if (cc.getContextIndex() != 0)
        return; // This is run entirely on one device
    
    // Set the particle parameters.
    
    ContextSelector selector(cc);
    int numParticles = force.getNumParticles();
    if (numParticles > 0) {
        if (!particleParams.isInitialized() || numParticles != particleParams.getSize())
            throw OpenMMException("updateParametersInContext: The number of Drude particles has changed");
        vector<mm_float4> paramVector(numParticles);
        for (int i = 0; i < numParticles; i++) {
            int p, p1, p2, p3, p4;
            double charge, polarizability, aniso12, aniso34;
            force.getParticleParameters(i, p, p1, p2, p3, p4, charge, polarizability, aniso12, aniso34);
            double a1 = (p2 == -1 ? 1 : aniso12);
            double a2 = (p3 == -1 || p4 == -1 ? 1 : aniso34);
            double a3 = 3-a1-a2;
            double k3 = ONE_4PI_EPS0*charge*charge/(polarizability*a3);
            double k1 = ONE_4PI_EPS0*charge*charge/(polarizability*a1) - k3;
            double k2 = ONE_4PI_EPS0*charge*charge/(polarizability*a2) - k3;
            if (p2 == -1)
                k1 = 0;
            if (p3 == -1 || p4 == -1)
                k2 = 0;
            paramVector[i] = mm_float4((float) k1, (float) k2, (float) k3, 0.0f);
        }
        particleParams.upload(paramVector);
    }
    
    // Set the pair parameters.
    
    int numPairs = force.getNumScreenedPairs();
    if (numPairs > 0) {
        if (!pairParams.isInitialized() || numPairs != pairParams.getSize())
            throw OpenMMException("updateParametersInContext: The number of screened pairs has changed");
        vector<mm_float2> paramVector(numPairs);
        for (int i = 0; i < numPairs; i++) {
            int drude1, drude2;
            double thole;
            force.getScreenedPairParameters(i, drude1, drude2, thole);
            int p, p1, p2, p3, p4;
            double charge1, charge2, polarizability1, polarizability2, aniso12, aniso34;
            force.getParticleParameters(drude1, p, p1, p2, p3, p4, charge1, polarizability1, aniso12, aniso34);
            force.getParticleParameters(drude2, p, p1, p2, p3, p4, charge2, polarizability2, aniso12, aniso34);
            double screeningScale = thole/pow(polarizability1*polarizability2, 1.0/6.0);
            double energyScale = ONE_4PI_EPS0*charge1*charge2;
            paramVector[i] = mm_float2((float) screeningScale, (float) energyScale);
        }
        pairParams.upload(paramVector);
    }
}

void CommonIntegrateDrudeLangevinStepKernel::initialize(const System& system, const DrudeLangevinIntegrator& integrator, const DrudeForce& force) {
    cc.initializeContexts();
    ContextSelector selector(cc);
    cc.getIntegrationUtilities().initRandomNumberGenerator((unsigned int) integrator.getRandomNumberSeed());
    
    // Identify particle pairs and ordinary particles.
    
    set<int> particles;
    vector<int> normalParticleVec;
    vector<mm_int2> pairParticleVec;
    for (int i = 0; i < system.getNumParticles(); i++)
        particles.insert(i);
    for (int i = 0; i < force.getNumParticles(); i++) {
        int p, p1, p2, p3, p4;
        double charge, polarizability, aniso12, aniso34;
        force.getParticleParameters(i, p, p1, p2, p3, p4, charge, polarizability, aniso12, aniso34);
        particles.erase(p);
        particles.erase(p1);
        pairParticleVec.push_back(mm_int2(p, p1));
    }
    normalParticleVec.insert(normalParticleVec.begin(), particles.begin(), particles.end());
    normalParticles.initialize<int>(cc, max((int) normalParticleVec.size(), 1), "drudeNormalParticles");
    pairParticles.initialize<mm_int2>(cc, max((int) pairParticleVec.size(), 1), "drudePairParticles");
    if (normalParticleVec.size() > 0)
        normalParticles.upload(normalParticleVec);
    if (pairParticleVec.size() > 0)
        pairParticles.upload(pairParticleVec);

    // Create kernels.
    
    map<string, string> defines;
    defines["NUM_ATOMS"] = cc.intToString(cc.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    defines["NUM_NORMAL_PARTICLES"] = cc.intToString(normalParticleVec.size());
    defines["NUM_PAIRS"] = cc.intToString(pairParticleVec.size());
    map<string, string> replacements;
    ComputeProgram program = cc.compileProgram(CommonDrudeKernelSources::drudeLangevin, defines);
    kernel1 = program->createKernel("integrateDrudeLangevinPart1");
    kernel2 = program->createKernel("integrateDrudeLangevinPart2");
    hardwallKernel = program->createKernel("applyHardWallConstraints");
    prevStepSize = -1.0;
}

void CommonIntegrateDrudeLangevinStepKernel::execute(ContextImpl& context, const DrudeLangevinIntegrator& integrator) {
    ContextSelector selector(cc);
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    int numAtoms = cc.getNumAtoms();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1->addArg(cc.getVelm());
        kernel1->addArg(cc.getLongForceBuffer());
        kernel1->addArg(integration.getPosDelta());
        kernel1->addArg(normalParticles);
        kernel1->addArg(pairParticles);
        kernel1->addArg(integration.getStepSize());
        for (int i = 0; i < 6; i++)
            kernel1->addArg();
        kernel1->addArg(integration.getRandom());
        kernel1->addArg();
        kernel2->addArg(cc.getPosq());
        if (cc.getUseMixedPrecision())
            kernel2->addArg(cc.getPosqCorrection());
        else
            kernel2->addArg(nullptr);
        kernel2->addArg(integration.getPosDelta());
        kernel2->addArg(cc.getVelm());
        kernel2->addArg(integration.getStepSize());
        hardwallKernel->addArg(cc.getPosq());
        if (cc.getUseMixedPrecision())
            hardwallKernel->addArg(cc.getPosqCorrection());
        else
            hardwallKernel->addArg(nullptr);
        hardwallKernel->addArg(cc.getVelm());
        hardwallKernel->addArg(pairParticles);
        hardwallKernel->addArg(integration.getStepSize());
        hardwallKernel->addArg();
        hardwallKernel->addArg();
    }
    
    // Compute integrator coefficients.
    
    double stepSize = integrator.getStepSize();
    double vscale = exp(-stepSize*integrator.getFriction());
    double fscale = (1-vscale)/integrator.getFriction()/(double) 0x100000000;
    double noisescale = sqrt(2*BOLTZ*integrator.getTemperature()*integrator.getFriction())*sqrt(0.5*(1-vscale*vscale)/integrator.getFriction());
    double vscaleDrude = exp(-stepSize*integrator.getDrudeFriction());
    double fscaleDrude = (1-vscaleDrude)/integrator.getDrudeFriction()/(double) 0x100000000;
    double noisescaleDrude = sqrt(2*BOLTZ*integrator.getDrudeTemperature()*integrator.getDrudeFriction())*sqrt(0.5*(1-vscaleDrude*vscaleDrude)/integrator.getDrudeFriction());
    double maxDrudeDistance = integrator.getMaxDrudeDistance();
    double hardwallscaleDrude = sqrt(BOLTZ*integrator.getDrudeTemperature());
    if (stepSize != prevStepSize) {
        if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
            mm_double2 ss = mm_double2(0, stepSize);
            integration.getStepSize().upload(&ss);
        }
        else {
            mm_float2 ss = mm_float2(0, (float) stepSize);
            integration.getStepSize().upload(&ss);
        }
        prevStepSize = stepSize;
    }
    if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
            kernel1->setArg(6, vscale);
            kernel1->setArg(7, fscale);
            kernel1->setArg(8, noisescale);
            kernel1->setArg(9, vscaleDrude);
            kernel1->setArg(10, fscaleDrude);
            kernel1->setArg(11, noisescaleDrude);
            hardwallKernel->setArg(5, maxDrudeDistance);
            hardwallKernel->setArg(6, hardwallscaleDrude);
    }
    else {
            kernel1->setArg(6, (float) vscale);
            kernel1->setArg(7, (float) fscale);
            kernel1->setArg(8, (float) noisescale);
            kernel1->setArg(9, (float) vscaleDrude);
            kernel1->setArg(10, (float) fscaleDrude);
            kernel1->setArg(11, (float) noisescaleDrude);
            hardwallKernel->setArg(5, (float) maxDrudeDistance);
            hardwallKernel->setArg(6, (float) hardwallscaleDrude);
    }

    // Call the first integration kernel.

    kernel1->setArg(13, integration.prepareRandomNumbers(normalParticles.getSize()+2*pairParticles.getSize()));
    kernel1->execute(numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    kernel2->execute(numAtoms);
    
    // Apply hard wall constraints.
    
    if (maxDrudeDistance > 0)
        hardwallKernel->execute(pairParticles.getSize());
    integration.computeVirtualSites();

    // Update the time and step count.

    cc.setTime(cc.getTime()+stepSize);
    cc.setStepCount(cc.getStepCount()+1);
    cc.reorderAtoms();
}

double CommonIntegrateDrudeLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const DrudeLangevinIntegrator& integrator) {
    return cc.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

CommonIntegrateDrudeSCFStepKernel::~CommonIntegrateDrudeSCFStepKernel() {
}

void CommonIntegrateDrudeSCFStepKernel::initialize(const System& system, const DrudeSCFIntegrator& integrator, const DrudeForce& force) {
    cc.initializeContexts();
    ContextSelector selector(cc);
    int numDrude = force.getNumParticles();
    drudeParams.initialize<mm_float4>(cc, numDrude, "drudeParams");
    drudeIndices.initialize<int>(cc, numDrude, "drudeIndices");
    drudeParents.initialize<mm_int4>(cc, numDrude, "drudeParents");
    vector<mm_float4> paramVec(numDrude);
    vector<mm_int4> parentVec(numDrude);
    drudeIndexVec.resize(numDrude);
    for (int i = 0; i < numDrude; i++) {
        int p, p1, p2, p3, p4;
        double charge, polarizability, aniso12, aniso34;
        force.getParticleParameters(i, p, p1, p2, p3, p4, charge, polarizability, aniso12, aniso34);
        double a1 = (p2 == -1 ? 1 : aniso12);
        double a2 = (p3 == -1 || p4 == -1 ? 1 : aniso34);
        double a3 = 3-a1-a2;
        double k3 = ONE_4PI_EPS0*charge*charge/(polarizability*a3);
        double k1 = ONE_4PI_EPS0*charge*charge/(polarizability*a1) - k3;
        double k2 = ONE_4PI_EPS0*charge*charge/(polarizability*a2) - k3;
        paramVec[i] = mm_float4((float) k1, (float) k2, (float) k3, 0.0f);
        drudeIndexVec[i] = p;
        parentVec[i] = mm_int4(p1, p2, p3, p4);
    }
    drudeParams.upload(paramVec);
    drudeIndices.upload(drudeIndexVec);
    drudeParents.upload(parentVec);

    // Create the kernels.
    
    ComputeProgram program = cc.compileProgram(CommonKernelSources::verlet);
    kernel1 = program->createKernel("integrateVerletPart1");
    kernel2 = program->createKernel("integrateVerletPart2");
    program = cc.compileProgram(CommonDrudeKernelSources::drudeSCF);
    minimizeKernel = program->createKernel("minimizeDrudePositions");
    prevStepSize = -1.0;
}

void CommonIntegrateDrudeSCFStepKernel::execute(ContextImpl& context, const DrudeSCFIntegrator& integrator) {
    ContextSelector selector(cc);
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    int numAtoms = cc.getNumAtoms();
    double dt = integrator.getStepSize();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1->addArg(numAtoms);
        kernel1->addArg(cc.getPaddedNumAtoms());
        kernel1->addArg(cc.getIntegrationUtilities().getStepSize());
        kernel1->addArg(cc.getPosq());
        kernel1->addArg(cc.getVelm());
        kernel1->addArg(cc.getLongForceBuffer());
        kernel1->addArg(integration.getPosDelta());
        if (cc.getUseMixedPrecision())
            kernel1->addArg(cc.getPosqCorrection());
        kernel2->addArg(numAtoms);
        kernel2->addArg(cc.getIntegrationUtilities().getStepSize());
        kernel2->addArg(cc.getPosq());
        kernel2->addArg(cc.getVelm());
        kernel2->addArg(integration.getPosDelta());
        if (cc.getUseMixedPrecision())
            kernel2->addArg(cc.getPosqCorrection());
        minimizeKernel->addArg((int) drudeParams.getSize());
        minimizeKernel->addArg(cc.getPaddedNumAtoms());
        minimizeKernel->addArg();
        minimizeKernel->addArg(cc.getPosq());
        minimizeKernel->addArg(cc.getLongForceBuffer());
        minimizeKernel->addArg(drudeParams);
        minimizeKernel->addArg(drudeIndices);
        minimizeKernel->addArg(drudeParents);
    }
    if (dt != prevStepSize) {
        if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
            vector<mm_double2> stepSizeVec(1);
            stepSizeVec[0] = mm_double2(dt, dt);
            cc.getIntegrationUtilities().getStepSize().upload(stepSizeVec);
        }
        else {
            vector<mm_float2> stepSizeVec(1);
            stepSizeVec[0] = mm_float2((float) dt, (float) dt);
            cc.getIntegrationUtilities().getStepSize().upload(stepSizeVec);
        }
        prevStepSize = dt;
    }

    // Call the first integration kernel.

    kernel1->execute(numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    kernel2->execute(numAtoms);

    // Update the positions of virtual sites and Drude particles.

    integration.computeVirtualSites();
    minimize(context, integrator.getMinimizationErrorTolerance());

    // Update the time and step count.

    cc.setTime(cc.getTime()+dt);
    cc.setStepCount(cc.getStepCount()+1);
    cc.reorderAtoms();
    
    // Reduce UI lag.
    
#ifdef WIN32
    cc.flushQueue();
#endif
}

double CommonIntegrateDrudeSCFStepKernel::computeKineticEnergy(ContextImpl& context, const DrudeSCFIntegrator& integrator) {
    return cc.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

void CommonIntegrateDrudeSCFStepKernel::minimize(ContextImpl& context, double tolerance) {
    minimizeKernel->setArg(2, (float) tolerance);
    long long* forces = (long long*) cc.getPinnedBuffer();
    double scale = 1/(double) 0x100000000;
    double lastForce = 0;
    int numDrude = drudeParams.getSize();
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    for (int iteration = 0; iteration < 50; iteration++) {
        context.calcForcesAndEnergy(true, false, context.getIntegrator().getIntegrationForceGroups());
        minimizeKernel->execute(drudeParams.getSize());
        cc.getLongForceBuffer().download(forces);
        double totalForce = 0;
        for (int i : drudeIndexVec) {
            Vec3 f(scale*forces[i], scale*forces[i+paddedNumAtoms], scale*forces[i+paddedNumAtoms*2]);
            totalForce += f.dot(f);
        }
        if (sqrt(totalForce/(3*numDrude)) < tolerance || (iteration > 0 && totalForce > 0.9*lastForce)) 
            break;
        lastForce = totalForce;
    }
}

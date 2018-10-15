/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2018 Stanford University and the Authors.      *
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

#include "CudaDrudeKernels.h"
#include "CudaDrudeKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "CudaBondedUtilities.h"
#include "CudaForceInfo.h"
#include "CudaIntegrationUtilities.h"
#include "CudaKernelSources.h"
#include "SimTKOpenMMRealType.h"
#include <set>

using namespace OpenMM;
using namespace std;

class CudaDrudeForceInfo : public CudaForceInfo {
public:
    CudaDrudeForceInfo(const DrudeForce& force) : force(force) {
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

void CudaCalcDrudeForceKernel::initialize(const System& system, const DrudeForce& force) {
    cu.setAsCurrent();
    if (cu.getContextIndex() != 0)
        return; // This is run entirely on one device
    int numParticles = force.getNumParticles();
    if (numParticles > 0) {
        // Create the harmonic interaction .
        
        vector<vector<int> > atoms(numParticles, vector<int>(5));
        particleParams.initialize<float4>(cu, numParticles, "drudeParticleParams");
        vector<float4> paramVector(numParticles);
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
                atoms[i][2] = 0;
                k1 = 0;
            }
            if (atoms[i][3] == -1 || atoms[i][4] == -1) {
                atoms[i][3] = 0;
                atoms[i][4] = 0;
                k2 = 0;
            }
            paramVector[i] = make_float4((float) k1, (float) k2, (float) k3, 0.0f);
        }
        particleParams.upload(paramVector);
        map<string, string> replacements;
        replacements["PARAMS"] = cu.getBondedUtilities().addArgument(particleParams.getDevicePointer(), "float4");
        cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaDrudeKernelSources::drudeParticleForce, replacements), force.getForceGroup());
    }
    int numPairs = force.getNumScreenedPairs();
    if (numPairs > 0) {
        // Create the screened interaction between dipole pairs.
        
        vector<vector<int> > atoms(numPairs, vector<int>(4));
        pairParams.initialize<float2>(cu, numPairs, "drudePairParams");
        vector<float2> paramVector(numPairs);
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
            paramVector[i] = make_float2((float) screeningScale, (float) energyScale);
        }
        pairParams.upload(paramVector);
        map<string, string> replacements;
        replacements["PARAMS"] = cu.getBondedUtilities().addArgument(pairParams.getDevicePointer(), "float2");
        cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaDrudeKernelSources::drudePairForce, replacements), force.getForceGroup());
    }
    cu.addForce(new CudaDrudeForceInfo(force));
}

double CudaCalcDrudeForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CudaCalcDrudeForceKernel::copyParametersToContext(ContextImpl& context, const DrudeForce& force) {
    if (cu.getContextIndex() != 0)
        return; // This is run entirely on one device
    
    // Set the particle parameters.
    
    int numParticles = force.getNumParticles();
    if (numParticles > 0) {
        if (!particleParams.isInitialized() || numParticles != particleParams.getSize())
            throw OpenMMException("updateParametersInContext: The number of Drude particles has changed");
        vector<float4> paramVector(numParticles);
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
            paramVector[i] = make_float4((float) k1, (float) k2, (float) k3, 0.0f);
        }
        particleParams.upload(paramVector);
    }
    
    // Set the pair parameters.
    
    int numPairs = force.getNumScreenedPairs();
    if (numPairs > 0) {
        if (!pairParams.isInitialized() || numPairs != pairParams.getSize())
            throw OpenMMException("updateParametersInContext: The number of screened pairs has changed");
        vector<float2> paramVector(numPairs);
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
            paramVector[i] = make_float2((float) screeningScale, (float) energyScale);
        }
        pairParams.upload(paramVector);
    }
}

void CudaIntegrateDrudeLangevinStepKernel::initialize(const System& system, const DrudeLangevinIntegrator& integrator, const DrudeForce& force) {
    cu.getPlatformData().initializeContexts(system);
    cu.getIntegrationUtilities().initRandomNumberGenerator((unsigned int) integrator.getRandomNumberSeed());
    
    // Identify particle pairs and ordinary particles.
    
    set<int> particles;
    vector<int> normalParticleVec;
    vector<int2> pairParticleVec;
    for (int i = 0; i < system.getNumParticles(); i++)
        particles.insert(i);
    for (int i = 0; i < force.getNumParticles(); i++) {
        int p, p1, p2, p3, p4;
        double charge, polarizability, aniso12, aniso34;
        force.getParticleParameters(i, p, p1, p2, p3, p4, charge, polarizability, aniso12, aniso34);
        particles.erase(p);
        particles.erase(p1);
        pairParticleVec.push_back(make_int2(p, p1));
    }
    normalParticleVec.insert(normalParticleVec.begin(), particles.begin(), particles.end());
    normalParticles.initialize<int>(cu, max((int) normalParticleVec.size(), 1), "drudeNormalParticles");
    pairParticles.initialize<int2>(cu, max((int) pairParticleVec.size(), 1), "drudePairParticles");
    if (normalParticleVec.size() > 0)
        normalParticles.upload(normalParticleVec);
    if (pairParticleVec.size() > 0)
        pairParticles.upload(pairParticleVec);

    // Create kernels.
    
    map<string, string> defines;
    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    defines["NUM_NORMAL_PARTICLES"] = cu.intToString(normalParticleVec.size());
    defines["NUM_PAIRS"] = cu.intToString(pairParticleVec.size());
    map<string, string> replacements;
    CUmodule module = cu.createModule(CudaKernelSources::vectorOps+CudaDrudeKernelSources::drudeLangevin, defines, "");
    kernel1 = cu.getKernel(module, "integrateDrudeLangevinPart1");
    kernel2 = cu.getKernel(module, "integrateDrudeLangevinPart2");
    hardwallKernel = cu.getKernel(module, "applyHardWallConstraints");
    prevStepSize = -1.0;
}

void CudaIntegrateDrudeLangevinStepKernel::execute(ContextImpl& context, const DrudeLangevinIntegrator& integrator) {
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    
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
        if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
            double2 ss = make_double2(0, stepSize);
            integration.getStepSize().upload(&ss);
        }
        else {
            float2 ss = make_float2(0, (float) stepSize);
            integration.getStepSize().upload(&ss);
        }
        prevStepSize = stepSize;
    }
    
    // Create appropriate pointer for the precision mode.
    
    float vscaleFloat = (float) vscale;
    float fscaleFloat = (float) fscale;
    float noisescaleFloat = (float) noisescale;
    float vscaleDrudeFloat = (float) vscaleDrude;
    float fscaleDrudeFloat = (float) fscaleDrude;
    float noisescaleDrudeFloat = (float) noisescaleDrude;
    float maxDrudeDistanceFloat =(float) maxDrudeDistance;
    float hardwallscaleDrudeFloat = (float) hardwallscaleDrude;
    void *vscalePtr, *fscalePtr, *noisescalePtr, *vscaleDrudePtr, *fscaleDrudePtr, *noisescaleDrudePtr, *maxDrudeDistancePtr, *hardwallscaleDrudePtr;
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
        vscalePtr = &vscale;
        fscalePtr = &fscale;
        noisescalePtr = &noisescale;
        vscaleDrudePtr = &vscaleDrude;
        fscaleDrudePtr = &fscaleDrude;
        noisescaleDrudePtr = &noisescaleDrude;
        maxDrudeDistancePtr = &maxDrudeDistance;
        hardwallscaleDrudePtr = &hardwallscaleDrude;
    }
    else {
        vscalePtr = &vscaleFloat;
        fscalePtr = &fscaleFloat;
        noisescalePtr = &noisescaleFloat;
        vscaleDrudePtr = &vscaleDrudeFloat;
        fscaleDrudePtr = &fscaleDrudeFloat;
        noisescaleDrudePtr = &noisescaleDrudeFloat;
        maxDrudeDistancePtr = &maxDrudeDistanceFloat;
        hardwallscaleDrudePtr = &hardwallscaleDrudeFloat;
    }

    // Call the first integration kernel.

    int randomIndex = integration.prepareRandomNumbers(normalParticles.getSize()+2*pairParticles.getSize());
    void* args1[] = {&cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &integration.getPosDelta().getDevicePointer(),
            &normalParticles.getDevicePointer(), &pairParticles.getDevicePointer(), &integration.getStepSize().getDevicePointer(),
            vscalePtr, fscalePtr, noisescalePtr, vscaleDrudePtr, fscaleDrudePtr, noisescaleDrudePtr, &integration.getRandom().getDevicePointer(), &randomIndex};
    cu.executeKernel(kernel1, args1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args2[] = {&cu.getPosq().getDevicePointer(), &posCorrection, &integration.getPosDelta().getDevicePointer(),
            &cu.getVelm().getDevicePointer(), &integration.getStepSize().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms);
    
    // Apply hard wall constraints.
    
    if (maxDrudeDistance > 0) {
        void* hardwallArgs[] = {&cu.getPosq().getDevicePointer(), &posCorrection, &cu.getVelm().getDevicePointer(),
                &pairParticles.getDevicePointer(), &integration.getStepSize().getDevicePointer(), maxDrudeDistancePtr, hardwallscaleDrudePtr};
        cu.executeKernel(hardwallKernel, hardwallArgs, pairParticles.getSize());
    }
    integration.computeVirtualSites();

    // Update the time and step count.

    cu.setTime(cu.getTime()+stepSize);
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
}

double CudaIntegrateDrudeLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const DrudeLangevinIntegrator& integrator) {
    return cu.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

CudaIntegrateDrudeSCFStepKernel::~CudaIntegrateDrudeSCFStepKernel() {
    if (minimizerPos != NULL)
        lbfgs_free(minimizerPos);
}

void CudaIntegrateDrudeSCFStepKernel::initialize(const System& system, const DrudeSCFIntegrator& integrator, const DrudeForce& force) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();

    // Identify Drude particles.
    
    for (int i = 0; i < force.getNumParticles(); i++) {
        int p, p1, p2, p3, p4;
        double charge, polarizability, aniso12, aniso34;
        force.getParticleParameters(i, p, p1, p2, p3, p4, charge, polarizability, aniso12, aniso34);
        drudeParticles.push_back(p);
    }
    
    // Initialize the energy minimizer.
    
    minimizerPos = lbfgs_malloc(drudeParticles.size()*3);
    if (minimizerPos == NULL)
        throw OpenMMException("DrudeSCFIntegrator: Failed to allocate memory");
    lbfgs_parameter_init(&minimizerParams);
    minimizerParams.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;    

    // Create the kernels.
    
    map<string, string> defines;
    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    CUmodule module = cu.createModule(CudaKernelSources::verlet, defines, "");
    kernel1 = cu.getKernel(module, "integrateVerletPart1");
    kernel2 = cu.getKernel(module, "integrateVerletPart2");
    prevStepSize = -1.0;
}

void CudaIntegrateDrudeSCFStepKernel::execute(ContextImpl& context, const DrudeSCFIntegrator& integrator) {
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    double dt = integrator.getStepSize();
    if (dt != prevStepSize) {
        if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
            vector<double2> stepSizeVec(1);
            stepSizeVec[0] = make_double2(dt, dt);
            cu.getIntegrationUtilities().getStepSize().upload(stepSizeVec);
        }
        else {
            vector<float2> stepSizeVec(1);
            stepSizeVec[0] = make_float2((float) dt, (float) dt);
            cu.getIntegrationUtilities().getStepSize().upload(stepSizeVec);
        }
        prevStepSize = dt;
    }

    // Call the first integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args1[] = {&numAtoms, &paddedNumAtoms, &cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(), &posCorrection,
            &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel1, args1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    void* args2[] = {&numAtoms, &cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(), &posCorrection,
            &cu.getVelm().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms);

    // Update the positions of virtual sites and Drude particles.

    integration.computeVirtualSites();
    minimize(context, integrator.getMinimizationErrorTolerance());

    // Update the time and step count.

    cu.setTime(cu.getTime()+dt);
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
}

double CudaIntegrateDrudeSCFStepKernel::computeKineticEnergy(ContextImpl& context, const DrudeSCFIntegrator& integrator) {
    return cu.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

struct MinimizerData {
    ContextImpl& context;
    CudaContext& cu;
    vector<int>& drudeParticles;
    MinimizerData(ContextImpl& context, CudaContext& cu, vector<int>& drudeParticles) : context(context), cu(cu), drudeParticles(drudeParticles) {}
};

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    MinimizerData* data = reinterpret_cast<MinimizerData*>(instance);
    ContextImpl& context = data->context;
    CudaContext& cu = data->cu;
    vector<int>& drudeParticles = data->drudeParticles;
    int numDrudeParticles = drudeParticles.size();

    // Set the particle positions.
    
    cu.getPosq().download(cu.getPinnedBuffer());
    if (cu.getUseDoublePrecision()) {
        double4* posq = (double4*) cu.getPinnedBuffer();
        for (int i = 0; i < numDrudeParticles; ++i) {
            double4& p = posq[drudeParticles[i]];
            p.x = x[3*i];
            p.y = x[3*i+1];
            p.z = x[3*i+2];
        }
    }
    else {
        float4* posq = (float4*) cu.getPinnedBuffer();
        for (int i = 0; i < numDrudeParticles; ++i) {
            float4& p = posq[drudeParticles[i]];
            p.x = x[3*i];
            p.y = x[3*i+1];
            p.z = x[3*i+2];
        }
    }
    cu.getPosq().upload(cu.getPinnedBuffer());

    // Compute the forces and energy for this configuration.

    double energy = context.calcForcesAndEnergy(true, true);
    long long* force = (long long*) cu.getPinnedBuffer();
    cu.getForce().download(force);
    double forceScale = -1.0/0x100000000;
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    for (int i = 0; i < numDrudeParticles; ++i) {
        int index = drudeParticles[i];
        g[3*i] = forceScale*force[index];
        g[3*i+1] = forceScale*force[index+paddedNumAtoms];
        g[3*i+2] = forceScale*force[index+paddedNumAtoms*2];
    }
    return energy;
}

void CudaIntegrateDrudeSCFStepKernel::minimize(ContextImpl& context, double tolerance) {
    // Record the initial positions.

    int numDrudeParticles = drudeParticles.size();
    cu.getPosq().download(cu.getPinnedBuffer());
    if (cu.getUseDoublePrecision()) {
        double4* posq = (double4*) cu.getPinnedBuffer();
        for (int i = 0; i < numDrudeParticles; ++i) {
            double4 p = posq[drudeParticles[i]];
            minimizerPos[3*i] = p.x;
            minimizerPos[3*i+1] = p.y;
            minimizerPos[3*i+2] = p.z;
        }
    }
    else {
        float4* posq = (float4*) cu.getPinnedBuffer();
        for (int i = 0; i < numDrudeParticles; ++i) {
            float4 p = posq[drudeParticles[i]];
            minimizerPos[3*i] = p.x;
            minimizerPos[3*i+1] = p.y;
            minimizerPos[3*i+2] = p.z;
        }
        minimizerParams.xtol = 1e-7;
    }
    
    // Determine a normalization constant for scaling the tolerance.
    
    double norm = 0.0;
    for (int i = 0; i < 3*numDrudeParticles; i++)
        norm += minimizerPos[i]*minimizerPos[i];
    norm /= numDrudeParticles;
    norm = (norm < 1 ? 1 : sqrt(norm));
    minimizerParams.epsilon = tolerance/norm;
    
    // Perform the minimization.

    lbfgsfloatval_t fx;
    MinimizerData data(context, cu, drudeParticles);
    lbfgs(numDrudeParticles*3, minimizerPos, &fx, evaluate, NULL, &data, &minimizerParams);
}
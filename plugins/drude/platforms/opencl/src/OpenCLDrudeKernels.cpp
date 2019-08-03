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

#include "OpenCLDrudeKernels.h"
#include "OpenCLDrudeKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "OpenCLBondedUtilities.h"
#include "OpenCLForceInfo.h"
#include "OpenCLIntegrationUtilities.h"
#include "OpenCLKernelSources.h"
#include "SimTKOpenMMRealType.h"
#include <set>

using namespace OpenMM;
using namespace std;

static void setPosqCorrectionArg(OpenCLContext& cl, cl::Kernel& kernel, int index) {
    if (cl.getUseMixedPrecision())
        kernel.setArg<cl::Buffer>(index, cl.getPosqCorrection().getDeviceBuffer());
    else
        kernel.setArg<void*>(index, NULL);
}

class OpenCLDrudeForceInfo : public OpenCLForceInfo {
public:
    OpenCLDrudeForceInfo(const DrudeForce& force) : OpenCLForceInfo(0), force(force) {
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

void OpenCLCalcDrudeForceKernel::initialize(const System& system, const DrudeForce& force) {
    if (cl.getContextIndex() != 0)
        return; // This is run entirely on one device
    int numParticles = force.getNumParticles();
    if (numParticles > 0) {
        // Create the harmonic interaction .
        
        vector<vector<int> > atoms(numParticles, vector<int>(5));
        particleParams.initialize<mm_float4>(cl, numParticles, "drudeParticleParams");
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
                atoms[i][2] = 0;
                k1 = 0;
            }
            if (atoms[i][3] == -1 || atoms[i][4] == -1) {
                atoms[i][3] = 0;
                atoms[i][4] = 0;
                k2 = 0;
            }
            paramVector[i] = mm_float4((float) k1, (float) k2, (float) k3, 0.0f);
        }
        particleParams.upload(paramVector);
        map<string, string> replacements;
        replacements["PARAMS"] = cl.getBondedUtilities().addArgument(particleParams.getDeviceBuffer(), "float4");
        cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLDrudeKernelSources::drudeParticleForce, replacements), force.getForceGroup());
    }
    int numPairs = force.getNumScreenedPairs();
    if (numPairs > 0) {
        // Create the screened interaction between dipole pairs.
        
        vector<vector<int> > atoms(numPairs, vector<int>(4));
        pairParams.initialize<mm_float2>(cl, numPairs, "drudePairParams");
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
        replacements["PARAMS"] = cl.getBondedUtilities().addArgument(pairParams.getDeviceBuffer(), "float2");
        cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLDrudeKernelSources::drudePairForce, replacements), force.getForceGroup());
    }
    cl.addForce(new OpenCLDrudeForceInfo(force));
}

double OpenCLCalcDrudeForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void OpenCLCalcDrudeForceKernel::copyParametersToContext(ContextImpl& context, const DrudeForce& force) {
    if (cl.getContextIndex() != 0)
        return; // This is run entirely on one device
    
    // Set the particle parameters.
    
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

void OpenCLIntegrateDrudeLangevinStepKernel::initialize(const System& system, const DrudeLangevinIntegrator& integrator, const DrudeForce& force) {
    cl.getPlatformData().initializeContexts(system);
    cl.getIntegrationUtilities().initRandomNumberGenerator((unsigned int) integrator.getRandomNumberSeed());
    
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
    normalParticles.initialize<int>(cl, max((int) normalParticleVec.size(), 1), "drudeNormalParticles");
    pairParticles.initialize<cl_int2>(cl, max((int) pairParticleVec.size(), 1), "drudePairParticles");
    if (normalParticleVec.size() > 0)
        normalParticles.upload(normalParticleVec);
    if (pairParticleVec.size() > 0)
        pairParticles.upload(pairParticleVec);

    // Create kernels.
    
    map<string, string> defines;
    defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
    defines["NUM_NORMAL_PARTICLES"] = cl.intToString(normalParticleVec.size());
    defines["NUM_PAIRS"] = cl.intToString(pairParticleVec.size());
    map<string, string> replacements;
    cl::Program program = cl.createProgram(OpenCLDrudeKernelSources::drudeLangevin, defines, "");
    kernel1 = cl::Kernel(program, "integrateDrudeLangevinPart1");
    kernel2 = cl::Kernel(program, "integrateDrudeLangevinPart2");
    hardwallKernel = cl::Kernel(program, "applyHardWallConstraints");
    prevStepSize = -1.0;
}

void OpenCLIntegrateDrudeLangevinStepKernel::execute(ContextImpl& context, const DrudeLangevinIntegrator& integrator) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int numAtoms = cl.getNumAtoms();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1.setArg<cl::Buffer>(0, cl.getVelm().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(1, cl.getForce().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(3, normalParticles.getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(4, pairParticles.getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(5, integration.getStepSize().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(12, integration.getRandom().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
        if (cl.getUseMixedPrecision())
            kernel2.setArg<cl::Buffer>(1, cl.getPosqCorrection().getDeviceBuffer());
        else
            kernel2.setArg<void*>(1, NULL);
        kernel2.setArg<cl::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(3, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(4, integration.getStepSize().getDeviceBuffer());
        hardwallKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
        if (cl.getUseMixedPrecision())
            hardwallKernel.setArg<cl::Buffer>(1, cl.getPosqCorrection().getDeviceBuffer());
        else
            hardwallKernel.setArg<void*>(1, NULL);
        hardwallKernel.setArg<cl::Buffer>(2, cl.getVelm().getDeviceBuffer());
        hardwallKernel.setArg<cl::Buffer>(3, pairParticles.getDeviceBuffer());
        hardwallKernel.setArg<cl::Buffer>(4, integration.getStepSize().getDeviceBuffer());
    }
    
    // Compute integrator coefficients.
    
    double stepSize = integrator.getStepSize();
    double vscale = exp(-stepSize*integrator.getFriction());
    double fscale = (1-vscale)/integrator.getFriction();
    double noisescale = sqrt(2*BOLTZ*integrator.getTemperature()*integrator.getFriction())*sqrt(0.5*(1-vscale*vscale)/integrator.getFriction());
    double vscaleDrude = exp(-stepSize*integrator.getDrudeFriction());
    double fscaleDrude = (1-vscaleDrude)/integrator.getDrudeFriction();
    double noisescaleDrude = sqrt(2*BOLTZ*integrator.getDrudeTemperature()*integrator.getDrudeFriction())*sqrt(0.5*(1-vscaleDrude*vscaleDrude)/integrator.getDrudeFriction());
    double maxDrudeDistance = integrator.getMaxDrudeDistance();
    double hardwallscaleDrude = sqrt(BOLTZ*integrator.getDrudeTemperature());
    if (stepSize != prevStepSize) {
        if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
            mm_double2 ss = mm_double2(0, stepSize);
            integration.getStepSize().upload(&ss);
        }
        else {
            mm_float2 ss = mm_float2(0, (float) stepSize);
            integration.getStepSize().upload(&ss);
        }
        prevStepSize = stepSize;
    }
    if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
            kernel1.setArg<cl_double>(6, vscale);
            kernel1.setArg<cl_double>(7, fscale);
            kernel1.setArg<cl_double>(8, noisescale);
            kernel1.setArg<cl_double>(9, vscaleDrude);
            kernel1.setArg<cl_double>(10, fscaleDrude);
            kernel1.setArg<cl_double>(11, noisescaleDrude);
            hardwallKernel.setArg<cl_double>(5, maxDrudeDistance);
            hardwallKernel.setArg<cl_double>(6, hardwallscaleDrude);
    }
    else {
            kernel1.setArg<cl_float>(6, (cl_float) vscale);
            kernel1.setArg<cl_float>(7, (cl_float) fscale);
            kernel1.setArg<cl_float>(8, (cl_float) noisescale);
            kernel1.setArg<cl_float>(9, (cl_float) vscaleDrude);
            kernel1.setArg<cl_float>(10, (cl_float) fscaleDrude);
            kernel1.setArg<cl_float>(11, (cl_float) noisescaleDrude);
            hardwallKernel.setArg<cl_float>(5, (cl_float) maxDrudeDistance);
            hardwallKernel.setArg<cl_float>(6, (cl_float) hardwallscaleDrude);
    }

    // Call the first integration kernel.

    kernel1.setArg<cl_uint>(13, integration.prepareRandomNumbers(normalParticles.getSize()+2*pairParticles.getSize()));
    cl.executeKernel(kernel1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    cl.executeKernel(kernel2, numAtoms);
    
    // Apply hard wall constraints.
    
    if (maxDrudeDistance > 0)
        cl.executeKernel(hardwallKernel, pairParticles.getSize());
    integration.computeVirtualSites();

    // Update the time and step count.

    cl.setTime(cl.getTime()+stepSize);
    cl.setStepCount(cl.getStepCount()+1);
    cl.reorderAtoms();
}

double OpenCLIntegrateDrudeLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const DrudeLangevinIntegrator& integrator) {
    return cl.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

OpenCLIntegrateDrudeSCFStepKernel::~OpenCLIntegrateDrudeSCFStepKernel() {
    if (minimizerPos != NULL)
        lbfgs_free(minimizerPos);
}

void OpenCLIntegrateDrudeSCFStepKernel::initialize(const System& system, const DrudeSCFIntegrator& integrator, const DrudeForce& force) {
    cl.getPlatformData().initializeContexts(system);

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
    
    cl::Program program = cl.createProgram(OpenCLKernelSources::verlet, "");
    kernel1 = cl::Kernel(program, "integrateVerletPart1");
    kernel2 = cl::Kernel(program, "integrateVerletPart2");
    prevStepSize = -1.0;
}

void OpenCLIntegrateDrudeSCFStepKernel::execute(ContextImpl& context, const DrudeSCFIntegrator& integrator) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int numAtoms = cl.getNumAtoms();
    double dt = integrator.getStepSize();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1.setArg<cl_int>(0, numAtoms);
        kernel1.setArg<cl::Buffer>(1, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
        setPosqCorrectionArg(cl, kernel1, 3);
        kernel1.setArg<cl::Buffer>(4, cl.getVelm().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(5, cl.getForce().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(6, integration.getPosDelta().getDeviceBuffer());
        kernel2.setArg<cl_int>(0, numAtoms);
        kernel2.setArg<cl::Buffer>(1, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
        setPosqCorrectionArg(cl, kernel2, 3);
        kernel2.setArg<cl::Buffer>(4, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(5, integration.getPosDelta().getDeviceBuffer());
    }
    if (dt != prevStepSize) {
        if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
            vector<mm_double2> stepSizeVec(1);
            stepSizeVec[0] = mm_double2(dt, dt);
            cl.getIntegrationUtilities().getStepSize().upload(stepSizeVec);
        }
        else {
            vector<mm_float2> stepSizeVec(1);
            stepSizeVec[0] = mm_float2((cl_float) dt, (cl_float) dt);
            cl.getIntegrationUtilities().getStepSize().upload(stepSizeVec);
        }
        prevStepSize = dt;
    }

    // Call the first integration kernel.

    cl.executeKernel(kernel1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    cl.executeKernel(kernel2, numAtoms);

    // Update the positions of virtual sites and Drude particles.

    integration.computeVirtualSites();
    minimize(context, integrator.getMinimizationErrorTolerance());

    // Update the time and step count.

    cl.setTime(cl.getTime()+dt);
    cl.setStepCount(cl.getStepCount()+1);
    cl.reorderAtoms();
    
    // Reduce UI lag.
    
#ifdef WIN32
    cl.getQueue().flush();
#endif
}

double OpenCLIntegrateDrudeSCFStepKernel::computeKineticEnergy(ContextImpl& context, const DrudeSCFIntegrator& integrator) {
    return cl.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

struct MinimizerData {
    ContextImpl& context;
    OpenCLContext& cl;
    vector<int>& drudeParticles;
    MinimizerData(ContextImpl& context, OpenCLContext& cl, vector<int>& drudeParticles) : context(context), cl(cl), drudeParticles(drudeParticles) {}
};

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    MinimizerData* data = reinterpret_cast<MinimizerData*>(instance);
    ContextImpl& context = data->context;
    OpenCLContext& cl = data->cl;
    vector<int>& drudeParticles = data->drudeParticles;
    int numDrudeParticles = drudeParticles.size();

    // Set the particle positions.
    
    cl.getPosq().download(cl.getPinnedBuffer());
    if (cl.getUseDoublePrecision()) {
        mm_double4* posq = (mm_double4*) cl.getPinnedBuffer();
        for (int i = 0; i < numDrudeParticles; ++i) {
            mm_double4& p = posq[drudeParticles[i]];
            p.x = x[3*i];
            p.y = x[3*i+1];
            p.z = x[3*i+2];
        }
    }
    else {
        mm_float4* posq = (mm_float4*) cl.getPinnedBuffer();
        for (int i = 0; i < numDrudeParticles; ++i) {
            mm_float4& p = posq[drudeParticles[i]];
            p.x = x[3*i];
            p.y = x[3*i+1];
            p.z = x[3*i+2];
        }
    }
    cl.getPosq().upload(cl.getPinnedBuffer());

    // Compute the forces and energy for this configuration.

    double energy = context.calcForcesAndEnergy(true, true);
    cl.getForce().download(cl.getPinnedBuffer());
    if (cl.getUseDoublePrecision()) {
        mm_double4* force = (mm_double4*) cl.getPinnedBuffer();
        for (int i = 0; i < numDrudeParticles; ++i) {
            int index = drudeParticles[i];
            g[3*i] = -force[index].x;
            g[3*i+1] = -force[index].y;
            g[3*i+2] = -force[index].z;
        }
    }
    else {
        mm_float4* force = (mm_float4*) cl.getPinnedBuffer();
        for (int i = 0; i < numDrudeParticles; ++i) {
            int index = drudeParticles[i];
            g[3*i] = -force[index].x;
            g[3*i+1] = -force[index].y;
            g[3*i+2] = -force[index].z;
        }
    }
    return energy;
}

void OpenCLIntegrateDrudeSCFStepKernel::minimize(ContextImpl& context, double tolerance) {
    // Record the initial positions.

    int numDrudeParticles = drudeParticles.size();
    cl.getPosq().download(cl.getPinnedBuffer());
    if (cl.getUseDoublePrecision()) {
        mm_double4* posq = (mm_double4*) cl.getPinnedBuffer();
        for (int i = 0; i < numDrudeParticles; ++i) {
            mm_double4 p = posq[drudeParticles[i]];
            minimizerPos[3*i] = p.x;
            minimizerPos[3*i+1] = p.y;
            minimizerPos[3*i+2] = p.z;
        }
    }
    else {
        mm_float4* posq = (mm_float4*) cl.getPinnedBuffer();
        for (int i = 0; i < numDrudeParticles; ++i) {
            mm_float4 p = posq[drudeParticles[i]];
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
    MinimizerData data(context, cl, drudeParticles);
    lbfgs(numDrudeParticles*3, minimizerPos, &fx, evaluate, NULL, &data, &minimizerParams);
}
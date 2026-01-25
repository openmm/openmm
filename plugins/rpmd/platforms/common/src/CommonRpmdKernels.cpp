/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2011-2021 Stanford University and the Authors.      *
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

#include "CommonRpmdKernels.h"
#include "CommonRpmdKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/PythonForceImpl.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/common/IntegrationUtilities.h"
#include "openmm/common/ExpressionUtilities.h"
#include "openmm/common/NonbondedUtilities.h"
#include "openmm/PythonForce.h"
#include "openmm/CMMotionRemover.h"
#include "SimTKOpenMMRealType.h"
#include <fstream>
#include <chrono>
#include <random>
#include <sstream>

using namespace OpenMM;
using namespace std;

static void appendDebugLog(const char* location, const char* message, const std::string& data, const char* hypothesisId, const char* runId) {
    std::ofstream out("/media/extradrive/Trajectories/openmm/.cursor/debug.log", std::ios::app);
    if (!out)
        return;
    const long long timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()).count();
    out << "{\"sessionId\":\"debug-session\",\"runId\":\"" << runId
        << "\",\"hypothesisId\":\"" << hypothesisId
        << "\",\"location\":\"" << location
        << "\",\"message\":\"" << message
        << "\",\"data\":" << data
        << ",\"timestamp\":" << timestamp << "}\n";
}


/**
 * Select a size for an FFT that is a multiple of 2, 3, 5, and 7.
 */
static int findFFTDimension(int minimum) {
    if (minimum < 1)
        return 1;
    while (true) {
        // Attempt to factor the current value.

        int unfactored = minimum;
        for (int factor = 2; factor < 8; factor++) {
            while (unfactored > 1 && unfactored%factor == 0)
                unfactored /= factor;
        }
        if (unfactored == 1)
            return minimum;
        minimum++;
    }
}

void CommonIntegrateRPMDStepKernel::initialize(const System& system, const RPMDIntegrator& integrator) {
    cc.initializeContexts();
    ContextSelector selector(cc);
    numCopies = integrator.getNumCopies();
    numParticles = system.getNumParticles();
    workgroupSize = numCopies;
    if (numCopies != findFFTDimension(numCopies))
        throw OpenMMException("RPMDIntegrator: the number of copies must be a multiple of powers of 2, 3, and 5.");
    int paddedParticles = cc.getPaddedNumAtoms();
    bool useDoublePrecision = (cc.getUseDoublePrecision() || cc.getUseMixedPrecision());
    int elementSize = (useDoublePrecision ? sizeof(mm_double4) : sizeof(mm_float4));
    
    // Classify particles into quantum and classical based on particle types
    const map<int, int>& particleTypes = integrator.getParticleTypes();
    const set<int>& quantumTypes = integrator.getQuantumParticleTypes();
    bool defaultQuantum = integrator.getDefaultQuantum();
    
    // Check if hybrid mode is actually being used
    bool hybridModeRequested = (!particleTypes.empty() || !quantumTypes.empty() || !defaultQuantum);
    
    vector<int> quantumParticlesList;
    vector<int> classicalParticlesList;
    
    if (hybridModeRequested) {
        for (int i = 0; i < numParticles; i++) {
            if (system.getParticleMass(i) == 0)
                continue;  // Skip massless particles
            
            // Determine particle type (default to 0 if not set)
            int type = 0;
            auto it = particleTypes.find(i);
            if (it != particleTypes.end())
                type = it->second;
            
            // Check if this type is quantum
            bool isQuantum = (type == 0) ? defaultQuantum : (quantumTypes.count(type) > 0);
            
            if (isQuantum)
                quantumParticlesList.push_back(i);
            else
                classicalParticlesList.push_back(i);
        }
    } else {
        // Default: all particles are quantum (standard RPMD)
        for (int i = 0; i < numParticles; i++) {
            if (system.getParticleMass(i) > 0)
                quantumParticlesList.push_back(i);
        }
    }
    
    numQuantumParticles = quantumParticlesList.size();
    numClassicalParticles = classicalParticlesList.size();
    hybridMode = (numClassicalParticles > 0);
    
    // Store particle indices for later use
    quantumParticleIndices = quantumParticlesList;
    classicalParticleIndices = classicalParticlesList;
    
    // UNIFORM MEMORY LAYOUT: All particles stored with all beads
    // Classical particles have all beads kept identical via sync kernel
    // This is simpler and more robust than sparse storage
    forces.initialize<long long>(cc, numCopies * paddedParticles * 3, "rpmdForces");
    positions.initialize(cc, numCopies * paddedParticles, elementSize, "rpmdPositions");
    velocities.initialize(cc, numCopies * paddedParticles, elementSize, "rpmdVelocities");
    
    // Create per-particle isQuantum flag array for GPU
    vector<int> isQuantumFlags(paddedParticles, 1);  // Default: quantum
    if (hybridMode) {
        // Mark classical particles
        for (int idx : classicalParticlesList) {
            isQuantumFlags[idx] = 0;
        }
        // Also mark massless particles as classical (no dynamics)
        for (int i = 0; i < numParticles; i++) {
            if (system.getParticleMass(i) == 0)
                isQuantumFlags[i] = 0;
        }
        
        // Allocate buffer for classical KE computation (for Bussi thermostat)
        int numClassicalWorkGroups = (numClassicalParticles + workgroupSize - 1) / workgroupSize;
        classicalKE.initialize<double>(cc, std::max(1, numClassicalWorkGroups), "classicalKE");
    }
    isQuantum.initialize<int>(cc, paddedParticles, "isQuantum");
    isQuantum.upload(isQuantumFlags);

    // #region agent log
    {
        std::ostringstream data;
        data << "{\"numCopies\":" << numCopies
             << ",\"numParticles\":" << numParticles
             << ",\"numQuantumParticles\":" << numQuantumParticles
             << ",\"numClassicalParticles\":" << numClassicalParticles
             << ",\"hybridMode\":" << (hybridMode ? 1 : 0)
             << ",\"storageLayout\":\"uniform\""
             << "}";
        appendDebugLog("CommonRpmdKernels.cpp:initialize", "uniform storage", data.str(), "uniform-v1", "run1");
    }
    // #endregion
    
    cc.getIntegrationUtilities().initRandomNumberGenerator((unsigned int) integrator.getRandomNumberSeed());
    
    // Fill in the posq and velm arrays with safe values to avoid a risk of nans.
    
    if (useDoublePrecision) {
        vector<mm_double4> temp(positions.getSize());
        for (int i = 0; i < positions.getSize(); i++)
            temp[i] = mm_double4(0, 0, 0, 0);
        positions.upload(temp);
        for (int i = 0; i < velocities.getSize(); i++)
            temp[i] = mm_double4(0, 0, 0, 1);
        velocities.upload(temp);
    }
    else {
        vector<mm_float4> temp(positions.getSize());
        for (int i = 0; i < positions.getSize(); i++)
            temp[i] = mm_float4(0, 0, 0, 0);
        positions.upload(temp);
        for (int i = 0; i < velocities.getSize(); i++)
            temp[i] = mm_float4(0, 0, 0, 1);
        velocities.upload(temp);
    }
    
    // Build a list of contractions.
    
    groupsNotContracted = -1;
    const map<int, int>& contractions = integrator.getContractions();
    int maxContractedCopies = 0;
    for (auto& c : contractions) {
        int group = c.first;
        int copies = c.second;
        if (group < 0 || group > 31)
            throw OpenMMException("RPMDIntegrator: Force group must be between 0 and 31");
        if (copies < 0 || copies > numCopies)
            throw OpenMMException("RPMDIntegrator: Number of copies for contraction cannot be greater than the total number of copies being simulated");
        if (copies != findFFTDimension(copies))
            throw OpenMMException("RPMDIntegrator: Number of copies for contraction must be a multiple of powers of 2, 3, and 5.");
        if (copies != numCopies) {
            if (groupsByCopies.find(copies) == groupsByCopies.end()) {
                groupsByCopies[copies] = 1<<group;
                if (copies > maxContractedCopies)
                    maxContractedCopies = copies;
            }
            else
                groupsByCopies[copies] |= 1<<group;
            groupsNotContracted -= 1<<group;
        }
    }
    groupsNotContracted &= integrator.getIntegrationForceGroups();
    if (maxContractedCopies > 0) {
        contractedForces.initialize<long long>(cc, maxContractedCopies*paddedParticles*3, "rpmdContractedForces");
        contractedPositions.initialize(cc, maxContractedCopies*paddedParticles, elementSize, "rpmdContractedPositions");
    }

    // Create kernels.
    
    map<string, string> defines;
    defines["NUM_ATOMS"] = cc.intToString(cc.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    defines["NUM_COPIES"] = cc.intToString(numCopies);
    defines["THREAD_BLOCK_SIZE"] = cc.intToString(workgroupSize);
    defines["HBAR"] = cc.doubleToString(1.054571628e-34*AVOGADRO/(1000*1e-12), true);
    defines["SCALE"] = cc.doubleToString(1.0/sqrt((double) numCopies), true);
    defines["M_PI"] = cc.doubleToString(M_PI, true);
    map<string, string> replacements;
    replacements["FFT_Q_FORWARD"] = createFFT(numCopies, "q", true);
    replacements["FFT_Q_BACKWARD"] = createFFT(numCopies, "q", false);
    replacements["FFT_V_FORWARD"] = createFFT(numCopies, "v", true);
    replacements["FFT_V_BACKWARD"] = createFFT(numCopies, "v", false);
    ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonRpmdKernelSources::rpmd, replacements), defines);
    pileKernel = program->createKernel("applyPileThermostat");
    stepKernel = program->createKernel("integrateStep");
    velocitiesKernel = program->createKernel("advanceVelocities");
    computeCentroidKEKernel = program->createKernel("computeCentroidKE");
    applyBussiScalingKernel = program->createKernel("applyBussiScaling");
    copyToContextKernel = program->createKernel("copyDataToContext");
    copyFromContextKernel = program->createKernel("copyDataFromContext");
    translateKernel = program->createKernel("applyCellTranslations");
    
    // Create hybrid mode kernels (uniform layout with sync)
    if (hybridMode) {
        // Use hybrid versions that check isQuantum flag per particle
        pileKernelHybrid = program->createKernel("applyPileThermostatHybrid");
        stepKernelHybrid = program->createKernel("integrateStepHybrid");
        velocitiesKernelHybrid = program->createKernel("advanceVelocitiesHybrid");
        
        // Sync kernel to keep classical beads identical
        ComputeKernel syncClassicalBeadsKernel = program->createKernel("syncClassicalBeads");
        
        // Classical thermostat kernels
        applyClassicalThermostatKernel = program->createKernel("applyClassicalThermostat");
        computeClassicalKEKernel = program->createKernel("computeClassicalKEHybrid");
        applyBussiClassicalScalingKernel = program->createKernel("applyBussiClassicalHybrid");
    }
    
    // Allocate buffer for centroid kinetic energy reduction
    int numWorkGroups = (numParticles + workgroupSize - 1) / workgroupSize;
    centroidKE.initialize<double>(cc, numWorkGroups, "centroidKE");
    
    // Create kernels for doing contractions.
    
    for (auto& g : groupsByCopies) {
        int copies = g.first;
        replacements.clear();
        replacements["NUM_CONTRACTED_COPIES"] = cc.intToString(copies);
        replacements["POS_SCALE"] = cc.doubleToString(1.0/numCopies, true);
        replacements["FORCE_SCALE"] = cc.doubleToString(0x100000000/(double) copies, true);
        replacements["FFT_Q_FORWARD"] = createFFT(numCopies, "q", true);
        replacements["FFT_Q_BACKWARD"] = createFFT(copies, "q", false);
        replacements["FFT_F_FORWARD"] = createFFT(copies, "f", true);
        replacements["FFT_F_BACKWARD"] = createFFT(numCopies, "f", false);
        program = cc.compileProgram(cc.replaceStrings(CommonRpmdKernelSources::rpmdContraction, replacements), defines);
        positionContractionKernels[copies] = program->createKernel("contractPositions");
        forceContractionKernels[copies] = program->createKernel("contractForces");
    }
}

void CommonIntegrateRPMDStepKernel::initializeKernels(ContextImpl& context) {
    hasInitializedKernels = true;
    pileKernel->addArg(velocities);
    pileKernel->addArg(cc.getIntegrationUtilities().getRandom());
    pileKernel->addArg();
    pileKernel->addArg();
    pileKernel->addArg();
    pileKernel->addArg();
    pileKernel->addArg(); // Add 7th argument for applyToCentroid
    stepKernel->addArg(positions);
    stepKernel->addArg(velocities);
    stepKernel->addArg(forces);
    stepKernel->addArg();
    stepKernel->addArg();
    velocitiesKernel->addArg(velocities);
    velocitiesKernel->addArg(forces);
    velocitiesKernel->addArg();
    translateKernel->addArg(positions);
    translateKernel->addArg(cc.getPosq());
    translateKernel->addArg(cc.getAtomIndexArray());
    translateKernel->addArg();
    copyToContextKernel->addArg(velocities);
    copyToContextKernel->addArg(cc.getVelm());
    copyToContextKernel->addArg();
    copyToContextKernel->addArg(cc.getPosq());
    copyToContextKernel->addArg(cc.getAtomIndexArray());
    copyToContextKernel->addArg();
    copyFromContextKernel->addArg(cc.getLongForceBuffer());
    copyFromContextKernel->addArg();
    copyFromContextKernel->addArg(cc.getVelm());
    copyFromContextKernel->addArg(velocities);
    copyFromContextKernel->addArg(cc.getPosq());
    copyFromContextKernel->addArg();
    copyFromContextKernel->addArg(cc.getAtomIndexArray());
    copyFromContextKernel->addArg();
    
    // Initialize Bussi thermostat kernels
    computeCentroidKEKernel->addArg(velocities);
    computeCentroidKEKernel->addArg(centroidKE);
    applyBussiScalingKernel->addArg(velocities);
    applyBussiScalingKernel->addArg(); // alpha (will be set at runtime)
    
    // Initialize hybrid mode kernels (uniform layout with isQuantum check)
    if (hybridMode) {
        // PILE thermostat (hybrid) - checks isQuantum per particle
        pileKernelHybrid->addArg(velocities);
        pileKernelHybrid->addArg(cc.getIntegrationUtilities().getRandom());
        pileKernelHybrid->addArg(isQuantum);
        pileKernelHybrid->addArg(); // randomIndex
        pileKernelHybrid->addArg(); // dt
        pileKernelHybrid->addArg(); // kT
        pileKernelHybrid->addArg(); // friction
        pileKernelHybrid->addArg(); // applyToCentroid
        
        // Integration step (hybrid) - quantum: FFT, classical: Verlet
        stepKernelHybrid->addArg(positions);
        stepKernelHybrid->addArg(velocities);
        stepKernelHybrid->addArg(forces);
        stepKernelHybrid->addArg(isQuantum);
        stepKernelHybrid->addArg(); // dt
        stepKernelHybrid->addArg(); // kT
        
        // Velocity advance (hybrid)
        velocitiesKernelHybrid->addArg(velocities);
        velocitiesKernelHybrid->addArg(forces);
        velocitiesKernelHybrid->addArg(isQuantum);
        velocitiesKernelHybrid->addArg(); // dt
        
        // Classical thermostat (Langevin) - uses isQuantum to skip quantum particles
        applyClassicalThermostatKernel->addArg(velocities);
        applyClassicalThermostatKernel->addArg(cc.getIntegrationUtilities().getRandom());
        applyClassicalThermostatKernel->addArg(isQuantum);
        applyClassicalThermostatKernel->addArg(); // randomIndex
        applyClassicalThermostatKernel->addArg(); // dt
        applyClassicalThermostatKernel->addArg(); // kT
        applyClassicalThermostatKernel->addArg(); // friction
        applyClassicalThermostatKernel->addArg(); // thermostatType
        
        // Classical KE for Bussi
        computeClassicalKEKernel->addArg(velocities);
        computeClassicalKEKernel->addArg(isQuantum);
        computeClassicalKEKernel->addArg(classicalKE);
        
        // Bussi scaling for classical
        applyBussiClassicalScalingKernel->addArg(velocities);
        applyBussiClassicalScalingKernel->addArg(isQuantum);
        applyBussiClassicalScalingKernel->addArg(); // scalingFactor
    }
    
    for (auto& g : groupsByCopies) {
        int copies = g.first;
        positionContractionKernels[copies]->addArg(positions);
        positionContractionKernels[copies]->addArg(contractedPositions);
        forceContractionKernels[copies]->addArg(forces);
        forceContractionKernels[copies]->addArg(contractedForces);
    }

    // #region agent log
    {
        std::ostringstream data;
        data << "{\"positionsSize\":" << positions.getSize()
             << ",\"velocitiesSize\":" << velocities.getSize()
             << ",\"forcesSize\":" << forces.getSize()
             << ",\"centroidKESize\":" << centroidKE.getSize()
             << ",\"groupsByCopies\":" << groupsByCopies.size() << "}";
        appendDebugLog("CommonRpmdKernels.cpp:initializeKernels", "kernel args set", data.str(), "H4", "pre-fix");
    }
    // #endregion
}

void CommonIntegrateRPMDStepKernel::execute(ContextImpl& context, const RPMDIntegrator& integrator, bool forcesAreValid) {
    ContextSelector selector(cc);
    if (!hasInitializedKernels)
        initializeKernels(context);
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    RPMDIntegrator::ThermostatType thermostatType = integrator.getThermostatType();
    bool applyThermostat = integrator.getApplyThermostat() && (thermostatType != RPMDIntegrator::NoneThermo);
    bool hybridMode = (numClassicalParticles > 0);

    // #region agent log
    {
        std::ostringstream data;
        data << "{\"thermostatType\":" << static_cast<int>(thermostatType)
             << ",\"applyThermostat\":" << (applyThermostat ? 1 : 0)
             << ",\"hybridMode\":" << (hybridMode ? 1 : 0)
             << ",\"numQuantumParticles\":" << numQuantumParticles
             << ",\"numClassicalParticles\":" << numClassicalParticles
             << ",\"hasInitializedKernels\":" << (hasInitializedKernels ? 1 : 0) << "}";
        appendDebugLog("CommonRpmdKernels.cpp:execute", "exec start", data.str(), "hybrid-v1", "run1");
    }
    // #endregion
    
    // Compute forces on all particles (all beads)
    if (!forcesAreValid)
        computeForces(context);
    
    bool useDoublePrecision = (cc.getUseDoublePrecision() || cc.getUseMixedPrecision());
    double dt = integrator.getStepSize();
    int applyToCentroid = (thermostatType == RPMDIntegrator::Pile) ? 1 : 0;
    
    if (hybridMode) {
        // ====================================================================
        // HYBRID MODE: Uniform storage with quantum/classical distinction
        // Quantum particles: full RPMD with FFT (PILE thermostat)
        // Classical particles: simple dynamics (Langevin/Bussi thermostat)
        // All particles stored with all beads; classical beads synced to stay identical
        // ====================================================================
        
        // Set kernel arguments for hybrid kernels
        pileKernelHybrid->setArg(3, integration.prepareRandomNumbers(numParticles*numCopies));
        if (useDoublePrecision) {
            pileKernelHybrid->setArg(4, dt);
            pileKernelHybrid->setArg(5, integrator.getTemperature()*BOLTZ);
            pileKernelHybrid->setArg(6, integrator.getFriction());
            pileKernelHybrid->setArg(7, applyToCentroid);
            stepKernelHybrid->setArg(4, dt);
            stepKernelHybrid->setArg(5, integrator.getTemperature()*BOLTZ);
            velocitiesKernelHybrid->setArg(3, dt);
        } else {
            pileKernelHybrid->setArg(4, (float) dt);
            pileKernelHybrid->setArg(5, (float) (integrator.getTemperature()*BOLTZ));
            pileKernelHybrid->setArg(6, (float) integrator.getFriction());
            pileKernelHybrid->setArg(7, applyToCentroid);
            stepKernelHybrid->setArg(4, (float) dt);
            stepKernelHybrid->setArg(5, (float) (integrator.getTemperature()*BOLTZ));
            velocitiesKernelHybrid->setArg(3, (float) dt);
        }
        
        // THERMOSTAT (FIRST HALF)
        if (applyThermostat) {
            if (thermostatType == RPMDIntegrator::PileG) {
                applyBussiCentroidThermostat(context.getSystem(), integrator, dt*0.5);
            }
            // Apply PILE to quantum particles only (hybrid kernel checks isQuantum)
            pileKernelHybrid->execute(numParticles*numCopies, workgroupSize);
            
            // Apply classical thermostat
            RPMDIntegrator::ClassicalThermostatType classicalThermostat = integrator.getClassicalThermostat();
            if (classicalThermostat == RPMDIntegrator::BussiClassical) {
                applyBussiClassicalThermostat(context.getSystem(), integrator, dt*0.5);
            } else if (classicalThermostat == RPMDIntegrator::LangevinClassical) {
                applyClassicalThermostatKernel->setArg(3, integration.prepareRandomNumbers(numParticles));
                if (useDoublePrecision) {
                    applyClassicalThermostatKernel->setArg(4, dt);
                    applyClassicalThermostatKernel->setArg(5, integrator.getTemperature()*BOLTZ);
                    applyClassicalThermostatKernel->setArg(6, integrator.getFriction());
                } else {
                    applyClassicalThermostatKernel->setArg(4, (float) dt);
                    applyClassicalThermostatKernel->setArg(5, (float) (integrator.getTemperature()*BOLTZ));
                    applyClassicalThermostatKernel->setArg(6, (float) integrator.getFriction());
                }
                applyClassicalThermostatKernel->setArg(7, 1);  // Langevin mode
                applyClassicalThermostatKernel->execute(numParticles);
            }
        }
        
        // INTEGRATION (hybrid kernel handles both quantum and classical)
        stepKernelHybrid->execute(numParticles*numCopies, workgroupSize);
        
        // FORCE COMPUTATION
        computeForces(context);
        
        // VELOCITY UPDATE (SECOND HALF-STEP)
        velocitiesKernelHybrid->execute(numParticles*numCopies, workgroupSize);
        
        // THERMOSTAT (SECOND HALF)
        if (applyThermostat) {
            if (thermostatType == RPMDIntegrator::PileG) {
                applyBussiCentroidThermostat(context.getSystem(), integrator, dt*0.5);
            }
            pileKernelHybrid->setArg(3, integration.prepareRandomNumbers(numParticles*numCopies));
            pileKernelHybrid->execute(numParticles*numCopies, workgroupSize);
            
            // Apply classical thermostat again
            RPMDIntegrator::ClassicalThermostatType classicalThermostat = integrator.getClassicalThermostat();
            if (classicalThermostat == RPMDIntegrator::BussiClassical) {
                applyBussiClassicalThermostat(context.getSystem(), integrator, dt*0.5);
            } else if (classicalThermostat == RPMDIntegrator::LangevinClassical) {
                applyClassicalThermostatKernel->setArg(3, integration.prepareRandomNumbers(numParticles));
                applyClassicalThermostatKernel->setArg(7, 1);
                applyClassicalThermostatKernel->execute(numParticles);
            }
        }
    } else {
        // ====================================================================
        // STANDARD RPMD MODE (all quantum, no sparse storage)
        // ====================================================================
        
        // Set kernel arguments for standard RPMD kernels
        pileKernel->setArg(2, integration.prepareRandomNumbers(numParticles*numCopies));
        if (useDoublePrecision) {
            pileKernel->setArg(3, dt);
            pileKernel->setArg(4, integrator.getTemperature()*BOLTZ);
            pileKernel->setArg(5, integrator.getFriction());
            pileKernel->setArg(6, applyToCentroid);
            stepKernel->setArg(3, dt);
            stepKernel->setArg(4, integrator.getTemperature()*BOLTZ);
            velocitiesKernel->setArg(2, dt);
        } else {
            pileKernel->setArg(3, (float) dt);
            pileKernel->setArg(4, (float) (integrator.getTemperature()*BOLTZ));
            pileKernel->setArg(5, (float) integrator.getFriction());
            pileKernel->setArg(6, applyToCentroid);
            stepKernel->setArg(3, (float) dt);
            stepKernel->setArg(4, (float) (integrator.getTemperature()*BOLTZ));
            velocitiesKernel->setArg(2, (float) dt);
        }
        
        // THERMOSTAT (FIRST HALF)
        if (applyThermostat) {
            if (thermostatType == RPMDIntegrator::PileG) {
                applyBussiCentroidThermostat(context.getSystem(), integrator, dt*0.5);
            }
            pileKernel->execute(numParticles*numCopies, workgroupSize);
        }
        
        // INTEGRATION
        stepKernel->execute(numParticles*numCopies, workgroupSize);
        
        // FORCE COMPUTATION
        computeForces(context);
        
        // VELOCITY UPDATE (SECOND HALF-STEP)
        velocitiesKernel->execute(numParticles*numCopies, workgroupSize);
        
        // THERMOSTAT (SECOND HALF)
        if (applyThermostat) {
            if (thermostatType == RPMDIntegrator::PileG) {
                applyBussiCentroidThermostat(context.getSystem(), integrator, dt*0.5);
            }
            pileKernel->setArg(2, integration.prepareRandomNumbers(numParticles*numCopies));
            pileKernel->execute(numParticles*numCopies, workgroupSize);
        }
    }

    // Update the time and step count
    cc.setTime(cc.getTime()+dt);
    cc.setStepCount(cc.getStepCount()+1);
    cc.reorderAtoms();
    if (cc.getAtomsWereReordered() && cc.getNonbondedUtilities().getUsePeriodic()) {
        // Atoms may have been translated into a different periodic box, so apply
        // the same translation to all the beads.
        translateKernel->setArg(3, numCopies-1);
        translateKernel->execute(cc.getNumAtoms());
    }
}

void CommonIntegrateRPMDStepKernel::computeForces(ContextImpl& context) {
    // Check if we can use batched force evaluation
    bool useBatchedEvaluation = false;
    const PythonForce* batchedPythonForce = nullptr;
    int batchedForceIndex = -1;
    const System& system = context.getSystem();
    
    // Check all forces in groupsNotContracted to see if any support batching
    for (int i = 0; i < system.getNumForces(); i++) {
        if ((groupsNotContracted & (1<<system.getForce(i).getForceGroup())) != 0) {
            if (const PythonForce* pythonForce = dynamic_cast<const PythonForce*>(&system.getForce(i))) {
                if (pythonForce->supportsBatchedEvaluation()) {
                    useBatchedEvaluation = true;
                    batchedPythonForce = pythonForce;
                    batchedForceIndex = i;
                    break;
                }
            }
        }
    }
    
    if (useBatchedEvaluation && batchedPythonForce != nullptr) {
        // NEW BATCHED PATH: Evaluate all copies at once
        // printf("DEBUG: Using batched RPMD evaluation for %d beads\n", numCopies);
        
        // 1. Collect all bead positions from GPU
        std::vector<std::vector<Vec3>> allBeadPositions(numCopies);
        for (int i = 0; i < numCopies; i++) {
            downloadPositionsFromGPU(i, allBeadPositions[i]);
        }
        
        // 2. Get the PythonForce implementation and call batched evaluation
        const std::vector<ForceImpl*>& forceImpls = context.getForceImpls();
        ForceImpl& forceImpl = *forceImpls[batchedForceIndex];
        PythonForceImpl& pyImpl = dynamic_cast<PythonForceImpl&>(forceImpl);
        
        // Allocate output forces
        std::vector<std::vector<Vec3>> allBeadForces(numCopies);
        for (auto& beadForces : allBeadForces) {
            beadForces.resize(numParticles);
        }
        
        // Call batched force evaluation (single Python call for all beads!)
        double energy = pyImpl.calcForcesAndEnergyBatched(context, allBeadPositions, allBeadForces);
        
        // 3. Upload all forces back to GPU in one shot
        uploadAllForcesToGPU(allBeadForces);
        
        // 4. Handle any non-PythonForce forces sequentially (if any exist)
        for (int i = 0; i < system.getNumForces(); i++) {
            if ((groupsNotContracted & (1<<system.getForce(i).getForceGroup())) != 0) {
                if (i != batchedForceIndex) {
                    // Skip CMMotionRemover - it doesn't compute forces
                    if (dynamic_cast<const CMMotionRemover*>(&system.getForce(i)) != nullptr) {
                        // printf("DEBUG: Skipping CMMotionRemover (index %d)\n", i);
                        continue;
                    }
                    
                    // Non-batched force: use original sequential path
                    // printf("DEBUG: Force %d (type %s) using sequential path\n", i, typeid(system.getForce(i)).name());
                    int forceGroupFlag = 1 << system.getForce(i).getForceGroup();
                    copyToContextKernel->setArg(2, positions);
                    copyFromContextKernel->setArg(1, forces);
                    copyFromContextKernel->setArg(5, positions);
                    for (int bead = 0; bead < numCopies; bead++) {
                        copyToContextKernel->setArg(5, bead);
                        copyToContextKernel->execute(cc.getNumAtoms());
                        context.computeVirtualSites();
                        context.updateContextState();
                        context.calcForcesAndEnergy(true, false, forceGroupFlag);
                        copyFromContextKernel->setArg(7, bead);
                        copyFromContextKernel->execute(cc.getNumAtoms());
                    }
                }
            }
        }
    } else {
        // ORIGINAL PATH: Compute forces from all groups that didn't have a specified contraction.
        // Uses standard copy kernels - uniform layout works for both hybrid and standard modes.
        // printf("DEBUG: Using sequential RPMD evaluation (%d beads)\n", numCopies);

        copyToContextKernel->setArg(2, positions);
        copyFromContextKernel->setArg(1, forces);
        copyFromContextKernel->setArg(5, positions);
        for (int i = 0; i < numCopies; i++) {
            copyToContextKernel->setArg(5, i);
            copyToContextKernel->execute(cc.getNumAtoms());
            context.computeVirtualSites();
            Vec3 initialBox[3];
            context.getPeriodicBoxVectors(initialBox[0], initialBox[1], initialBox[2]);
            context.updateContextState();
            Vec3 finalBox[3];
            context.getPeriodicBoxVectors(finalBox[0], finalBox[1], finalBox[2]);
            if (initialBox[0] != finalBox[0] || initialBox[1] != finalBox[1] || initialBox[2] != finalBox[2])
                throw OpenMMException("Standard barostats cannot be used with RPMDIntegrator.  Use RPMDMonteCarloBarostat instead.");
            context.calcForcesAndEnergy(true, false, groupsNotContracted);
            copyFromContextKernel->setArg(7, i);
            copyFromContextKernel->execute(cc.getNumAtoms());
        }
    }
    
    // Now loop over contractions and compute forces from them.
    
    if (groupsByCopies.size() > 0) {
        copyToContextKernel->setArg(2, contractedPositions);
        copyFromContextKernel->setArg(1, contractedForces);
        copyFromContextKernel->setArg(5, contractedPositions);
        for (auto& g : groupsByCopies) {
            int copies = g.first;
            int groupFlags = g.second;

            // Find the contracted positions.

           positionContractionKernels[copies]->execute(numParticles*numCopies, workgroupSize);

            // Compute forces.

            for (int i = 0; i < copies; i++) {
                copyToContextKernel->setArg(5, i);
                copyToContextKernel->execute(cc.getNumAtoms());
                context.computeVirtualSites();
                context.calcForcesAndEnergy(true, false, groupFlags);
                copyFromContextKernel->setArg(7, i);
                copyFromContextKernel->execute(cc.getNumAtoms());
            }

            // Apply the forces to the original copies.

            forceContractionKernels[copies]->execute(numParticles*numCopies, workgroupSize);
        }
    }
    if (groupsByCopies.size() > 0) {
        // Ensure the Context contains the positions from the last copy, since we'll assume that later.
        
        copyToContextKernel->setArg(2, positions);
        copyToContextKernel->setArg(5, numCopies-1);
        copyToContextKernel->execute(cc.getNumAtoms());
    }
}

double CommonIntegrateRPMDStepKernel::computeKineticEnergy(ContextImpl& context, const RPMDIntegrator& integrator) {
    return cc.getIntegrationUtilities().computeKineticEnergy(0);
}

void CommonIntegrateRPMDStepKernel::applyBussiCentroidThermostat(const System& system, const RPMDIntegrator& integrator, double halfdt) {
    // Bussi stochastic velocity rescaling thermostat for centroid mode ONLY
    // Reference: Bussi, Donadio, Parrinello, J. Chem. Phys. 126, 014101 (2007)
    //
    // This implementation uses GPU kernels to compute centroid kinetic energy
    // and apply the Bussi scaling factor, keeping all data on the GPU.
    
    const double kT = BOLTZ * integrator.getTemperature();
    const double c1 = exp(-halfdt * integrator.getCentroidFriction());
    bool useDoublePrecision = (cc.getUseDoublePrecision() || cc.getUseMixedPrecision());
    
    // Step 1: Compute centroid kinetic energy on GPU
    // Arguments already set in initializeKernels()
    computeCentroidKEKernel->execute(numParticles, workgroupSize);
    
    // Step 2: Download kinetic energy and compute scaling factor on CPU
    int numWorkGroups = (numParticles + workgroupSize - 1) / workgroupSize;
    vector<double> keData(numWorkGroups);
    centroidKE.download(keData);
    
    double totalKE = 0.0;
    for (int i = 0; i < numWorkGroups; i++)
        totalKE += keData[i];
    
    // Count degrees of freedom
    int ndof = 0;
    for (int i = 0; i < numParticles; i++) {
        if (system.getParticleMass(i) > 0.0)
            ndof += 3;
    }

    // #region agent log
    {
        std::ostringstream data;
        data << "{\"numParticles\":" << numParticles
             << ",\"workgroupSize\":" << workgroupSize
             << ",\"numWorkGroups\":" << numWorkGroups
             << ",\"centroidKESize\":" << centroidKE.getSize()
             << ",\"totalKE\":" << totalKE
             << ",\"ndof\":" << ndof
             << ",\"useDoublePrecision\":" << (useDoublePrecision ? 1 : 0) << "}";
        appendDebugLog("CommonRpmdKernels.cpp:applyBussiCentroidThermostat", "centroid KE", data.str(), "H1", "pre-fix");
    }
    // #endregion

    if (totalKE <= 0.0 || ndof == 0)
        return;
    
    // Calculate Bussi rescaling factor using Gaussian random numbers
    double K_target = 0.5 * ndof * kT;
    
    // Generate Gaussian random numbers on CPU using standard library
    static thread_local std::mt19937 rng;
    static thread_local std::normal_distribution<double> normal(0.0, 1.0);
    
    double R1 = normal(rng);
    double R_gamma = 0.0;
    for (int i = 0; i < ndof - 1; i++) {
        double rnd = normal(rng);
        R_gamma += rnd * rnd;
    }
    
    double ratio = K_target / totalKE;
    double alpha2 = c1 + ratio * (1.0 - c1) * (R1 * R1 + R_gamma)
                    + 2.0 * R1 * sqrt(ratio * c1 * (1.0 - c1));
    
    if (alpha2 <= 0.0)
        return;
    
    double alpha = sqrt(alpha2);
    if (R1 + sqrt(2.0 * totalKE / K_target * c1 / (1.0 - c1)) < 0.0)
        alpha = -alpha;
    
    // Step 3: Apply scaling on GPU
    // velocities already set in initializeKernels(), only need to update alpha
    if (useDoublePrecision)
        applyBussiScalingKernel->setArg(1, alpha);
    else
        applyBussiScalingKernel->setArg(1, (float) alpha);
    
    applyBussiScalingKernel->execute(numParticles);
}

void CommonIntegrateRPMDStepKernel::applyBussiClassicalThermostat(const System& system, const RPMDIntegrator& integrator, double halfdt) {
    // Bussi stochastic velocity rescaling thermostat for CLASSICAL particles only
    // Reference: Bussi, Donadio, Parrinello, J. Chem. Phys. 126, 014101 (2007)
    
    if (!hybridMode || numClassicalParticles == 0)
        return;
    
    const double kT = BOLTZ * integrator.getTemperature();
    const double c1 = exp(-halfdt * integrator.getFriction());  // Use main friction for classical
    bool useDoublePrecision = (cc.getUseDoublePrecision() || cc.getUseMixedPrecision());
    
    // Step 1: Compute classical kinetic energy on GPU (uses isQuantum to identify classical particles)
    int numWorkGroups = (numParticles + workgroupSize - 1) / workgroupSize;
    computeClassicalKEKernel->execute(numParticles, workgroupSize);
    
    // Step 2: Download kinetic energy and compute scaling factor on CPU
    vector<double> keData(numWorkGroups);
    classicalKE.download(keData);
    
    double totalKE = 0.0;
    for (int i = 0; i < numWorkGroups; i++)
        totalKE += keData[i];
    
    // Count degrees of freedom for classical particles
    int ndof = numClassicalParticles * 3;  // Classical particles always have 3 DOF each
    
    if (totalKE <= 0.0 || ndof == 0)
        return;
    
    // Calculate Bussi rescaling factor using Gaussian random numbers
    double K_target = 0.5 * ndof * kT;
    
    // Generate Gaussian random numbers on CPU
    static thread_local std::mt19937 rng;
    static thread_local std::normal_distribution<double> normal(0.0, 1.0);
    
    double R1 = normal(rng);
    double R_gamma = 0.0;
    for (int i = 0; i < ndof - 1; i++) {
        double rnd = normal(rng);
        R_gamma += rnd * rnd;
    }
    
    double ratio = K_target / totalKE;
    double alpha2 = c1 + ratio * (1.0 - c1) * (R1 * R1 + R_gamma)
                    + 2.0 * R1 * sqrt(ratio * c1 * (1.0 - c1));
    
    if (alpha2 <= 0.0)
        return;
    
    double alpha = sqrt(alpha2);
    if (R1 + sqrt(2.0 * totalKE / K_target * c1 / (1.0 - c1)) < 0.0)
        alpha = -alpha;
    
    // Step 3: Apply scaling on GPU (uses isQuantum to identify classical particles)
    if (useDoublePrecision)
        applyBussiClassicalScalingKernel->setArg(2, alpha);
    else
        applyBussiClassicalScalingKernel->setArg(2, (float) alpha);
    
    applyBussiClassicalScalingKernel->execute(numParticles);
}

void CommonIntegrateRPMDStepKernel::setPositions(int copy, const vector<Vec3>& pos) {
    if (!positions.isInitialized())
        throw OpenMMException("RPMDIntegrator: Cannot set positions before the integrator is added to a Context");
    if (pos.size() != numParticles)
        throw OpenMMException("RPMDIntegrator: wrong number of values passed to setPositions()");

    // Adjust the positions based on the current cell offsets.
    
    const vector<int>& order = cc.getAtomIndex();
    Vec3 a, b, c;
    cc.getPeriodicBoxVectors(a, b, c);
    vector<Vec3> offsetPos(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        mm_int4 offset = cc.getPosCellOffsets()[i];
        offsetPos[order[i]] = pos[order[i]] + Vec3(offset.x*a[0], offset.y*b[1], offset.z*c[2]);
    }

    // Record the positions (uniform layout - all particles × all beads)
    ContextSelector selector(cc);
    if (cc.getUseDoublePrecision()) {
        vector<mm_double4> posq(cc.getPaddedNumAtoms());
        cc.getPosq().download(posq);
        for (int i = 0; i < numParticles; i++)
            posq[i] = mm_double4(offsetPos[i][0], offsetPos[i][1], offsetPos[i][2], posq[i].w);
        positions.uploadSubArray(&posq[0], copy*cc.getPaddedNumAtoms(), numParticles);
    }
    else if (cc.getUseMixedPrecision()) {
        vector<mm_float4> posqf(cc.getPaddedNumAtoms());
        cc.getPosq().download(posqf);
        vector<mm_double4> posq(cc.getPaddedNumAtoms());
        for (int i = 0; i < numParticles; i++)
            posq[i] = mm_double4(offsetPos[i][0], offsetPos[i][1], offsetPos[i][2], posqf[i].w);
        positions.uploadSubArray(&posq[0], copy*cc.getPaddedNumAtoms(), numParticles);
    }
    else {
        vector<mm_float4> posq(cc.getPaddedNumAtoms());
        cc.getPosq().download(posq);
        for (int i = 0; i < numParticles; i++)
            posq[i] = mm_float4((float) offsetPos[i][0], (float) offsetPos[i][1], (float) offsetPos[i][2], posq[i].w);
        positions.uploadSubArray(&posq[0], copy*cc.getPaddedNumAtoms(), numParticles);
    }
}

void CommonIntegrateRPMDStepKernel::setVelocities(int copy, const vector<Vec3>& vel) {
    if (!velocities.isInitialized())
        throw OpenMMException("RPMDIntegrator: Cannot set velocities before the integrator is added to a Context");
    if (vel.size() != numParticles)
        throw OpenMMException("RPMDIntegrator: wrong number of values passed to setVelocities()");
    ContextSelector selector(cc);
    
    // Uniform layout - all particles × all beads
    if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
        vector<mm_double4> velm(cc.getPaddedNumAtoms());
        cc.getVelm().download(velm);
        for (int i = 0; i < numParticles; i++)
            velm[i] = mm_double4(vel[i][0], vel[i][1], vel[i][2], velm[i].w);
        velocities.uploadSubArray(&velm[0], copy*cc.getPaddedNumAtoms(), numParticles);
    }
    else {
        vector<mm_float4> velm(cc.getPaddedNumAtoms());
        cc.getVelm().download(velm);
        for (int i = 0; i < numParticles; i++)
            velm[i] = mm_float4((float) vel[i][0], (float) vel[i][1], (float) vel[i][2], velm[i].w);
        velocities.uploadSubArray(&velm[0], copy*cc.getPaddedNumAtoms(), numParticles);
    }
}

void CommonIntegrateRPMDStepKernel::copyToContext(int copy, ContextImpl& context) {
    ContextSelector selector(cc);
    if (!hasInitializedKernels)
        initializeKernels(context);
    copyToContextKernel->setArg(2, positions);
    copyToContextKernel->setArg(5, copy);
    copyToContextKernel->execute(cc.getNumAtoms());
}

string CommonIntegrateRPMDStepKernel::createFFT(int size, const string& variable, bool forward) {
    stringstream source;
    int stage = 0;
    int L = size;
    int m = 1;
    string sign = (forward ? "1.0f" : "-1.0f");
    string multReal = (forward ? "multiplyComplexRealPart" : "multiplyComplexRealPartConj");
    string multImag = (forward ? "multiplyComplexImagPart" : "multiplyComplexImagPartConj");

    source<<"{\n";
    source<<"LOCAL_ARG mixed3* real0 = "<<variable<<"real;\n";
    source<<"LOCAL_ARG mixed3* imag0 = "<<variable<<"imag;\n";
    source<<"LOCAL_ARG mixed3* real1 = &temp[blockStart];\n";
    source<<"LOCAL_ARG mixed3* imag1 = &temp[blockStart+LOCAL_SIZE];\n";

    // Factor size, generating an appropriate block of code for each factor.

    while (L > 1) {
        int input = stage%2;
        int output = 1-input;
        int radix;
        if (L%5 == 0)
            radix = 5;
        else if (L%4 == 0)
            radix = 4;
        else if (L%3 == 0)
            radix = 3;
        else if (L%2 == 0)
            radix = 2;
        else
            throw OpenMMException("Illegal size for FFT: "+cc.intToString(size));
        source<<"{\n";
        L = L/radix;
        source<<"// Pass "<<(stage+1)<<" (radix "<<radix<<")\n";
        source<<"if (indexInBlock < "<<(L*m)<<") {\n";
        source<<"int i = indexInBlock;\n";
        source<<"int j = i/"<<m<<";\n";
        if (radix == 5) {
            source<<"mixed3 c0r = real"<<input<<"[i];\n";
            source<<"mixed3 c0i = imag"<<input<<"[i];\n";
            source<<"mixed3 c1r = real"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"mixed3 c1i = imag"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"mixed3 c2r = real"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"mixed3 c2i = imag"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"mixed3 c3r = real"<<input<<"[i+"<<(3*L*m)<<"];\n";
            source<<"mixed3 c3i = imag"<<input<<"[i+"<<(3*L*m)<<"];\n";
            source<<"mixed3 c4r = real"<<input<<"[i+"<<(4*L*m)<<"];\n";
            source<<"mixed3 c4i = imag"<<input<<"[i+"<<(4*L*m)<<"];\n";
            source<<"mixed3 d0r = c1r+c4r;\n";
            source<<"mixed3 d0i = c1i+c4i;\n";
            source<<"mixed3 d1r = c2r+c3r;\n";
            source<<"mixed3 d1i = c2i+c3i;\n";
            source<<"mixed3 d2r = "<<cc.doubleToString(sin(0.4*M_PI), true)<<"*(c1r-c4r);\n";
            source<<"mixed3 d2i = "<<cc.doubleToString(sin(0.4*M_PI), true)<<"*(c1i-c4i);\n";
            source<<"mixed3 d3r = "<<cc.doubleToString(sin(0.4*M_PI), true)<<"*(c2r-c3r);\n";
            source<<"mixed3 d3i = "<<cc.doubleToString(sin(0.4*M_PI), true)<<"*(c2i-c3i);\n";
            source<<"mixed3 d4r = d0r+d1r;\n";
            source<<"mixed3 d4i = d0i+d1i;\n";
            source<<"mixed3 d5r = "<<cc.doubleToString(0.25*sqrt(5.0), true)<<"*(d0r-d1r);\n";
            source<<"mixed3 d5i = "<<cc.doubleToString(0.25*sqrt(5.0), true)<<"*(d0i-d1i);\n";
            source<<"mixed3 d6r = c0r-0.25f*d4r;\n";
            source<<"mixed3 d6i = c0i-0.25f*d4i;\n";
            source<<"mixed3 d7r = d6r+d5r;\n";
            source<<"mixed3 d7i = d6i+d5i;\n";
            source<<"mixed3 d8r = d6r-d5r;\n";
            source<<"mixed3 d8i = d6i-d5i;\n";
            string coeff = cc.doubleToString(sin(0.2*M_PI)/sin(0.4*M_PI), true);
            source<<"mixed3 d9r = "<<sign<<"*(d2i+"<<coeff<<"*d3i);\n";
            source<<"mixed3 d9i = "<<sign<<"*(-d2r-"<<coeff<<"*d3r);\n";
            source<<"mixed3 d10r = "<<sign<<"*("<<coeff<<"*d2i-d3i);\n";
            source<<"mixed3 d10i = "<<sign<<"*(d3r-"<<coeff<<"*d2r);\n";
            source<<"real"<<output<<"[i+4*j*"<<m<<"] = c0r+d4r;\n";
            source<<"imag"<<output<<"[i+4*j*"<<m<<"] = c0i+d4i;\n";
            source<<"real"<<output<<"[i+(4*j+1)*"<<m<<"] = "<<multReal<<"(w[j*"<<size<<"/"<<(5*L)<<"], d7r+d9r, d7i+d9i);\n";
            source<<"imag"<<output<<"[i+(4*j+1)*"<<m<<"] = "<<multImag<<"(w[j*"<<size<<"/"<<(5*L)<<"], d7r+d9r, d7i+d9i);\n";
            source<<"real"<<output<<"[i+(4*j+2)*"<<m<<"] = "<<multReal<<"(w[j*"<<(2*size)<<"/"<<(5*L)<<"], d8r+d10r, d8i+d10i);\n";
            source<<"imag"<<output<<"[i+(4*j+2)*"<<m<<"] = "<<multImag<<"(w[j*"<<(2*size)<<"/"<<(5*L)<<"], d8r+d10r, d8i+d10i);\n";
            source<<"real"<<output<<"[i+(4*j+3)*"<<m<<"] = "<<multReal<<"(w[j*"<<(3*size)<<"/"<<(5*L)<<"], d8r-d10r, d8i-d10i);\n";
            source<<"imag"<<output<<"[i+(4*j+3)*"<<m<<"] = "<<multImag<<"(w[j*"<<(3*size)<<"/"<<(5*L)<<"], d8r-d10r, d8i-d10i);\n";
            source<<"real"<<output<<"[i+(4*j+4)*"<<m<<"] = "<<multReal<<"(w[j*"<<(4*size)<<"/"<<(5*L)<<"], d7r-d9r, d7i-d9i);\n";
            source<<"imag"<<output<<"[i+(4*j+4)*"<<m<<"] = "<<multImag<<"(w[j*"<<(4*size)<<"/"<<(5*L)<<"], d7r-d9r, d7i-d9i);\n";
        }
        else if (radix == 4) {
            source<<"mixed3 c0r = real"<<input<<"[i];\n";
            source<<"mixed3 c0i = imag"<<input<<"[i];\n";
            source<<"mixed3 c1r = real"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"mixed3 c1i = imag"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"mixed3 c2r = real"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"mixed3 c2i = imag"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"mixed3 c3r = real"<<input<<"[i+"<<(3*L*m)<<"];\n";
            source<<"mixed3 c3i = imag"<<input<<"[i+"<<(3*L*m)<<"];\n";
            source<<"mixed3 d0r = c0r+c2r;\n";
            source<<"mixed3 d0i = c0i+c2i;\n";
            source<<"mixed3 d1r = c0r-c2r;\n";
            source<<"mixed3 d1i = c0i-c2i;\n";
            source<<"mixed3 d2r = c1r+c3r;\n";
            source<<"mixed3 d2i = c1i+c3i;\n";
            source<<"mixed3 d3r = "<<sign<<"*(c1i-c3i);\n";
            source<<"mixed3 d3i = "<<sign<<"*(c3r-c1r);\n";
            source<<"real"<<output<<"[i+3*j*"<<m<<"] = d0r+d2r;\n";
            source<<"imag"<<output<<"[i+3*j*"<<m<<"] = d0i+d2i;\n";
            source<<"real"<<output<<"[i+(3*j+1)*"<<m<<"] = "<<multReal<<"(w[j*"<<size<<"/"<<(4*L)<<"], d1r+d3r, d1i+d3i);\n";
            source<<"imag"<<output<<"[i+(3*j+1)*"<<m<<"] = "<<multImag<<"(w[j*"<<size<<"/"<<(4*L)<<"], d1r+d3r, d1i+d3i);\n";
            source<<"real"<<output<<"[i+(3*j+2)*"<<m<<"] = "<<multReal<<"(w[j*"<<(2*size)<<"/"<<(4*L)<<"], d0r-d2r, d0i-d2i);\n";
            source<<"imag"<<output<<"[i+(3*j+2)*"<<m<<"] = "<<multImag<<"(w[j*"<<(2*size)<<"/"<<(4*L)<<"], d0r-d2r, d0i-d2i);\n";
            source<<"real"<<output<<"[i+(3*j+3)*"<<m<<"] = "<<multReal<<"(w[j*"<<(3*size)<<"/"<<(4*L)<<"], d1r-d3r, d1i-d3i);\n";
            source<<"imag"<<output<<"[i+(3*j+3)*"<<m<<"] = "<<multImag<<"(w[j*"<<(3*size)<<"/"<<(4*L)<<"], d1r-d3r, d1i-d3i);\n";
        }
        else if (radix == 3) {
            source<<"mixed3 c0r = real"<<input<<"[i];\n";
            source<<"mixed3 c0i = imag"<<input<<"[i];\n";
            source<<"mixed3 c1r = real"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"mixed3 c1i = imag"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"mixed3 c2r = real"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"mixed3 c2i = imag"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"mixed3 d0r = c1r+c2r;\n";
            source<<"mixed3 d0i = c1i+c2i;\n";
            source<<"mixed3 d1r = c0r-0.5f*d0r;\n";
            source<<"mixed3 d1i = c0i-0.5f*d0i;\n";
            source<<"mixed3 d2r = "<<sign<<"*"<<cc.doubleToString(sin(M_PI/3.0), true)<<"*(c1i-c2i);\n";
            source<<"mixed3 d2i = "<<sign<<"*"<<cc.doubleToString(sin(M_PI/3.0), true)<<"*(c2r-c1r);\n";
            source<<"real"<<output<<"[i+2*j*"<<m<<"] = c0r+d0r;\n";
            source<<"imag"<<output<<"[i+2*j*"<<m<<"] = c0i+d0i;\n";
            source<<"real"<<output<<"[i+(2*j+1)*"<<m<<"] = "<<multReal<<"(w[j*"<<size<<"/"<<(3*L)<<"], d1r+d2r, d1i+d2i);\n";
            source<<"imag"<<output<<"[i+(2*j+1)*"<<m<<"] = "<<multImag<<"(w[j*"<<size<<"/"<<(3*L)<<"], d1r+d2r, d1i+d2i);\n";
            source<<"real"<<output<<"[i+(2*j+2)*"<<m<<"] = "<<multReal<<"(w[j*"<<(2*size)<<"/"<<(3*L)<<"], d1r-d2r, d1i-d2i);\n";
            source<<"imag"<<output<<"[i+(2*j+2)*"<<m<<"] = "<<multImag<<"(w[j*"<<(2*size)<<"/"<<(3*L)<<"], d1r-d2r, d1i-d2i);\n";
        }
        else if (radix == 2) {
            source<<"mixed3 c0r = real"<<input<<"[i];\n";
            source<<"mixed3 c0i = imag"<<input<<"[i];\n";
            source<<"mixed3 c1r = real"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"mixed3 c1i = imag"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"real"<<output<<"[i+j*"<<m<<"] = c0r+c1r;\n";
            source<<"imag"<<output<<"[i+j*"<<m<<"] = c0i+c1i;\n";
            source<<"real"<<output<<"[i+(j+1)*"<<m<<"] = "<<multReal<<"(w[j*"<<size<<"/"<<(2*L)<<"], c0r-c1r, c0i-c1i);\n";
            source<<"imag"<<output<<"[i+(j+1)*"<<m<<"] = "<<multImag<<"(w[j*"<<size<<"/"<<(2*L)<<"], c0r-c1r, c0i-c1i);\n";
        }
        source<<"}\n";
        m = m*radix;
        source<<"SYNC_THREADS;\n";
        source<<"}\n";
        ++stage;
    }

    // Create the kernel.

    if (stage%2 == 1) {
        source<<"real0[indexInBlock] = real1[indexInBlock];\n";
        source<<"imag0[indexInBlock] = imag1[indexInBlock];\n";
        source<<"SYNC_THREADS;\n";
    }
    source<<"}\n";
    return source.str();
}

void CommonIntegrateRPMDStepKernel::downloadPositionsFromGPU(int beadIndex, std::vector<Vec3>& outPositions) {
    // positions array layout: stored as mm_float4 or mm_double4 depending on precision
    // [bead0_atom0, bead0_atom1, ..., bead1_atom0, ...]
    outPositions.resize(numParticles);
    
    int paddedParticles = cc.getPaddedNumAtoms();
    bool useDoublePrecision = (cc.getUseDoublePrecision() || cc.getUseMixedPrecision());
    
    // Download all positions from GPU
    if (useDoublePrecision) {
        vector<mm_double4> posData(paddedParticles * numCopies);
        positions.download(posData);  // Pass vector directly, not .data()
        
        // Extract this bead's slice (only non-padded particles)
        int offset = beadIndex * paddedParticles;
        for (int i = 0; i < numParticles; i++) {
            outPositions[i] = Vec3(posData[offset + i].x, posData[offset + i].y, posData[offset + i].z);
        }
    } else {
        vector<mm_float4> posData(paddedParticles * numCopies);
        positions.download(posData);  // Pass vector directly, not .data()
        
        // Extract this bead's slice (only non-padded particles)
        int offset = beadIndex * paddedParticles;
        for (int i = 0; i < numParticles; i++) {
            outPositions[i] = Vec3(posData[offset + i].x, posData[offset + i].y, posData[offset + i].z);
        }
    }
}

void CommonIntegrateRPMDStepKernel::uploadForcesToGPU(int beadIndex, const std::vector<Vec3>& inForces) {
    // This function is deprecated - use uploadAllForcesToGPU instead
    throw OpenMMException("uploadForcesToGPU should not be called - use uploadAllForcesToGPU");
}

void CommonIntegrateRPMDStepKernel::uploadAllForcesToGPU(const std::vector<std::vector<Vec3>>& allBeadForces) {
    // Upload all bead forces at once to avoid multiple download/upload cycles
    // forces array layout: [bead0_atom0_x, bead0_atom0_y, bead0_atom0_z, bead0_atom1_x, ...]
    // forces is stored as long long (fixed point for atomics)
    
    int paddedParticles = cc.getPaddedNumAtoms();
    int totalSize = paddedParticles * numCopies * 3;
    vector<long long> forceData(totalSize, 0);  // Initialize to zero (important for padded atoms)
    
    double forceScale = (cc.getUseDoublePrecision() || cc.getUseMixedPrecision() ? 1.0 : 0x100000000);
    
    // Pack all bead forces into flat array
    for (int b = 0; b < numCopies; b++) {
        if (allBeadForces[b].size() != (size_t)numParticles) {
            throw OpenMMException("uploadAllForcesToGPU: force vector size mismatch for bead " + std::to_string(b));
        }
        
        int offset = b * paddedParticles * 3;  // Use paddedParticles, not numParticles
        for (int i = 0; i < numParticles; i++) {
            forceData[offset + i*3 + 0] = (long long)(allBeadForces[b][i][0] * forceScale);
            forceData[offset + i*3 + 1] = (long long)(allBeadForces[b][i][1] * forceScale);
            forceData[offset + i*3 + 2] = (long long)(allBeadForces[b][i][2] * forceScale);
        }
        // Padded atoms (from numParticles to paddedParticles) remain zero
    }
    
    // Upload to GPU in one shot - pass vector directly, not .data()
    forces.upload(forceData);
}

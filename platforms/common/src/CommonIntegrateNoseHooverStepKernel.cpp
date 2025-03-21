/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
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

#include "openmm/common/CommonIntegrateNoseHooverStepKernel.h"
#include "openmm/common/CommonKernelUtilities.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "CommonKernelSources.h"
#include "SimTKOpenMMRealType.h"

using namespace OpenMM;
using namespace std;

void CommonIntegrateNoseHooverStepKernel::initialize(const System& system, const NoseHooverIntegrator& integrator) {
    cc.initializeContexts();
    ContextSelector selector(cc);
    bool useDouble = cc.getUseDoublePrecision() || cc.getUseMixedPrecision();
    map<string, string> defines;
    defines["BOLTZ"] = cc.doubleToString(BOLTZ, true);
    ComputeProgram program = cc.compileProgram(CommonKernelSources::noseHooverIntegrator, defines);
    kernel1 = program->createKernel("integrateNoseHooverMiddlePart1");
    kernel2 = program->createKernel("integrateNoseHooverMiddlePart2");
    kernel3 = program->createKernel("integrateNoseHooverMiddlePart3");
    kernel4 = program->createKernel("integrateNoseHooverMiddlePart4");
    if (useDouble) {
        oldDelta.initialize<mm_double4>(cc, cc.getPaddedNumAtoms(), "oldDelta");
    } else {
        oldDelta.initialize<mm_float4>(cc, cc.getPaddedNumAtoms(), "oldDelta");
    }
    kernelHardWall = program->createKernel("integrateNoseHooverHardWall");
    prevMaxPairDistance = -1.0f;
    maxPairDistanceBuffer.initialize<float>(cc, 1, "maxPairDistanceBuffer");

    int workGroupSize = std::min(cc.getMaxThreadBlockSize(), 512);
    defines["WORK_GROUP_SIZE"] = std::to_string(workGroupSize);

    defines["BEGIN_YS_LOOP"] = "const real arr[1] = {1.0};"
                               "for(int i=0;i<1;++i) {"
                               "const real ys = arr[i];";
    defines["END_YS_LOOP"] = "}";
    program = cc.compileProgram(CommonKernelSources::noseHooverChain, defines);
    propagateKernels[1] = program->createKernel("propagateNoseHooverChain");

    defines["BEGIN_YS_LOOP"] = "const real arr[3] = {0.828981543588751, -0.657963087177502, 0.828981543588751};"
                               "for(int i=0;i<3;++i) {"
                               "const real ys = arr[i];";
    program = cc.compileProgram(CommonKernelSources::noseHooverChain, defines);
    propagateKernels[3] = program->createKernel("propagateNoseHooverChain");

    defines["BEGIN_YS_LOOP"] = "const real arr[5] = {0.2967324292201065, 0.2967324292201065, -0.186929716880426, 0.2967324292201065, 0.2967324292201065};"
                               "for(int i=0;i<5;++i) {"
                               "const real ys = arr[i];";
    program = cc.compileProgram(CommonKernelSources::noseHooverChain, defines);
    propagateKernels[5] = program->createKernel("propagateNoseHooverChain");

    defines["BEGIN_YS_LOOP"] = "const real arr[7] = {0.784513610477560, 0.235573213359357, -1.17767998417887, 1.31518632068391,-1.17767998417887, 0.235573213359357, 0.784513610477560};"
                               "for(int i=0;i<7;++i) {"
                               "const real ys = arr[i];";
    program = cc.compileProgram(CommonKernelSources::noseHooverChain, defines);
    propagateKernels[7] = program->createKernel("propagateNoseHooverChain");
    program = cc.compileProgram(CommonKernelSources::noseHooverChain, defines);
    reduceEnergyKernel = program->createKernel("reduceEnergyPair");

    computeHeatBathEnergyKernel = program->createKernel("computeHeatBathEnergy");
    computeAtomsKineticEnergyKernel = program->createKernel("computeAtomsKineticEnergy");
    computePairsKineticEnergyKernel = program->createKernel("computePairsKineticEnergy");
    scaleAtomsVelocitiesKernel = program->createKernel("scaleAtomsVelocities");
    scalePairsVelocitiesKernel = program->createKernel("scalePairsVelocities");
    int energyBufferSize = cc.getEnergyBuffer().getSize();
    if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision())
        energyBuffer.initialize<mm_double2>(cc, energyBufferSize, "energyBuffer");
    else
        energyBuffer.initialize<mm_float2>(cc, energyBufferSize, "energyBuffer");
}

void CommonIntegrateNoseHooverStepKernel::execute(ContextImpl& context, const NoseHooverIntegrator& integrator) {
    ContextSelector selector(cc);
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    double dt = integrator.getStepSize();
    cc.getIntegrationUtilities().setNextStepSize(dt);

    // If the atom reordering has occured, the forces from the previous step are permuted and thus invalid.
    // They need to be either sorted or recomputed; here we choose the latter.
    if (cc.getAtomsWereReordered())
        context.calcForcesAndEnergy(true, false, integrator.getIntegrationForceGroups());

    const auto& atomList = integrator.getAllThermostatedIndividualParticles();
    const auto& pairList = integrator.getAllThermostatedPairs();
    int numAtoms = atomList.size();
    int numPairs = pairList.size();
    int numParticles = numAtoms + 2*numPairs;
    float maxPairDistance = integrator.getMaximumPairDistance();
    // Make sure atom and pair metadata is uploaded and has the correct dimensions
    if (prevMaxPairDistance != maxPairDistance) {
        std::vector<float> tmp(1, maxPairDistance);
        maxPairDistanceBuffer.upload(tmp);
        prevMaxPairDistance = maxPairDistance;
    }
    if (numAtoms !=0 && (!atomListBuffer.isInitialized() || atomListBuffer.getSize() != numAtoms)) {
        if (atomListBuffer.isInitialized())
            atomListBuffer.resize(atomList.size());
        else
            atomListBuffer.initialize<int>(cc, atomList.size(), "atomListBuffer");
        atomListBuffer.upload(atomList);
    }
    if (numPairs !=0 && (!pairListBuffer.isInitialized() || pairListBuffer.getSize() != numPairs)) {
        if (pairListBuffer.isInitialized()) {
            pairListBuffer.resize(pairList.size());
            pairTemperatureBuffer.resize(pairList.size());
        }
        else {
            pairListBuffer.initialize<mm_int2>(cc, pairList.size(), "pairListBuffer");
            pairTemperatureBuffer.initialize<float>(cc, pairList.size(), "pairTemperatureBuffer");
        }
        std::vector<mm_int2> tmp;
        std::vector<float> tmp2;
        for(const auto &pair : pairList) {
            tmp.push_back(mm_int2(std::get<0>(pair), std::get<1>(pair)));
            tmp2.push_back(std::get<2>(pair));
        }
        pairListBuffer.upload(tmp);
        pairTemperatureBuffer.upload(tmp2);
    }
    int totalAtoms = cc.getNumAtoms();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1->addArg(numAtoms);
        kernel1->addArg(numPairs);
        kernel1->addArg(paddedNumAtoms);
        kernel1->addArg(cc.getVelm());
        kernel1->addArg(cc.getLongForceBuffer());
        kernel1->addArg(integration.getStepSize());
        kernel1->addArg(numAtoms > 0 ? atomListBuffer : cc.getEnergyBuffer()); // The array is not used if num == 0
        kernel1->addArg(numPairs > 0 ? pairListBuffer : cc.getEnergyBuffer()); // The array is not used if num == 0
        kernel2->addArg(totalAtoms);
        kernel2->addArg(cc.getVelm());
        kernel2->addArg(integration.getPosDelta());
        kernel2->addArg(oldDelta);
        kernel2->addArg(integration.getStepSize());
        kernel3->addArg(totalAtoms);
        kernel3->addArg(cc.getVelm());
        kernel3->addArg(integration.getPosDelta());
        kernel3->addArg(oldDelta);
        kernel3->addArg(integration.getStepSize());
        kernel4->addArg(totalAtoms);
        kernel4->addArg(cc.getPosq());
        kernel4->addArg(cc.getVelm());
        kernel4->addArg(integration.getPosDelta());
        kernel4->addArg(oldDelta);
        kernel4->addArg(integration.getStepSize());
        if (cc.getUseMixedPrecision())
            kernel4->addArg(cc.getPosqCorrection());
        if (numPairs > 0) {
            kernelHardWall->addArg(numPairs);
            kernelHardWall->addArg(maxPairDistanceBuffer);
            kernelHardWall->addArg(integration.getStepSize());
            kernelHardWall->addArg(cc.getPosq());
            kernelHardWall->addArg(cc.getVelm());
            kernelHardWall->addArg(pairListBuffer);
            kernelHardWall->addArg(pairTemperatureBuffer);
            if (cc.getUseMixedPrecision())
                kernelHardWall->addArg(cc.getPosqCorrection());
        }
    }

    /*
     * Carry out the LF-middle integration (c.f. J. Phys. Chem. A 2019, 123, 6056−6079)
     */
    // Velocity update
    kernel1->execute(std::max(numAtoms, numPairs));
    integration.applyVelocityConstraints(integrator.getConstraintTolerance());
    // Position update
    kernel2->execute(numParticles);
    // Apply the thermostat
    int numChains = integrator.getNumThermostats();
    for(int chain = 0; chain < numChains; ++chain) {
        const auto &thermostatChain = integrator.getThermostat(chain);
        auto KEs = computeMaskedKineticEnergy(context, thermostatChain, false);
        auto scaleFactors = propagateChain(context, thermostatChain, KEs, dt);
        scaleVelocities(context, thermostatChain, scaleFactors);
    }
    // Position update
    kernel3->execute(numParticles);
    integration.applyConstraints(integrator.getConstraintTolerance());
    // Apply constraint forces
    kernel4->execute(numAtoms);
    // Make sure any Drude-like particles have not wandered too far from home
    if (numPairs > 0) kernelHardWall->execute(numPairs);
    integration.computeVirtualSites();

    // Update the time and step count.
    cc.setTime(cc.getTime()+dt);
    cc.setStepCount(cc.getStepCount()+1);
    cc.reorderAtoms();

    // Reduce UI lag.

    flushPeriodically(cc);
}

double CommonIntegrateNoseHooverStepKernel::computeKineticEnergy(ContextImpl& context, const NoseHooverIntegrator& integrator) {
    return cc.getIntegrationUtilities().computeKineticEnergy(0);
}


std::pair<double, double> CommonIntegrateNoseHooverStepKernel::propagateChain(ContextImpl& context, const NoseHooverChain &nhc, std::pair<double, double> kineticEnergies, double timeStep) {
    bool useDouble = cc.getUseDoublePrecision() || cc.getUseMixedPrecision();
    int chainID = nhc.getChainID();
    int nAtoms = nhc.getThermostatedAtoms().size();
    int nPairs = nhc.getThermostatedPairs().size();
    int chainLength = nhc.getChainLength();
    int numYS = nhc.getNumYoshidaSuzukiTimeSteps();
    int numMTS = nhc.getNumMultiTimeSteps();

    if (numYS != 1 && numYS != 3 && numYS != 5 && numYS != 7) {
        throw OpenMMException("Number of Yoshida Suzuki time steps has to be 1, 3, 5, or 7.");
    }

    if (!scaleFactorBuffer.isInitialized() || scaleFactorBuffer.getSize() == 0) {
        if (useDouble) {
            std::vector<mm_double2> zeros{{0,0}};
            if (scaleFactorBuffer.isInitialized())
                scaleFactorBuffer.resize(1);
            else
                scaleFactorBuffer.initialize<mm_double2>(cc, 1, "scaleFactorBuffer");
            scaleFactorBuffer.upload(zeros);
        }
        else {
            std::vector<mm_float2> zeros{{0,0}};
            if (scaleFactorBuffer.isInitialized())
                scaleFactorBuffer.resize(1);
            else
                scaleFactorBuffer.initialize<mm_float2>(cc, 1, "scaleFactorBuffer");
            scaleFactorBuffer.upload(zeros);
        }
    }
    if (!chainForces.isInitialized() || !chainMasses.isInitialized()) {
        if (useDouble) {
            std::vector<double> zeros(chainLength,0);
            if (chainForces.isInitialized()) {
                chainMasses.resize(chainLength);
                chainForces.resize(chainLength);
            }
            else {
                chainMasses.initialize<double>(cc, chainLength, "chainMasses");
                chainForces.initialize<double>(cc, chainLength, "chainForces");
            }
            chainMasses.upload(zeros);
            chainForces.upload(zeros);
        }
        else {
            std::vector<float> zeros(chainLength,0);
            if (chainForces.isInitialized()) {
                chainMasses.resize(chainLength);
                chainForces.resize(chainLength);
            }
            else {
                chainMasses.initialize<float>(cc, chainLength, "chainMasses");
                chainForces.initialize<float>(cc, chainLength, "chainForces");
            }
            chainMasses.upload(zeros);
            chainForces.upload(zeros);
        }
    }
    if (chainForces.getSize() < chainLength)
        chainForces.resize(chainLength);
    if (chainMasses.getSize() < chainLength)
        chainMasses.resize(chainLength);


    // N.B. We ignore the incoming kineticEnergy and grab it from the device buffer instead
    if (nAtoms) {
        if (!chainState.count(2*chainID))
            chainState[2*chainID] = ComputeArray();
        if (!chainState.at(2*chainID).isInitialized() || chainState.at(2*chainID).getSize() != chainLength) {
            // We need to upload the Common array
            if (useDouble) {
                if (chainState.at(2*chainID).isInitialized())
                    chainState.at(2*chainID).resize(chainLength);
                else
                    chainState.at(2*chainID).initialize<mm_double2>(cc, chainLength, "chainState" + std::to_string(2*chainID));
                std::vector<mm_double2> zeros(chainLength, mm_double2(0.0, 0.0));
                chainState.at(2*chainID).upload(zeros.data());
            }
            else {
                if (chainState.at(2*chainID).isInitialized())
                    chainState.at(2*chainID).resize(chainLength);
                else
                    chainState.at(2*chainID).initialize<mm_float2>(cc, chainLength, "chainState" + std::to_string(2*chainID));
                std::vector<mm_float2> zeros(chainLength, mm_float2(0.0f, 0.0f));
                chainState.at(2*chainID).upload(zeros.data());
            }
        }
    }

    if (nPairs) {
        if (!chainState.count(2*chainID+1))
            chainState[2*chainID+1] = ComputeArray();
        if (!chainState.at(2*chainID+1).isInitialized() || chainState.at(2*chainID+1).getSize() != chainLength) {
            // We need to upload the Common array
            if (useDouble) {
                if (chainState.at(2*chainID+1).isInitialized())
                    chainState.at(2*chainID+1).resize(chainLength);
                else
                    chainState.at(2*chainID+1).initialize<mm_double2>(cc, chainLength, "chainState" + std::to_string(2*chainID+1));
                std::vector<mm_double2> zeros(chainLength, mm_double2(0.0, 0.0));
                chainState.at(2*chainID+1).upload(zeros.data());
            }
            else {
                if (chainState.at(2*chainID+1).isInitialized())
                    chainState.at(2*chainID+1).resize(chainLength);
                else
                    chainState.at(2*chainID+1).initialize<mm_float2>(cc, chainLength, "chainState" + std::to_string(2*chainID+1));
                std::vector<mm_float2> zeros(chainLength, mm_float2(0.0f, 0.0f));
                chainState.at(2*chainID+1).upload(zeros.data());
            }
        }
    }

    if (!hasInitializedPropagateKernel) {
        hasInitializedPropagateKernel = true;
        propagateKernels[numYS]->addArg(); // ChainState
        propagateKernels[numYS]->addArg(kineticEnergyBuffer);
        propagateKernels[numYS]->addArg(scaleFactorBuffer);
        propagateKernels[numYS]->addArg(chainMasses);
        propagateKernels[numYS]->addArg(chainForces);
        propagateKernels[numYS]->addArg(); // ChainType
        propagateKernels[numYS]->addArg(chainLength);
        propagateKernels[numYS]->addArg(numMTS);
        propagateKernels[numYS]->addArg(); // numDoFs
        propagateKernels[numYS]->addArg((float)timeStep);
        propagateKernels[numYS]->addArg(); // kT
        propagateKernels[numYS]->addArg(); // frequency
    }

    if (nAtoms) {
        int chainType = 0;
        double temperature = nhc.getTemperature();
        float frequency = nhc.getCollisionFrequency();
        double kT = BOLTZ * temperature;
        int numDOFs = nhc.getNumDegreesOfFreedom();
        propagateKernels[numYS]->setArg(0, chainState[2*chainID]);
        propagateKernels[numYS]->setArg(5, chainType);
        propagateKernels[numYS]->setArg(8, numDOFs);
        if (useDouble) {
            propagateKernels[numYS]->setArg(10, kT);
        } else {
            propagateKernels[numYS]->setArg(10, (float)kT);
        }
        propagateKernels[numYS]->setArg(11, frequency);
        propagateKernels[numYS]->execute(1, 1);
    }
    if (nPairs) {
        int chainType = 1;
        double relativeTemperature = nhc.getRelativeTemperature();
        float relativeFrequency = nhc.getRelativeCollisionFrequency();
        double kT = BOLTZ * relativeTemperature;
        int ndf = 3*nPairs;
        propagateKernels[numYS]->setArg(0, chainState[2*chainID+1]);
        propagateKernels[numYS]->setArg(5, chainType);
        propagateKernels[numYS]->setArg(8, ndf);
        if (useDouble) {
            propagateKernels[numYS]->setArg(10, kT);
        } else {
            propagateKernels[numYS]->setArg(10, (float)kT);
        }
        propagateKernels[numYS]->setArg(11, relativeFrequency);
        propagateKernels[numYS]->execute(1, 1);
    }
    return {0, 0};
}

double CommonIntegrateNoseHooverStepKernel::computeHeatBathEnergy(ContextImpl& context, const NoseHooverChain &nhc) {
    ContextSelector selector(cc);
    bool useDouble = cc.getUseDoublePrecision() || cc.getUseMixedPrecision();
    int chainID = nhc.getChainID();
    int chainLength = nhc.getChainLength();

    bool absChainIsValid = chainState.count(2*chainID) != 0 &&
                           chainState[2*chainID].isInitialized() &&
                           chainState[2*chainID].getSize() == chainLength;
    bool relChainIsValid = chainState.count(2*chainID+1) != 0 &&
                           chainState[2*chainID+1].isInitialized() &&
                           chainState[2*chainID+1].getSize() == chainLength;

    if (!absChainIsValid && !relChainIsValid) return 0.0;

    if (!heatBathEnergy.isInitialized() || heatBathEnergy.getSize() == 0) {
        if (useDouble) {
            std::vector<double> one(1);
            heatBathEnergy.initialize<double>(cc, 1, "heatBathEnergy");
            heatBathEnergy.upload(one);
        }
        else {
            std::vector<float> one(1);
            heatBathEnergy.initialize<float>(cc, 1, "heatBathEnergy");
            heatBathEnergy.upload(one);
        }
    }

    cc.clearBuffer(heatBathEnergy);

    if(!hasInitializedHeatBathEnergyKernel) {
        hasInitializedHeatBathEnergyKernel = true;
        computeHeatBathEnergyKernel->addArg(heatBathEnergy);
        computeHeatBathEnergyKernel->addArg(chainLength);
        computeHeatBathEnergyKernel->addArg(); // numDOFs
        computeHeatBathEnergyKernel->addArg(); // kT
        computeHeatBathEnergyKernel->addArg(); // frequency
        computeHeatBathEnergyKernel->addArg(); // chainstate
    }

    if (absChainIsValid) {
        int numDOFs = nhc.getNumDegreesOfFreedom();
        double temperature = nhc.getTemperature();
        float frequency = nhc.getCollisionFrequency();
        double kT = BOLTZ * temperature;

        computeHeatBathEnergyKernel->setArg(2, numDOFs);
        if (useDouble) {
            computeHeatBathEnergyKernel->setArg(3, kT);
        } else {
            computeHeatBathEnergyKernel->setArg(3, (float)kT);
        }
        computeHeatBathEnergyKernel->setArg(4, frequency);
        computeHeatBathEnergyKernel->setArg(5, chainState[2*chainID]);
        computeHeatBathEnergyKernel->execute(1, 1);
    }
    if (relChainIsValid) {
        int numDOFs = 3 * nhc.getThermostatedPairs().size();
        double temperature = nhc.getRelativeTemperature();
        float frequency = nhc.getRelativeCollisionFrequency();
        double kT = BOLTZ * temperature;

        computeHeatBathEnergyKernel->setArg(2, numDOFs);
        if (useDouble) {
            computeHeatBathEnergyKernel->setArg(3, kT);
        } else {
            computeHeatBathEnergyKernel->setArg(3, (float)kT);
        }
        computeHeatBathEnergyKernel->setArg(4, frequency);
        computeHeatBathEnergyKernel->setArg(5, chainState[2*chainID+1]);
        computeHeatBathEnergyKernel->execute(1, 1);
    }


    void * pinnedBuffer = cc.getPinnedBuffer();
    heatBathEnergy.download(pinnedBuffer);
    if (useDouble)
        return *((double*) pinnedBuffer);
    else
        return *((float*) pinnedBuffer);
}

std::pair<double, double> CommonIntegrateNoseHooverStepKernel::computeMaskedKineticEnergy(ContextImpl& context, const NoseHooverChain &nhc, bool downloadValue) {
    ContextSelector selector(cc);
    bool useDouble = cc.getUseDoublePrecision() || cc.getUseMixedPrecision();
    int chainID = nhc.getChainID();
    const auto & nhcAtoms = nhc.getThermostatedAtoms();
    const auto & nhcPairs = nhc.getThermostatedPairs();
    int nAtoms = nhcAtoms.size();
    int nPairs = nhcPairs.size();
    if (nAtoms) {
        if (!atomlists.count(chainID)) {
            // We need to upload the Common array
            atomlists[chainID] = ComputeArray();
            atomlists[chainID].initialize<int>(cc, nAtoms, "atomlist" + std::to_string(chainID));
            atomlists[chainID].upload(nhcAtoms);
        }
        if (atomlists[chainID].getSize() != nAtoms) {
            throw OpenMMException("Number of atoms changed. Cannot be handled by the same Nose-Hoover thermostat.");
        }
    }
    if (nPairs) {
        if (!pairlists.count(chainID)) {
            // We need to upload the Common array
            pairlists[chainID] = ComputeArray();
            pairlists[chainID].initialize<mm_int2>(cc, nPairs, "pairlist" + std::to_string(chainID));
            std::vector<mm_int2> int2vec;
            for(const auto &p : nhcPairs) int2vec.push_back(mm_int2(p.first, p.second));
            pairlists[chainID].upload(int2vec);
        }
        if (pairlists[chainID].getSize() != nPairs) {
            throw OpenMMException("Number of thermostated pairs changed. Cannot be handled by the same Nose-Hoover thermostat.");
        }
    }
    if (!kineticEnergyBuffer.isInitialized() || kineticEnergyBuffer.getSize() == 0) {
        if (useDouble) {
            std::vector<mm_double2> zeros{{0,0}};
            kineticEnergyBuffer.initialize<mm_double2>(cc, 1, "kineticEnergyBuffer");
            kineticEnergyBuffer.upload(zeros);
        }
        else {
            std::vector<mm_float2> zeros{{0,0}};
            kineticEnergyBuffer.initialize<mm_float2>(cc, 1, "kineticEnergyBuffer");
            kineticEnergyBuffer.upload(zeros);
        }
    }

    int workGroupSize = std::min(cc.getMaxThreadBlockSize(), 512);
    if (!hasInitializedKineticEnergyKernel) {
        hasInitializedKineticEnergyKernel = true;
        computeAtomsKineticEnergyKernel->addArg(energyBuffer);
        computeAtomsKineticEnergyKernel->addArg(); // nAtoms
        computeAtomsKineticEnergyKernel->addArg(cc.getVelm());
        computeAtomsKineticEnergyKernel->addArg(); // atom list

        computePairsKineticEnergyKernel->addArg(energyBuffer);
        computePairsKineticEnergyKernel->addArg(); // nPairs
        computePairsKineticEnergyKernel->addArg(cc.getVelm());
        computePairsKineticEnergyKernel->addArg(); // pair list

        reduceEnergyKernel->addArg(energyBuffer);
        reduceEnergyKernel->addArg(kineticEnergyBuffer);
        reduceEnergyKernel->addArg((int) energyBuffer.getSize());
    }

    cc.clearBuffer(energyBuffer);
    if (nAtoms) {
        computeAtomsKineticEnergyKernel->setArg(1, nAtoms);
        computeAtomsKineticEnergyKernel->setArg(3, atomlists[chainID]);
        computeAtomsKineticEnergyKernel->execute(nAtoms);
    }
    if (nPairs) {
        computePairsKineticEnergyKernel->setArg(1, nPairs);
        computePairsKineticEnergyKernel->setArg(3, pairlists[chainID]);
        computePairsKineticEnergyKernel->execute(nPairs);
    }
    reduceEnergyKernel->execute(workGroupSize, workGroupSize);

    std::pair<double, double> KEs = {0, 0};
    if (downloadValue) {
        if (useDouble) {
            mm_double2 tmp;
            kineticEnergyBuffer.download(&tmp);
            KEs.first = tmp.x;
            KEs.second = tmp.y;
        }
        else {
            mm_float2 tmp;
            kineticEnergyBuffer.download(&tmp);
            KEs.first = tmp.x;
            KEs.second = tmp.y;
        }
    }
    return KEs;
}

void CommonIntegrateNoseHooverStepKernel::scaleVelocities(ContextImpl& context, const NoseHooverChain &nhc, std::pair<double, double> scaleFactor) {
    // For now we assume that the atoms and pairs info is valid, because compute{Atoms|Pairs}KineticEnergy must have been
    // called before this kernel.  If that ever ceases to be true, some sanity checks are needed here.

    int chainID = nhc.getChainID();
    int nAtoms = nhc.getThermostatedAtoms().size();
    int nPairs = nhc.getThermostatedPairs().size();
    if (!hasInitializedScaleVelocitiesKernel) {
        hasInitializedScaleVelocitiesKernel = true;
        scaleAtomsVelocitiesKernel->addArg(scaleFactorBuffer);
        scaleAtomsVelocitiesKernel->addArg(); // nAtoms
        scaleAtomsVelocitiesKernel->addArg(cc.getVelm());
        scaleAtomsVelocitiesKernel->addArg(); // atom list

        scalePairsVelocitiesKernel->addArg(scaleFactorBuffer);
        scalePairsVelocitiesKernel->addArg(); // nPairs
        scalePairsVelocitiesKernel->addArg(cc.getVelm());
        scalePairsVelocitiesKernel->addArg(); // pair list
    }
    if (nAtoms) {
        scaleAtomsVelocitiesKernel->setArg(1, nAtoms);
        scaleAtomsVelocitiesKernel->setArg(3, atomlists[chainID]);
        scaleAtomsVelocitiesKernel->execute(nAtoms);
    }
    if (nPairs) {
        scalePairsVelocitiesKernel->setArg(1, nPairs);
        scalePairsVelocitiesKernel->setArg(3, pairlists[chainID]);
        scalePairsVelocitiesKernel->execute(nPairs);
    }
}

void CommonIntegrateNoseHooverStepKernel::createCheckpoint(ContextImpl& context, ostream& stream) const {
    ContextSelector selector(cc);
    int numChains = chainState.size();
    bool useDouble = cc.getUseDoublePrecision() || cc.getUseMixedPrecision();
    stream.write((char*) &numChains, sizeof(int));
    for (auto& state : chainState){
        int chainID = state.first;
        int chainLength = state.second.getSize();
        stream.write((char*) &chainID, sizeof(int));
        stream.write((char*) &chainLength, sizeof(int));
        if (useDouble) {
            vector<mm_double2> stateVec;
            state.second.download(stateVec);
            stream.write((char*) stateVec.data(), sizeof(mm_double2)*chainLength);
        }
        else {
            vector<mm_float2> stateVec;
            state.second.download(stateVec);
            stream.write((char*) stateVec.data(), sizeof(mm_float2)*chainLength);
        }
    }
}

void CommonIntegrateNoseHooverStepKernel::loadCheckpoint(ContextImpl& context, istream& stream) {
    ContextSelector selector(cc);
    int numChains;
    bool useDouble = cc.getUseDoublePrecision() || cc.getUseMixedPrecision();
    stream.read((char*) &numChains, sizeof(int));
    chainState.clear();
    for (int i = 0; i < numChains; i++) {
        int chainID, chainLength;
        stream.read((char*) &chainID, sizeof(int));
        stream.read((char*) &chainLength, sizeof(int));
        if (useDouble) {
            chainState[chainID] = ComputeArray();
            chainState[chainID].initialize<mm_double2>(cc, chainLength, "chainState" + to_string(chainID));
            vector<mm_double2> stateVec(chainLength);
            stream.read((char*) &stateVec[0], sizeof(mm_double2)*chainLength);
            chainState[chainID].upload(stateVec);
        }
        else {
            chainState[chainID] = ComputeArray();
            chainState[chainID].initialize<mm_float2>(cc, chainLength, "chainState" + to_string(chainID));
            vector<mm_float2> stateVec(chainLength);
            stream.read((char*) &stateVec[0], sizeof(mm_float2)*chainLength);
            chainState[chainID].upload(stateVec);
        }
    }
}

void CommonIntegrateNoseHooverStepKernel::getChainStates(ContextImpl& context, vector<vector<double> >& positions, vector<vector<double> >& velocities) const {
    ContextSelector selector(cc);
    int numChains = chainState.size();
    bool useDouble = cc.getUseDoublePrecision() || cc.getUseMixedPrecision();
    positions.clear();
    velocities.clear();
    positions.resize(numChains);
    velocities.resize(numChains);
    for (int i = 0; i < numChains; i++) {
        const ComputeArray& state = chainState.at(i);
        if (useDouble) {
            vector<mm_double2> stateVec;
            state.download(stateVec);
            for (int j = 0; j < stateVec.size(); j++) {
                positions[i].push_back(stateVec[j].x);
                velocities[i].push_back(stateVec[j].y);
            }
        }
        else {
            vector<mm_float2> stateVec;
            state.download(stateVec);
            for (int j = 0; j < stateVec.size(); j++) {
                positions[i].push_back((float) stateVec[j].x);
                velocities[i].push_back((float) stateVec[j].y);
            }
        }
    }
}

void CommonIntegrateNoseHooverStepKernel::setChainStates(ContextImpl& context, const vector<vector<double> >& positions, const vector<vector<double> >& velocities) {
    ContextSelector selector(cc);
    int numChains = positions.size();
    bool useDouble = cc.getUseDoublePrecision() || cc.getUseMixedPrecision();
    chainState.clear();
    for (int i = 0; i < numChains; i++) {
        int chainLength = positions[i].size();
        chainState[i] = ComputeArray();
        if (useDouble) {
            chainState[i].initialize<mm_double2>(cc, chainLength, "chainState"+cc.intToString(i));
            vector<mm_double2> stateVec;
            for (int j = 0; j < chainLength; j++)
                stateVec.push_back(mm_double2(positions[i][j], velocities[i][j]));
            chainState[i].upload(stateVec);
        }
        else {
            chainState[i].initialize<mm_float2>(cc, chainLength, "chainState"+cc.intToString(i));
            vector<mm_float2> stateVec;
            for (int j = 0; j < chainLength; j++)
                stateVec.push_back(mm_float2((float) positions[i][j], (float) velocities[i][j]));
            chainState[i].upload(stateVec);
        }
    }
}

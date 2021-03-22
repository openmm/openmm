/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
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
#include "openmm/common/IntegrationUtilities.h"
#include "openmm/common/ExpressionUtilities.h"
#include "openmm/common/NonbondedUtilities.h"
#include "SimTKOpenMMRealType.h"

using namespace OpenMM;
using namespace std;


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
    numCopies = integrator.getNumCopies();
    numParticles = system.getNumParticles();
    workgroupSize = numCopies;
    if (numCopies != findFFTDimension(numCopies))
        throw OpenMMException("RPMDIntegrator: the number of copies must be a multiple of powers of 2, 3, and 5.");
    int paddedParticles = cc.getPaddedNumAtoms();
    bool useDoublePrecision = (cc.getUseDoublePrecision() || cc.getUseMixedPrecision());
    int elementSize = (useDoublePrecision ? sizeof(mm_double4) : sizeof(mm_float4));
    forces.initialize<long long>(cc, numCopies*paddedParticles*3, "rpmdForces");
    positions.initialize(cc, numCopies*paddedParticles, elementSize, "rpmdPositions");
    velocities.initialize(cc, numCopies*paddedParticles, elementSize, "rpmdVelocities");
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
    defines["HBAR"] = cc.doubleToString(1.054571628e-34*AVOGADRO/(1000*1e-12));
    defines["SCALE"] = cc.doubleToString(1.0/sqrt((double) numCopies));
    defines["M_PI"] = cc.doubleToString(M_PI);
    map<string, string> replacements;
    replacements["FFT_Q_FORWARD"] = createFFT(numCopies, "q", true);
    replacements["FFT_Q_BACKWARD"] = createFFT(numCopies, "q", false);
    replacements["FFT_V_FORWARD"] = createFFT(numCopies, "v", true);
    replacements["FFT_V_BACKWARD"] = createFFT(numCopies, "v", false);
    ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonRpmdKernelSources::rpmd, replacements), defines);
    pileKernel = program->createKernel("applyPileThermostat");
    stepKernel = program->createKernel("integrateStep");
    velocitiesKernel = program->createKernel("advanceVelocities");
    copyToContextKernel = program->createKernel("copyDataToContext");
    copyFromContextKernel = program->createKernel("copyDataFromContext");
    translateKernel = program->createKernel("applyCellTranslations");
    
    // Create kernels for doing contractions.
    
    for (auto& g : groupsByCopies) {
        int copies = g.first;
        replacements.clear();
        replacements["NUM_CONTRACTED_COPIES"] = cc.intToString(copies);
        replacements["POS_SCALE"] = cc.doubleToString(1.0/numCopies);
        replacements["FORCE_SCALE"] = cc.doubleToString(0x100000000/(double) copies);
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
    for (auto& g : groupsByCopies) {
        int copies = g.first;
        positionContractionKernels[copies]->addArg(positions);
        positionContractionKernels[copies]->addArg(contractedPositions);
        forceContractionKernels[copies]->addArg(forces);
        forceContractionKernels[copies]->addArg(contractedForces);
    }
}

void CommonIntegrateRPMDStepKernel::execute(ContextImpl& context, const RPMDIntegrator& integrator, bool forcesAreValid) {
    cc.setAsCurrent();
    if (!hasInitializedKernels)
        initializeKernels(context);
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    
    // Loop over copies and compute the force on each one.
    
    if (!forcesAreValid)
        computeForces(context);
    
    // Apply the PILE-L thermostat.
    
    bool useDoublePrecision = (cc.getUseDoublePrecision() || cc.getUseMixedPrecision());
    double dt = integrator.getStepSize();
    pileKernel->setArg(2, integration.prepareRandomNumbers(numParticles*numCopies));
    if (useDoublePrecision) {
        pileKernel->setArg(3, dt);
        pileKernel->setArg(4, integrator.getTemperature()*BOLTZ);
        pileKernel->setArg(5, integrator.getFriction());
        stepKernel->setArg(3, dt);
        stepKernel->setArg(4, integrator.getTemperature()*BOLTZ);
        velocitiesKernel->setArg(2, dt);
    }
    else {
        pileKernel->setArg(3, (float) dt);
        pileKernel->setArg(4, (float) (integrator.getTemperature()*BOLTZ));
        pileKernel->setArg(5, (float) integrator.getFriction());
        stepKernel->setArg(3, (float) dt);
        stepKernel->setArg(4, (float) (integrator.getTemperature()*BOLTZ));
        velocitiesKernel->setArg(2, (float) dt);
    }
    if (integrator.getApplyThermostat())
        pileKernel->execute(numParticles*numCopies, workgroupSize);

    // Update positions and velocities.
    
    stepKernel->execute(numParticles*numCopies, workgroupSize);

    // Calculate forces based on the updated positions.
    
    computeForces(context);
    
    // Update velocities.

    velocitiesKernel->execute(numParticles*numCopies, workgroupSize);

    // Apply the PILE-L thermostat again.

    if (integrator.getApplyThermostat()) {
        pileKernel->setArg(2, integration.prepareRandomNumbers(numParticles*numCopies));
        pileKernel->execute(numParticles*numCopies, workgroupSize);
    }

    // Update the time and step count.

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
    // Compute forces from all groups that didn't have a specified contraction.

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

    // Record the positions.

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
            source<<"mixed3 d2r = "<<cc.doubleToString(sin(0.4*M_PI))<<"*(c1r-c4r);\n";
            source<<"mixed3 d2i = "<<cc.doubleToString(sin(0.4*M_PI))<<"*(c1i-c4i);\n";
            source<<"mixed3 d3r = "<<cc.doubleToString(sin(0.4*M_PI))<<"*(c2r-c3r);\n";
            source<<"mixed3 d3i = "<<cc.doubleToString(sin(0.4*M_PI))<<"*(c2i-c3i);\n";
            source<<"mixed3 d4r = d0r+d1r;\n";
            source<<"mixed3 d4i = d0i+d1i;\n";
            source<<"mixed3 d5r = "<<cc.doubleToString(0.25*sqrt(5.0))<<"*(d0r-d1r);\n";
            source<<"mixed3 d5i = "<<cc.doubleToString(0.25*sqrt(5.0))<<"*(d0i-d1i);\n";
            source<<"mixed3 d6r = c0r-0.25f*d4r;\n";
            source<<"mixed3 d6i = c0i-0.25f*d4i;\n";
            source<<"mixed3 d7r = d6r+d5r;\n";
            source<<"mixed3 d7i = d6i+d5i;\n";
            source<<"mixed3 d8r = d6r-d5r;\n";
            source<<"mixed3 d8i = d6i-d5i;\n";
            string coeff = cc.doubleToString(sin(0.2*M_PI)/sin(0.4*M_PI));
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
            source<<"mixed3 d2r = "<<sign<<"*"<<cc.doubleToString(sin(M_PI/3.0))<<"*(c1i-c2i);\n";
            source<<"mixed3 d2i = "<<sign<<"*"<<cc.doubleToString(sin(M_PI/3.0))<<"*(c2r-c1r);\n";
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
        source<<"SYNC_WARPS;\n";
    }
    source<<"}\n";
    return source.str();
}

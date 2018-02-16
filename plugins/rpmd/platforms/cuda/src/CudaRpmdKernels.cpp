/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2018 Stanford University and the Authors.      *
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

#include "CudaRpmdKernels.h"
#include "CudaRpmdKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "CudaIntegrationUtilities.h"
#include "CudaExpressionUtilities.h"
#include "CudaKernelSources.h"
#include "CudaNonbondedUtilities.h"
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

void CudaIntegrateRPMDStepKernel::initialize(const System& system, const RPMDIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    numCopies = integrator.getNumCopies();
    numParticles = system.getNumParticles();
    workgroupSize = numCopies;
    if (numCopies != findFFTDimension(numCopies))
        throw OpenMMException("RPMDIntegrator: the number of copies must be a multiple of powers of 2, 3, and 5.");
    int paddedParticles = cu.getPaddedNumAtoms();
    bool useDoublePrecision = (cu.getUseDoublePrecision() || cu.getUseMixedPrecision());
    int elementSize = (useDoublePrecision ? sizeof(double4) : sizeof(float4));
    forces.initialize<long long>(cu, numCopies*paddedParticles*3, "rpmdForces");
    positions.initialize(cu, numCopies*paddedParticles, elementSize, "rpmdPositions");
    velocities.initialize(cu, numCopies*paddedParticles, elementSize, "rpmdVelocities");
    cu.getIntegrationUtilities().initRandomNumberGenerator((unsigned int) integrator.getRandomNumberSeed());
    
    // Fill in the posq and velm arrays with safe values to avoid a risk of nans.
    
    if (useDoublePrecision) {
        vector<double4> temp(positions.getSize());
        for (int i = 0; i < positions.getSize(); i++)
            temp[i] = make_double4(0, 0, 0, 0);
        positions.upload(temp);
        for (int i = 0; i < velocities.getSize(); i++)
            temp[i] = make_double4(0, 0, 0, 1);
        velocities.upload(temp);
    }
    else {
        vector<float4> temp(positions.getSize());
        for (int i = 0; i < positions.getSize(); i++)
            temp[i] = make_float4(0, 0, 0, 0);
        positions.upload(temp);
        for (int i = 0; i < velocities.getSize(); i++)
            temp[i] = make_float4(0, 0, 0, 1);
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
    if (maxContractedCopies > 0) {
        contractedForces.initialize<long long>(cu, maxContractedCopies*paddedParticles*3, "rpmdContractedForces");
        contractedPositions.initialize(cu, maxContractedCopies*paddedParticles, elementSize, "rpmdContractedPositions");
    }

    // Create kernels.
    
    map<string, string> defines;
    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    defines["NUM_COPIES"] = cu.intToString(numCopies);
    defines["THREAD_BLOCK_SIZE"] = cu.intToString(workgroupSize);
    defines["HBAR"] = cu.doubleToString(1.054571628e-34*AVOGADRO/(1000*1e-12));
    defines["SCALE"] = cu.doubleToString(1.0/sqrt((double) numCopies));
    defines["M_PI"] = cu.doubleToString(M_PI);
    map<string, string> replacements;
    replacements["FFT_Q_FORWARD"] = createFFT(numCopies, "q", true);
    replacements["FFT_Q_BACKWARD"] = createFFT(numCopies, "q", false);
    replacements["FFT_V_FORWARD"] = createFFT(numCopies, "v", true);
    replacements["FFT_V_BACKWARD"] = createFFT(numCopies, "v", false);
    CUmodule module = cu.createModule(cu.replaceStrings(CudaKernelSources::vectorOps+CudaRpmdKernelSources::rpmd, replacements), defines, "");
    pileKernel = cu.getKernel(module, "applyPileThermostat");
    stepKernel = cu.getKernel(module, "integrateStep");
    velocitiesKernel = cu.getKernel(module, "advanceVelocities");
    copyToContextKernel = cu.getKernel(module, "copyDataToContext");
    copyFromContextKernel = cu.getKernel(module, "copyDataFromContext");
    translateKernel = cu.getKernel(module, "applyCellTranslations");
    
    // Create kernels for doing contractions.
    
    for (auto& g : groupsByCopies) {
        int copies = g.first;
        replacements.clear();
        replacements["NUM_CONTRACTED_COPIES"] = cu.intToString(copies);
        replacements["POS_SCALE"] = cu.doubleToString(1.0/numCopies);
        replacements["FORCE_SCALE"] = cu.doubleToString(0x100000000/(double) copies);
        replacements["FFT_Q_FORWARD"] = createFFT(numCopies, "q", true);
        replacements["FFT_Q_BACKWARD"] = createFFT(copies, "q", false);
        replacements["FFT_F_FORWARD"] = createFFT(copies, "f", true);
        replacements["FFT_F_BACKWARD"] = createFFT(numCopies, "f", false);
        module = cu.createModule(cu.replaceStrings(CudaKernelSources::vectorOps+CudaRpmdKernelSources::rpmdContraction, replacements), defines, "");
        positionContractionKernels[copies] = cu.getKernel(module, "contractPositions");
        forceContractionKernels[copies] = cu.getKernel(module, "contractForces");
    }
}

void CudaIntegrateRPMDStepKernel::execute(ContextImpl& context, const RPMDIntegrator& integrator, bool forcesAreValid) {
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    
    // Loop over copies and compute the force on each one.
    
    if (!forcesAreValid)
        computeForces(context);
    
    // Apply the PILE-L thermostat.
    
    bool useDoublePrecision = (cu.getUseDoublePrecision() || cu.getUseMixedPrecision());
    double dt = integrator.getStepSize();
    float dtFloat = (float) dt;
    void* dtPtr = (useDoublePrecision ? (void*) &dt : (void*) &dtFloat);
    double kT = integrator.getTemperature()*BOLTZ;
    float kTFloat = (float) kT;
    void* kTPtr = (useDoublePrecision ? (void*) &kT : (void*) &kTFloat);
    double friction = integrator.getFriction();
    float frictionFloat = (float) friction;
    void* frictionPtr = (useDoublePrecision ? (void*) &friction : (void*) &frictionFloat);
    int randomIndex = integration.prepareRandomNumbers(numParticles*numCopies);
    void* pileArgs[] = {&velocities.getDevicePointer(), &integration.getRandom().getDevicePointer(), &randomIndex, dtPtr, kTPtr, frictionPtr};
    if (integrator.getApplyThermostat())
        cu.executeKernel(pileKernel, pileArgs, numParticles*numCopies, workgroupSize);

    // Update positions and velocities.
    
    void* stepArgs[] = {&positions.getDevicePointer(), &velocities.getDevicePointer(), &forces.getDevicePointer(), dtPtr, kTPtr};
    cu.executeKernel(stepKernel, stepArgs, numParticles*numCopies, workgroupSize);

    // Calculate forces based on the updated positions.
    
    computeForces(context);
    
    // Update velocities.

    void* velocitiesArgs[] = {&velocities.getDevicePointer(), &forces.getDevicePointer(), dtPtr};
    cu.executeKernel(velocitiesKernel, velocitiesArgs, numParticles*numCopies, workgroupSize);

    // Apply the PILE-L thermostat again.

    if (integrator.getApplyThermostat()) {
        randomIndex = integration.prepareRandomNumbers(numParticles*numCopies);
        cu.executeKernel(pileKernel, pileArgs, numParticles*numCopies, workgroupSize);
    }

    // Update the time and step count.

    cu.setTime(cu.getTime()+dt);
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
    if (cu.getAtomsWereReordered() && cu.getNonbondedUtilities().getUsePeriodic()) {
        // Atoms may have been translated into a different periodic box, so apply
        // the same translation to all the beads.

        int i = numCopies-1;
        void* args[] = {&positions.getDevicePointer(), &cu.getPosq().getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer(), &i};
        cu.executeKernel(translateKernel, args, cu.getNumAtoms());
    }
}

void CudaIntegrateRPMDStepKernel::computeForces(ContextImpl& context) {
    // Compute forces from all groups that didn't have a specified contraction.

    for (int i = 0; i < numCopies; i++) {
        void* copyToContextArgs[] = {&velocities.getDevicePointer(), &cu.getVelm().getDevicePointer(), &positions.getDevicePointer(),
                &cu.getPosq().getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer(), &i};
        cu.executeKernel(copyToContextKernel, copyToContextArgs, cu.getNumAtoms());
        context.computeVirtualSites();
        Vec3 initialBox[3];
        context.getPeriodicBoxVectors(initialBox[0], initialBox[1], initialBox[2]);
        context.updateContextState();
        Vec3 finalBox[3];
        context.getPeriodicBoxVectors(finalBox[0], finalBox[1], finalBox[2]);
        if (initialBox[0] != finalBox[0] || initialBox[1] != finalBox[1] || initialBox[2] != finalBox[2])
            throw OpenMMException("Standard barostats cannot be used with RPMDIntegrator.  Use RPMDMonteCarloBarostat instead.");
        context.calcForcesAndEnergy(true, false, groupsNotContracted);
        void* copyFromContextArgs[] = {&cu.getForce().getDevicePointer(), &forces.getDevicePointer(), &cu.getVelm().getDevicePointer(),
                &velocities.getDevicePointer(), &cu.getPosq().getDevicePointer(), &positions.getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer(), &i};
        cu.executeKernel(copyFromContextKernel, copyFromContextArgs, cu.getNumAtoms());
    }
    
    // Now loop over contractions and compute forces from them.
    
    for (auto& g : groupsByCopies) {
        int copies = g.first;
        int groupFlags = g.second;
        
        // Find the contracted positions.
        
        void* contractPosArgs[] = {&positions.getDevicePointer(), &contractedPositions.getDevicePointer()};
        cu.executeKernel(positionContractionKernels[copies], contractPosArgs, numParticles*numCopies, workgroupSize);

        // Compute forces.

        for (int i = 0; i < copies; i++) {
            void* copyToContextArgs[] = {&velocities.getDevicePointer(), &cu.getVelm().getDevicePointer(), &contractedPositions.getDevicePointer(),
                    &cu.getPosq().getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer(), &i};
            cu.executeKernel(copyToContextKernel, copyToContextArgs, cu.getNumAtoms());
            context.computeVirtualSites();
            context.calcForcesAndEnergy(true, false, groupFlags);
            void* copyFromContextArgs[] = {&cu.getForce().getDevicePointer(), &contractedForces.getDevicePointer(), &cu.getVelm().getDevicePointer(),
                   &velocities.getDevicePointer(), &cu.getPosq().getDevicePointer(), &contractedPositions.getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer(), &i};
            cu.executeKernel(copyFromContextKernel, copyFromContextArgs, cu.getNumAtoms());
        }
        
        // Apply the forces to the original copies.
        
        void* contractForceArgs[] = {&forces.getDevicePointer(), &contractedForces.getDevicePointer()};
        cu.executeKernel(forceContractionKernels[copies], contractForceArgs, numParticles*numCopies, workgroupSize);
    }
    if (groupsByCopies.size() > 0) {
        // Ensure the Context contains the positions from the last copy, since we'll assume that later.
        
        int i = numCopies-1;
        void* copyToContextArgs[] = {&velocities.getDevicePointer(), &cu.getVelm().getDevicePointer(), &positions.getDevicePointer(),
                &cu.getPosq().getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer(), &i};
        cu.executeKernel(copyToContextKernel, copyToContextArgs, cu.getNumAtoms());
    }
}

double CudaIntegrateRPMDStepKernel::computeKineticEnergy(ContextImpl& context, const RPMDIntegrator& integrator) {
    return cu.getIntegrationUtilities().computeKineticEnergy(0);
}

void CudaIntegrateRPMDStepKernel::setPositions(int copy, const vector<Vec3>& pos) {
    if (!positions.isInitialized())
        throw OpenMMException("RPMDIntegrator: Cannot set positions before the integrator is added to a Context");
    if (pos.size() != numParticles)
        throw OpenMMException("RPMDIntegrator: wrong number of values passed to setPositions()");

    // Adjust the positions based on the current cell offsets.
    
    const vector<int>& order = cu.getAtomIndex();
    double4 periodicBoxSize = cu.getPeriodicBoxSize();
    vector<Vec3> offsetPos(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        int4 offset = cu.getPosCellOffsets()[i];
        offsetPos[order[i]] = pos[order[i]] + Vec3(offset.x*periodicBoxSize.x, offset.y*periodicBoxSize.y, offset.z*periodicBoxSize.z);
    }

    // Record the positions.

    CUresult result;
    if (cu.getUseDoublePrecision()) {
        vector<double4> posq(cu.getPaddedNumAtoms());
        cu.getPosq().download(posq);
        for (int i = 0; i < numParticles; i++)
            posq[i] = make_double4(offsetPos[i][0], offsetPos[i][1], offsetPos[i][2], posq[i].w);
        result = cuMemcpyHtoD(positions.getDevicePointer()+copy*cu.getPaddedNumAtoms()*sizeof(double4), &posq[0], numParticles*sizeof(double4));
    }
    else if (cu.getUseMixedPrecision()) {
        vector<float4> posqf(cu.getPaddedNumAtoms());
        cu.getPosq().download(posqf);
        vector<double4> posq(cu.getPaddedNumAtoms());
        for (int i = 0; i < numParticles; i++)
            posq[i] = make_double4(offsetPos[i][0], offsetPos[i][1], offsetPos[i][2], posqf[i].w);
        result = cuMemcpyHtoD(positions.getDevicePointer()+copy*cu.getPaddedNumAtoms()*sizeof(double4), &posq[0], numParticles*sizeof(double4));
    }
    else {
        vector<float4> posq(cu.getPaddedNumAtoms());
        cu.getPosq().download(posq);
        for (int i = 0; i < numParticles; i++)
            posq[i] = make_float4((float) offsetPos[i][0], (float) offsetPos[i][1], (float) offsetPos[i][2], posq[i].w);
        result = cuMemcpyHtoD(positions.getDevicePointer()+copy*cu.getPaddedNumAtoms()*sizeof(float4), &posq[0], numParticles*sizeof(float4));
    }
    if (result != CUDA_SUCCESS) {
        std::stringstream str;
        str<<"Error uploading array "<<positions.getName()<<": "<<CudaContext::getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
}

void CudaIntegrateRPMDStepKernel::setVelocities(int copy, const vector<Vec3>& vel) {
    if (!velocities.isInitialized())
        throw OpenMMException("RPMDIntegrator: Cannot set velocities before the integrator is added to a Context");
    if (vel.size() != numParticles)
        throw OpenMMException("RPMDIntegrator: wrong number of values passed to setVelocities()");
    CUresult result;
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
        vector<double4> velm(cu.getPaddedNumAtoms());
        cu.getVelm().download(velm);
        for (int i = 0; i < numParticles; i++)
            velm[i] = make_double4(vel[i][0], vel[i][1], vel[i][2], velm[i].w);
        result = cuMemcpyHtoD(velocities.getDevicePointer()+copy*cu.getPaddedNumAtoms()*sizeof(double4), &velm[0], numParticles*sizeof(double4));
    }
    else {
        vector<float4> velm(cu.getPaddedNumAtoms());
        cu.getVelm().download(velm);
        for (int i = 0; i < numParticles; i++)
            velm[i] = make_float4((float) vel[i][0], (float) vel[i][1], (float) vel[i][2], velm[i].w);
        result = cuMemcpyHtoD(velocities.getDevicePointer()+copy*cu.getPaddedNumAtoms()*sizeof(float4), &velm[0], numParticles*sizeof(float4));
    }
    if (result != CUDA_SUCCESS) {
        std::stringstream str;
        str<<"Error uploading array "<<velocities.getName()<<": "<<CudaContext::getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
}

void CudaIntegrateRPMDStepKernel::copyToContext(int copy, ContextImpl& context) {
    void* copyArgs[] = {&velocities.getDevicePointer(), &cu.getVelm().getDevicePointer(), &positions.getDevicePointer(),
            &cu.getPosq().getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer(), &copy};
    cu.executeKernel(copyToContextKernel, copyArgs, cu.getNumAtoms());
}

string CudaIntegrateRPMDStepKernel::createFFT(int size, const string& variable, bool forward) {
    stringstream source;
    int stage = 0;
    int L = size;
    int m = 1;
    string sign = (forward ? "1.0f" : "-1.0f");
    string multReal = (forward ? "multiplyComplexRealPart" : "multiplyComplexRealPartConj");
    string multImag = (forward ? "multiplyComplexImagPart" : "multiplyComplexImagPartConj");

    source<<"{\n";
    source<<"mixed3* real0 = "<<variable<<"real;\n";
    source<<"mixed3* imag0 = "<<variable<<"imag;\n";
    source<<"mixed3* real1 = &temp[blockStart];\n";
    source<<"mixed3* imag1 = &temp[blockStart+blockDim.x];\n";

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
            throw OpenMMException("Illegal size for FFT: "+cu.intToString(size));
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
            source<<"mixed3 d2r = "<<cu.doubleToString(sin(0.4*M_PI))<<"*(c1r-c4r);\n";
            source<<"mixed3 d2i = "<<cu.doubleToString(sin(0.4*M_PI))<<"*(c1i-c4i);\n";
            source<<"mixed3 d3r = "<<cu.doubleToString(sin(0.4*M_PI))<<"*(c2r-c3r);\n";
            source<<"mixed3 d3i = "<<cu.doubleToString(sin(0.4*M_PI))<<"*(c2i-c3i);\n";
            source<<"mixed3 d4r = d0r+d1r;\n";
            source<<"mixed3 d4i = d0i+d1i;\n";
            source<<"mixed3 d5r = "<<cu.doubleToString(0.25*sqrt(5.0))<<"*(d0r-d1r);\n";
            source<<"mixed3 d5i = "<<cu.doubleToString(0.25*sqrt(5.0))<<"*(d0i-d1i);\n";
            source<<"mixed3 d6r = c0r-0.25f*d4r;\n";
            source<<"mixed3 d6i = c0i-0.25f*d4i;\n";
            source<<"mixed3 d7r = d6r+d5r;\n";
            source<<"mixed3 d7i = d6i+d5i;\n";
            source<<"mixed3 d8r = d6r-d5r;\n";
            source<<"mixed3 d8i = d6i-d5i;\n";
            string coeff = cu.doubleToString(sin(0.2*M_PI)/sin(0.4*M_PI));
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
            source<<"mixed3 d2r = "<<sign<<"*"<<cu.doubleToString(sin(M_PI/3.0))<<"*(c1i-c2i);\n";
            source<<"mixed3 d2i = "<<sign<<"*"<<cu.doubleToString(sin(M_PI/3.0))<<"*(c2r-c1r);\n";
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
        source<<"__syncthreads();\n";
        source<<"}\n";
        ++stage;
    }

    // Create the kernel.

    if (stage%2 == 1) {
        source<<"real0[indexInBlock] = real1[indexInBlock];\n";
        source<<"imag0[indexInBlock] = imag1[indexInBlock];\n";
    }
    source<<"}\n";
    return source.str();
}

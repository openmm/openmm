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

#include "OpenCLRpmdKernels.h"
#include "OpenCLRpmdKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "OpenCLIntegrationUtilities.h"
#include "OpenCLExpressionUtilities.h"
#include "OpenCLFFT3D.h"
#include "OpenCLNonbondedUtilities.h"
#include "SimTKOpenMMRealType.h"

using namespace OpenMM;
using namespace std;

void OpenCLIntegrateRPMDStepKernel::initialize(const System& system, const RPMDIntegrator& integrator) {
    cl.getPlatformData().initializeContexts(system);
    numCopies = integrator.getNumCopies();
    numParticles = system.getNumParticles();
    workgroupSize = numCopies;
    if (numCopies != OpenCLFFT3D::findLegalDimension(numCopies))
        throw OpenMMException("RPMDIntegrator: the number of copies must be a multiple of powers of 2, 3, and 5.");
    int paddedParticles = cl.getPaddedNumAtoms();
    int forceElementSize = (cl.getUseDoublePrecision() ? sizeof(mm_double4) : sizeof(mm_float4));
    forces.initialize(cl, numCopies*paddedParticles, forceElementSize, "rpmdForces");
    bool useDoublePrecision = (cl.getUseDoublePrecision() || cl.getUseMixedPrecision());
    int elementSize = (useDoublePrecision ? sizeof(mm_double4) : sizeof(mm_float4));
    positions.initialize(cl, numCopies*paddedParticles, elementSize, "rpmdPositions");
    velocities.initialize(cl, numCopies*paddedParticles, elementSize, "rpmdVelocities");
    cl.getIntegrationUtilities().initRandomNumberGenerator((unsigned int) integrator.getRandomNumberSeed());
    
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
        if (copies != OpenCLFFT3D::findLegalDimension(copies))
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
        contractedForces.initialize(cl, maxContractedCopies*paddedParticles, forceElementSize, "rpmdContractedForces");
        contractedPositions.initialize(cl, maxContractedCopies*paddedParticles, elementSize, "rpmdContractedPositions");
    }

    // Create kernels.
    
    map<string, string> defines;
    defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
    defines["NUM_COPIES"] = cl.intToString(numCopies);
    defines["THREAD_BLOCK_SIZE"] = cl.intToString(workgroupSize);
    defines["HBAR"] = cl.doubleToString(1.054571628e-34*AVOGADRO/(1000*1e-12));
    defines["SCALE"] = cl.doubleToString(1.0/sqrt((double) numCopies));
    defines["M_PI"] = cl.doubleToString(M_PI);
    map<string, string> replacements;
    replacements["FFT_Q_FORWARD"] = createFFT(numCopies, "q", true);
    replacements["FFT_Q_BACKWARD"] = createFFT(numCopies, "q", false);
    replacements["FFT_V_FORWARD"] = createFFT(numCopies, "v", true);
    replacements["FFT_V_BACKWARD"] = createFFT(numCopies, "v", false);
    cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLRpmdKernelSources::rpmd, replacements), defines, "");
    pileKernel = cl::Kernel(program, "applyPileThermostat");
    stepKernel = cl::Kernel(program, "integrateStep");
    velocitiesKernel = cl::Kernel(program, "advanceVelocities");
    copyToContextKernel = cl::Kernel(program, "copyDataToContext");
    copyFromContextKernel = cl::Kernel(program, "copyDataFromContext");
    translateKernel = cl::Kernel(program, "applyCellTranslations");
    
    // Create kernels for doing contractions.
    
    for (auto& g : groupsByCopies) {
        int copies = g.first;
        replacements.clear();
        replacements["NUM_CONTRACTED_COPIES"] = cl.intToString(copies);
        replacements["POS_SCALE"] = cl.doubleToString(1.0/numCopies);
        replacements["FORCE_SCALE"] = cl.doubleToString(1.0/copies);
        replacements["FFT_Q_FORWARD"] = createFFT(numCopies, "q", true);
        replacements["FFT_Q_BACKWARD"] = createFFT(copies, "q", false);
        replacements["FFT_F_FORWARD"] = createFFT(copies, "f", true);
        replacements["FFT_F_BACKWARD"] = createFFT(numCopies, "f", false);
        program = cl.createProgram(cl.replaceStrings(OpenCLRpmdKernelSources::rpmdContraction, replacements), defines, "");
        positionContractionKernels[copies] = cl::Kernel(program, "contractPositions");
        forceContractionKernels[copies] = cl::Kernel(program, "contractForces");
    }
}

void OpenCLIntegrateRPMDStepKernel::initializeKernels(ContextImpl& context) {
    hasInitializedKernel = true;
    pileKernel.setArg<cl::Buffer>(0, velocities.getDeviceBuffer());
    stepKernel.setArg<cl::Buffer>(0, positions.getDeviceBuffer());
    stepKernel.setArg<cl::Buffer>(1, velocities.getDeviceBuffer());
    stepKernel.setArg<cl::Buffer>(2, forces.getDeviceBuffer());
    velocitiesKernel.setArg<cl::Buffer>(0, velocities.getDeviceBuffer());
    velocitiesKernel.setArg<cl::Buffer>(1, forces.getDeviceBuffer());
    translateKernel.setArg<cl::Buffer>(0, positions.getDeviceBuffer());
    translateKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
    translateKernel.setArg<cl::Buffer>(2, cl.getAtomIndexArray().getDeviceBuffer());
    copyToContextKernel.setArg<cl::Buffer>(0, velocities.getDeviceBuffer());
    copyToContextKernel.setArg<cl::Buffer>(1, cl.getVelm().getDeviceBuffer());
    copyToContextKernel.setArg<cl::Buffer>(3, cl.getPosq().getDeviceBuffer());
    copyToContextKernel.setArg<cl::Buffer>(4, cl.getAtomIndexArray().getDeviceBuffer());
    copyFromContextKernel.setArg<cl::Buffer>(0, cl.getForce().getDeviceBuffer());
    copyFromContextKernel.setArg<cl::Buffer>(2, cl.getVelm().getDeviceBuffer());
    copyFromContextKernel.setArg<cl::Buffer>(3, velocities.getDeviceBuffer());
    copyFromContextKernel.setArg<cl::Buffer>(4, cl.getPosq().getDeviceBuffer());
    copyFromContextKernel.setArg<cl::Buffer>(6, cl.getAtomIndexArray().getDeviceBuffer());
    for (auto& g : groupsByCopies) {
        int copies = g.first;
        positionContractionKernels[copies].setArg<cl::Buffer>(0, positions.getDeviceBuffer());
        positionContractionKernels[copies].setArg<cl::Buffer>(1, contractedPositions.getDeviceBuffer());
        forceContractionKernels[copies].setArg<cl::Buffer>(0, forces.getDeviceBuffer());
        forceContractionKernels[copies].setArg<cl::Buffer>(1, contractedForces.getDeviceBuffer());
    }
}

void OpenCLIntegrateRPMDStepKernel::execute(ContextImpl& context, const RPMDIntegrator& integrator, bool forcesAreValid) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    if (!hasInitializedKernel)
        initializeKernels(context);
    
    // Loop over copies and compute the force on each one.
    
    if (!forcesAreValid)
        computeForces(context);
    
    // Apply the PILE-L thermostat.
    
    bool useDoublePrecision = (cl.getUseDoublePrecision() || cl.getUseMixedPrecision());
    const double dt = integrator.getStepSize();
    pileKernel.setArg<cl_uint>(2, integration.prepareRandomNumbers(numParticles*numCopies));
    pileKernel.setArg<cl::Buffer>(1, integration.getRandom().getDeviceBuffer()); // Do this *after* prepareRandomNumbers(), which might rebuild the array.
    if (useDoublePrecision) {
        pileKernel.setArg<cl_double>(3, dt);
        pileKernel.setArg<cl_double>(4, integrator.getTemperature()*BOLTZ);
        pileKernel.setArg<cl_double>(5, integrator.getFriction());
        stepKernel.setArg<cl_double>(3, dt);
        stepKernel.setArg<cl_double>(4, integrator.getTemperature()*BOLTZ);
        velocitiesKernel.setArg<cl_double>(2, dt);
    }
    else {
        pileKernel.setArg<cl_float>(3, (cl_float) dt);
        pileKernel.setArg<cl_float>(4, (cl_float) (integrator.getTemperature()*BOLTZ));
        pileKernel.setArg<cl_float>(5, (cl_float) integrator.getFriction());
        stepKernel.setArg<cl_float>(3, (cl_float) dt);
        stepKernel.setArg<cl_float>(4, (cl_float) (integrator.getTemperature()*BOLTZ));
        velocitiesKernel.setArg<cl_float>(2, (cl_float) dt);
    }
    if (integrator.getApplyThermostat())
        cl.executeKernel(pileKernel, numParticles*numCopies, workgroupSize);

    // Update positions and velocities.
    
    cl.executeKernel(stepKernel, numParticles*numCopies, workgroupSize);

    // Calculate forces based on the updated positions.
    
    computeForces(context);
    
    // Update velocities.
    cl.executeKernel(velocitiesKernel, numParticles*numCopies, workgroupSize);

    // Apply the PILE-L thermostat again.

    if (integrator.getApplyThermostat()) {
        pileKernel.setArg<cl_uint>(2, integration.prepareRandomNumbers(numParticles*numCopies));
        cl.executeKernel(pileKernel, numParticles*numCopies, workgroupSize);
    }

    // Update the time and step count.

    cl.setTime(cl.getTime()+dt);
    cl.setStepCount(cl.getStepCount()+1);
    cl.reorderAtoms();
    if (cl.getAtomsWereReordered() && cl.getNonbondedUtilities().getUsePeriodic()) {
        // Atoms may have been translated into a different periodic box, so apply
        // the same translation to all the beads.

        translateKernel.setArg<cl_int>(3, numCopies-1);
        cl.executeKernel(translateKernel, cl.getNumAtoms());
    }
}

void OpenCLIntegrateRPMDStepKernel::computeForces(ContextImpl& context) {
    // Compute forces from all groups that didn't have a specified contraction.

    copyToContextKernel.setArg<cl::Buffer>(2, positions.getDeviceBuffer());
    copyFromContextKernel.setArg<cl::Buffer>(1, forces.getDeviceBuffer());
    copyFromContextKernel.setArg<cl::Buffer>(5, positions.getDeviceBuffer());
    for (int i = 0; i < numCopies; i++) {
        copyToContextKernel.setArg<cl_int>(5, i);
        cl.executeKernel(copyToContextKernel, cl.getNumAtoms());
        context.computeVirtualSites();
        Vec3 initialBox[3];
        context.getPeriodicBoxVectors(initialBox[0], initialBox[1], initialBox[2]);
        context.updateContextState();
        Vec3 finalBox[3];
        context.getPeriodicBoxVectors(finalBox[0], finalBox[1], finalBox[2]);
        if (initialBox[0] != finalBox[0] || initialBox[1] != finalBox[1] || initialBox[2] != finalBox[2])
            throw OpenMMException("Standard barostats cannot be used with RPMDIntegrator.  Use RPMDMonteCarloBarostat instead.");
        context.calcForcesAndEnergy(true, false, groupsNotContracted);
        copyFromContextKernel.setArg<cl_int>(7, i);
        cl.executeKernel(copyFromContextKernel, cl.getNumAtoms());
    }
    
    // Now loop over contractions and compute forces from them.
    
    if (groupsByCopies.size() > 0) {
        copyToContextKernel.setArg<cl::Buffer>(2, contractedPositions.getDeviceBuffer());
        copyFromContextKernel.setArg<cl::Buffer>(1, contractedForces.getDeviceBuffer());
        copyFromContextKernel.setArg<cl::Buffer>(5, contractedPositions.getDeviceBuffer());
        for (auto& g : groupsByCopies) {
            int copies = g.first;
            int groupFlags = g.second;

            // Find the contracted positions.

            cl.executeKernel(positionContractionKernels[copies], numParticles*numCopies, workgroupSize);

            // Compute forces.

            for (int i = 0; i < copies; i++) {
                copyToContextKernel.setArg<cl_int>(5, i);
                cl.executeKernel(copyToContextKernel, cl.getNumAtoms());
                context.computeVirtualSites();
                context.calcForcesAndEnergy(true, false, groupFlags);
                copyFromContextKernel.setArg<cl_int>(7, i);
                cl.executeKernel(copyFromContextKernel, cl.getNumAtoms());
            }

            // Apply the forces to the original copies.

            cl.executeKernel(forceContractionKernels[copies], numParticles*numCopies, workgroupSize);
        }
    }
    if (groupsByCopies.size() > 0) {
        // Ensure the Context contains the positions from the last copy, since we'll assume that later.
        
        copyToContextKernel.setArg<cl::Buffer>(2, positions.getDeviceBuffer());
        copyToContextKernel.setArg<cl_int>(5, numCopies-1);
        cl.executeKernel(copyToContextKernel, cl.getNumAtoms());
    }
}

double OpenCLIntegrateRPMDStepKernel::computeKineticEnergy(ContextImpl& context, const RPMDIntegrator& integrator) {
    return cl.getIntegrationUtilities().computeKineticEnergy(0);
}

void OpenCLIntegrateRPMDStepKernel::setPositions(int copy, const vector<Vec3>& pos) {
    if (!positions.isInitialized())
        throw OpenMMException("RPMDIntegrator: Cannot set positions before the integrator is added to a Context");
    if (pos.size() != numParticles)
        throw OpenMMException("RPMDIntegrator: wrong number of values passed to setPositions()");

    // Adjust the positions based on the current cell offsets.
    
    const vector<int>& order = cl.getAtomIndex();
    mm_double4 periodicBoxSize = cl.getPeriodicBoxSizeDouble();
    vector<Vec3> offsetPos(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        mm_int4 offset = cl.getPosCellOffsets()[i];
        offsetPos[order[i]] = pos[order[i]] + Vec3(offset.x*periodicBoxSize.x, offset.y*periodicBoxSize.y, offset.z*periodicBoxSize.z);
    }

    // Record the positions.

    if (cl.getUseDoublePrecision()) {
        vector<mm_double4> posq(cl.getPaddedNumAtoms());
        cl.getPosq().download(posq);
        for (int i = 0; i < numParticles; i++)
            posq[i] = mm_double4(offsetPos[i][0], offsetPos[i][1], offsetPos[i][2], posq[i].w);
        cl.getQueue().enqueueWriteBuffer(positions.getDeviceBuffer(), CL_TRUE, copy*cl.getPaddedNumAtoms()*sizeof(mm_double4), numParticles*sizeof(mm_double4), &posq[0]);
    }
    else if (cl.getUseMixedPrecision()) {
        vector<mm_float4> posqf(cl.getPaddedNumAtoms());
        cl.getPosq().download(posqf);
        vector<mm_double4> posq(cl.getPaddedNumAtoms());
        for (int i = 0; i < numParticles; i++)
            posq[i] = mm_double4(offsetPos[i][0], offsetPos[i][1], offsetPos[i][2], posqf[i].w);
        cl.getQueue().enqueueWriteBuffer(positions.getDeviceBuffer(), CL_TRUE, copy*cl.getPaddedNumAtoms()*sizeof(mm_double4), numParticles*sizeof(mm_double4), &posq[0]);
    }
    else {
        vector<mm_float4> posq(cl.getPaddedNumAtoms());
        cl.getPosq().download(posq);
        for (int i = 0; i < numParticles; i++)
            posq[i] = mm_float4((cl_float) offsetPos[i][0], (cl_float) offsetPos[i][1], (cl_float) offsetPos[i][2], posq[i].w);
        cl.getQueue().enqueueWriteBuffer(positions.getDeviceBuffer(), CL_TRUE, copy*cl.getPaddedNumAtoms()*sizeof(mm_float4), numParticles*sizeof(mm_float4), &posq[0]);
    }
}

void OpenCLIntegrateRPMDStepKernel::setVelocities(int copy, const vector<Vec3>& vel) {
    if (!velocities.isInitialized())
        throw OpenMMException("RPMDIntegrator: Cannot set velocities before the integrator is added to a Context");
    if (vel.size() != numParticles)
        throw OpenMMException("RPMDIntegrator: wrong number of values passed to setVelocities()");
    if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
        vector<mm_double4> velm(cl.getPaddedNumAtoms());
        cl.getVelm().download(velm);
        for (int i = 0; i < numParticles; i++)
            velm[i] = mm_double4(vel[i][0], vel[i][1], vel[i][2], velm[i].w);
        cl.getQueue().enqueueWriteBuffer(velocities.getDeviceBuffer(), CL_TRUE, copy*cl.getPaddedNumAtoms()*sizeof(mm_double4), numParticles*sizeof(mm_double4), &velm[0]);
    }
    else {
        vector<mm_float4> velm(cl.getPaddedNumAtoms());
        cl.getVelm().download(velm);
        for (int i = 0; i < numParticles; i++)
            velm[i] = mm_float4((cl_float) vel[i][0], (cl_float) vel[i][1], (cl_float) vel[i][2], velm[i].w);
        cl.getQueue().enqueueWriteBuffer(velocities.getDeviceBuffer(), CL_TRUE, copy*cl.getPaddedNumAtoms()*sizeof(mm_float4), numParticles*sizeof(mm_float4), &velm[0]);
    }
}

void OpenCLIntegrateRPMDStepKernel::copyToContext(int copy, ContextImpl& context) {
    if (!hasInitializedKernel)
        initializeKernels(context);
    copyToContextKernel.setArg<cl::Buffer>(2, positions.getDeviceBuffer());
    copyToContextKernel.setArg<cl_int>(5, copy);
    cl.executeKernel(copyToContextKernel, cl.getNumAtoms());
}

string OpenCLIntegrateRPMDStepKernel::createFFT(int size, const string& variable, bool forward) {
    stringstream source;
    int stage = 0;
    int L = size;
    int m = 1;
    string sign = (forward ? "1.0f" : "-1.0f");
    string multReal = (forward ? "multiplyComplexRealPart" : "multiplyComplexRealPartConj");
    string multImag = (forward ? "multiplyComplexImagPart" : "multiplyComplexImagPartConj");

    source<<"{\n";
    source<<"__local mixed4* real0 = "<<variable<<"real;\n";
    source<<"__local mixed4* imag0 = "<<variable<<"imag;\n";
    source<<"__local mixed4* real1 = &temp[blockStart];\n";
    source<<"__local mixed4* imag1 = &temp[blockStart+get_local_size(0)];\n";

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
            throw OpenMMException("Illegal size for FFT: "+cl.intToString(size));
        source<<"{\n";
        L = L/radix;
        source<<"// Pass "<<(stage+1)<<" (radix "<<radix<<")\n";
        source<<"if (indexInBlock < "<<(L*m)<<") {\n";
        source<<"int i = indexInBlock;\n";
        source<<"int j = i/"<<m<<";\n";
        if (radix == 5) {
            source<<"mixed4 c0r = real"<<input<<"[i];\n";
            source<<"mixed4 c0i = imag"<<input<<"[i];\n";
            source<<"mixed4 c1r = real"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"mixed4 c1i = imag"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"mixed4 c2r = real"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"mixed4 c2i = imag"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"mixed4 c3r = real"<<input<<"[i+"<<(3*L*m)<<"];\n";
            source<<"mixed4 c3i = imag"<<input<<"[i+"<<(3*L*m)<<"];\n";
            source<<"mixed4 c4r = real"<<input<<"[i+"<<(4*L*m)<<"];\n";
            source<<"mixed4 c4i = imag"<<input<<"[i+"<<(4*L*m)<<"];\n";
            source<<"mixed4 d0r = c1r+c4r;\n";
            source<<"mixed4 d0i = c1i+c4i;\n";
            source<<"mixed4 d1r = c2r+c3r;\n";
            source<<"mixed4 d1i = c2i+c3i;\n";
            source<<"mixed4 d2r = "<<cl.doubleToString(sin(0.4*M_PI))<<"*(c1r-c4r);\n";
            source<<"mixed4 d2i = "<<cl.doubleToString(sin(0.4*M_PI))<<"*(c1i-c4i);\n";
            source<<"mixed4 d3r = "<<cl.doubleToString(sin(0.4*M_PI))<<"*(c2r-c3r);\n";
            source<<"mixed4 d3i = "<<cl.doubleToString(sin(0.4*M_PI))<<"*(c2i-c3i);\n";
            source<<"mixed4 d4r = d0r+d1r;\n";
            source<<"mixed4 d4i = d0i+d1i;\n";
            source<<"mixed4 d5r = "<<cl.doubleToString(0.25*sqrt(5.0))<<"*(d0r-d1r);\n";
            source<<"mixed4 d5i = "<<cl.doubleToString(0.25*sqrt(5.0))<<"*(d0i-d1i);\n";
            source<<"mixed4 d6r = c0r-0.25f*d4r;\n";
            source<<"mixed4 d6i = c0i-0.25f*d4i;\n";
            source<<"mixed4 d7r = d6r+d5r;\n";
            source<<"mixed4 d7i = d6i+d5i;\n";
            source<<"mixed4 d8r = d6r-d5r;\n";
            source<<"mixed4 d8i = d6i-d5i;\n";
            string coeff = cl.doubleToString(sin(0.2*M_PI)/sin(0.4*M_PI));
            source<<"mixed4 d9r = "<<sign<<"*(d2i+"<<coeff<<"*d3i);\n";
            source<<"mixed4 d9i = "<<sign<<"*(-d2r-"<<coeff<<"*d3r);\n";
            source<<"mixed4 d10r = "<<sign<<"*("<<coeff<<"*d2i-d3i);\n";
            source<<"mixed4 d10i = "<<sign<<"*(d3r-"<<coeff<<"*d2r);\n";
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
            source<<"mixed4 c0r = real"<<input<<"[i];\n";
            source<<"mixed4 c0i = imag"<<input<<"[i];\n";
            source<<"mixed4 c1r = real"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"mixed4 c1i = imag"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"mixed4 c2r = real"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"mixed4 c2i = imag"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"mixed4 c3r = real"<<input<<"[i+"<<(3*L*m)<<"];\n";
            source<<"mixed4 c3i = imag"<<input<<"[i+"<<(3*L*m)<<"];\n";
            source<<"mixed4 d0r = c0r+c2r;\n";
            source<<"mixed4 d0i = c0i+c2i;\n";
            source<<"mixed4 d1r = c0r-c2r;\n";
            source<<"mixed4 d1i = c0i-c2i;\n";
            source<<"mixed4 d2r = c1r+c3r;\n";
            source<<"mixed4 d2i = c1i+c3i;\n";
            source<<"mixed4 d3r = "<<sign<<"*(c1i-c3i);\n";
            source<<"mixed4 d3i = "<<sign<<"*(c3r-c1r);\n";
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
            source<<"mixed4 c0r = real"<<input<<"[i];\n";
            source<<"mixed4 c0i = imag"<<input<<"[i];\n";
            source<<"mixed4 c1r = real"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"mixed4 c1i = imag"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"mixed4 c2r = real"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"mixed4 c2i = imag"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"mixed4 d0r = c1r+c2r;\n";
            source<<"mixed4 d0i = c1i+c2i;\n";
            source<<"mixed4 d1r = c0r-0.5f*d0r;\n";
            source<<"mixed4 d1i = c0i-0.5f*d0i;\n";
            source<<"mixed4 d2r = "<<sign<<"*"<<cl.doubleToString(sin(M_PI/3.0))<<"*(c1i-c2i);\n";
            source<<"mixed4 d2i = "<<sign<<"*"<<cl.doubleToString(sin(M_PI/3.0))<<"*(c2r-c1r);\n";
            source<<"real"<<output<<"[i+2*j*"<<m<<"] = c0r+d0r;\n";
            source<<"imag"<<output<<"[i+2*j*"<<m<<"] = c0i+d0i;\n";
            source<<"real"<<output<<"[i+(2*j+1)*"<<m<<"] = "<<multReal<<"(w[j*"<<size<<"/"<<(3*L)<<"], d1r+d2r, d1i+d2i);\n";
            source<<"imag"<<output<<"[i+(2*j+1)*"<<m<<"] = "<<multImag<<"(w[j*"<<size<<"/"<<(3*L)<<"], d1r+d2r, d1i+d2i);\n";
            source<<"real"<<output<<"[i+(2*j+2)*"<<m<<"] = "<<multReal<<"(w[j*"<<(2*size)<<"/"<<(3*L)<<"], d1r-d2r, d1i-d2i);\n";
            source<<"imag"<<output<<"[i+(2*j+2)*"<<m<<"] = "<<multImag<<"(w[j*"<<(2*size)<<"/"<<(3*L)<<"], d1r-d2r, d1i-d2i);\n";
        }
        else if (radix == 2) {
            source<<"mixed4 c0r = real"<<input<<"[i];\n";
            source<<"mixed4 c0i = imag"<<input<<"[i];\n";
            source<<"mixed4 c1r = real"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"mixed4 c1i = imag"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"real"<<output<<"[i+j*"<<m<<"] = c0r+c1r;\n";
            source<<"imag"<<output<<"[i+j*"<<m<<"] = c0i+c1i;\n";
            source<<"real"<<output<<"[i+(j+1)*"<<m<<"] = "<<multReal<<"(w[j*"<<size<<"/"<<(2*L)<<"], c0r-c1r, c0i-c1i);\n";
            source<<"imag"<<output<<"[i+(j+1)*"<<m<<"] = "<<multImag<<"(w[j*"<<size<<"/"<<(2*L)<<"], c0r-c1r, c0i-c1i);\n";
        }
        source<<"}\n";
        m = m*radix;
        source<<"barrier(CLK_LOCAL_MEM_FENCE);\n";
        source<<"}\n";
        ++stage;
    }

    // Create the kernel.

    if (stage%2 == 1) {
        source<<"real0[indexInBlock] = real1[indexInBlock];\n";
        source<<"imag0[indexInBlock] = imag1[indexInBlock];\n";
        source<<"barrier(CLK_LOCAL_MEM_FENCE);\n";
    }
    source<<"}\n";
    return source.str();
}

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
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
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"

using namespace OpenMM;
using namespace std;

OpenCLIntegrateRPMDStepKernel::~OpenCLIntegrateRPMDStepKernel() {
    if (forces != NULL)
        delete forces;
    if (positions != NULL)
        delete positions;
    if (velocities != NULL)
        delete velocities;
}
void OpenCLIntegrateRPMDStepKernel::initialize(const System& system, const RPMDIntegrator& integrator) {
    cl.getPlatformData().initializeContexts(system);
    numCopies = integrator.getNumCopies();
    numParticles = system.getNumParticles();
    workgroupSize = numCopies;
    while (workgroupSize <= 128-numCopies)
        workgroupSize += numCopies;
    if (numCopies != OpenCLFFT3D::findLegalDimension(numCopies))
        throw OpenMMException("RPMDIntegrator: the number of copies must be a multiple of powers of 2, 3, and 5.");
    int paddedParticles = cl.getPaddedNumAtoms();
    forces = new OpenCLArray<mm_float4>(cl, numCopies*paddedParticles, "rpmdForces");
    positions = new OpenCLArray<mm_float4>(cl, numCopies*paddedParticles, "rpmdPositions");
    velocities = new OpenCLArray<mm_float4>(cl, numCopies*paddedParticles, "rpmdVelocities");
    cl.getIntegrationUtilities().initRandomNumberGenerator((unsigned int) integrator.getRandomNumberSeed());
    
    // Fill in the posq and velm arrays with safe values to avoid a risk of nans.
    
    vector<mm_float4> temp(positions->getSize());
    for (int i = 0; i < positions->getSize(); i++)
        temp[i] = mm_float4(0, 0, 0, 0);
    positions->upload(temp);
    for (int i = 0; i < velocities->getSize(); i++)
        temp[i] = mm_float4(0, 0, 0, 1);
    velocities->upload(temp);

    // Create kernels.
    
    map<string, string> defines;
    defines["NUM_ATOMS"] = OpenCLExpressionUtilities::intToString(cl.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = OpenCLExpressionUtilities::intToString(cl.getPaddedNumAtoms());
    defines["NUM_COPIES"] = OpenCLExpressionUtilities::intToString(numCopies);
    defines["HBAR"] = OpenCLExpressionUtilities::doubleToString(1.054571628e-34*AVOGADRO/(1000*1e-12));
    defines["SCALE"] = OpenCLExpressionUtilities::doubleToString(1.0/sqrt((double) numCopies));
    defines["M_PI"] = OpenCLExpressionUtilities::doubleToString(M_PI);
    map<string, string> replacements;
    replacements["FFT_Q_FORWARD"] = createFFT(numCopies, "q", true);
    replacements["FFT_Q_BACKWARD"] = createFFT(numCopies, "q", false);
    replacements["FFT_V_FORWARD"] = createFFT(numCopies, "v", true);
    replacements["FFT_V_BACKWARD"] = createFFT(numCopies, "v", false);
    cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLRpmdKernelSources::rpmd, replacements), defines, "");
    pileKernel = cl::Kernel(program, "applyPileThermostat");
    stepKernel = cl::Kernel(program, "integrateStep");
    velocitiesKernel = cl::Kernel(program, "advanceVelocities");
    copyToContextKernel = cl::Kernel(program, "copyToContext");
    copyFromContextKernel = cl::Kernel(program, "copyFromContext");
}

void OpenCLIntegrateRPMDStepKernel::execute(ContextImpl& context, const RPMDIntegrator& integrator, bool forcesAreValid) {
    const System& system = context.getSystem();
    const int paddedParticles = cl.getPaddedNumAtoms();
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();

    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        pileKernel.setArg<cl::Buffer>(0, velocities->getDeviceBuffer());
        pileKernel.setArg(1, 2*workgroupSize*sizeof(mm_float4), NULL);
        pileKernel.setArg(2, 2*workgroupSize*sizeof(mm_float4), NULL);
        pileKernel.setArg(3, numCopies*sizeof(mm_float2), NULL);
        pileKernel.setArg<cl::Buffer>(4, integration.getRandom().getDeviceBuffer());
        stepKernel.setArg<cl::Buffer>(0, positions->getDeviceBuffer());
        stepKernel.setArg<cl::Buffer>(1, velocities->getDeviceBuffer());
        stepKernel.setArg<cl::Buffer>(2, forces->getDeviceBuffer());
        stepKernel.setArg(3, 2*workgroupSize*sizeof(mm_float4), NULL);
        stepKernel.setArg(4, 2*workgroupSize*sizeof(mm_float4), NULL);
        stepKernel.setArg(5, 2*workgroupSize*sizeof(mm_float4), NULL);
        stepKernel.setArg(6, numCopies*sizeof(mm_float2), NULL);
        velocitiesKernel.setArg<cl::Buffer>(0, velocities->getDeviceBuffer());
        velocitiesKernel.setArg<cl::Buffer>(1, forces->getDeviceBuffer());
    }
    
    // Loop over copies and compute the force on each one.
    
    copyToContextKernel.setArg<cl::Buffer>(0, positions->getDeviceBuffer());
    copyToContextKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
    copyToContextKernel.setArg<cl::Buffer>(2, cl.getAtomIndex().getDeviceBuffer());
    copyFromContextKernel.setArg<cl::Buffer>(0, cl.getForce().getDeviceBuffer());
    copyFromContextKernel.setArg<cl::Buffer>(1, forces->getDeviceBuffer());
    copyFromContextKernel.setArg<cl::Buffer>(2, cl.getAtomIndex().getDeviceBuffer());
    if (!forcesAreValid) {
        for (int i = 0; i < numCopies; i++) {
            copyToContextKernel.setArg<cl_int>(3, i);
            cl.executeKernel(copyToContextKernel, cl.getNumAtoms());
            context.calcForcesAndEnergy(true, false);
            copyFromContextKernel.setArg<cl_int>(3, i);
            cl.executeKernel(copyFromContextKernel, cl.getNumAtoms());
        }
    }
    
    // Apply the PILE-L thermostat.
    
    const double dt = integrator.getStepSize();
    pileKernel.setArg<cl_uint>(5, integration.prepareRandomNumbers(numParticles*numCopies));
    pileKernel.setArg<cl_float>(6, dt);
    pileKernel.setArg<cl_float>(7, integrator.getTemperature()*BOLTZ);
    pileKernel.setArg<cl_float>(8, integrator.getFriction());
    cl.executeKernel(pileKernel, numParticles*numCopies, workgroupSize);

    // Update positions and velocities.
    
    stepKernel.setArg<cl_float>(7, dt);
    stepKernel.setArg<cl_float>(8, integrator.getTemperature()*BOLTZ);
    cl.executeKernel(stepKernel, numParticles*numCopies, workgroupSize);

    // Calculate forces based on the updated positions.
    
    for (int i = 0; i < numCopies; i++) {
        copyToContextKernel.setArg<cl_int>(3, i);
        cl.executeKernel(copyToContextKernel, cl.getNumAtoms());
        context.calcForcesAndEnergy(true, false);
        copyFromContextKernel.setArg<cl_int>(3, i);
        cl.executeKernel(copyFromContextKernel, cl.getNumAtoms());
    }
    
    // Update velocities.
    velocitiesKernel.setArg<cl_float>(2, dt);
    cl.executeKernel(velocitiesKernel, numParticles*numCopies, workgroupSize);

    // Apply the PILE-L thermostat again.

    pileKernel.setArg<cl_uint>(5, integration.prepareRandomNumbers(numParticles*numCopies));
    cl.executeKernel(pileKernel, numParticles*numCopies, workgroupSize);

    // Update the time and step count.

    cl.setTime(cl.getTime()+dt);
    cl.setStepCount(cl.getStepCount()+1);
}

void OpenCLIntegrateRPMDStepKernel::setPositions(int copy, const vector<Vec3>& pos) {
    if (positions == NULL)
        throw OpenMMException("RPMDIntegrator: Cannot set positions before the integrator is added to a Context");
    if (pos.size() != numParticles)
        throw OpenMMException("RPMDIntegrator: wrong number of values passed to setPositions()");
    vector<mm_float4> posq(numParticles);
    for (int i = 0; i < numParticles; i++)
        posq[i] = mm_float4(pos[i][0], pos[i][1], pos[i][2], cl.getPosq()[i].w);
    cl.getQueue().enqueueWriteBuffer(positions->getDeviceBuffer(), CL_TRUE, copy*cl.getPaddedNumAtoms()*sizeof(mm_float4), numParticles*sizeof(mm_float4), &posq[0]);
}

void OpenCLIntegrateRPMDStepKernel::setVelocities(int copy, const vector<Vec3>& vel) {
    if (velocities == NULL)
        throw OpenMMException("RPMDIntegrator: Cannot set velocities before the integrator is added to a Context");
    if (vel.size() != numParticles)
        throw OpenMMException("RPMDIntegrator: wrong number of values passed to setVelocities()");
    vector<mm_float4> velm(numParticles);
    for (int i = 0; i < numParticles; i++)
        velm[i] = mm_float4(vel[i][0], vel[i][1], vel[i][2], cl.getVelm()[i].w);
    cl.getQueue().enqueueWriteBuffer(velocities->getDeviceBuffer(), CL_TRUE, copy*cl.getPaddedNumAtoms()*sizeof(mm_float4), numParticles*sizeof(mm_float4), &velm[0]);
}

void OpenCLIntegrateRPMDStepKernel::copyToContext(int copy, ContextImpl& context) {
    copyToContextKernel.setArg<cl::Buffer>(0, positions->getDeviceBuffer());
    copyToContextKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
    copyToContextKernel.setArg<cl::Buffer>(2, cl.getAtomIndex().getDeviceBuffer());
    copyToContextKernel.setArg<cl_int>(3, copy);
    cl.executeKernel(copyToContextKernel, cl.getNumAtoms());
    copyToContextKernel.setArg<cl::Buffer>(0, velocities->getDeviceBuffer());
    copyToContextKernel.setArg<cl::Buffer>(1, cl.getVelm().getDeviceBuffer());
    cl.executeKernel(copyToContextKernel, cl.getNumAtoms());
}

string OpenCLIntegrateRPMDStepKernel::createFFT(int size, const string& variable, bool forward) {
    stringstream source;
    int unfactored = size;
    int stage = 0;
    int L = size;
    int m = 1;
    string sign = (forward ? "1.0f" : "-1.0f");
    string multReal = (forward ? "multiplyComplexRealPart" : "multiplyComplexRealPartConj");
    string multImag = (forward ? "multiplyComplexImagPart" : "multiplyComplexImagPartConj");

    source<<"{\n";
    source<<"__local float4* real0 = "<<variable<<"real;\n";
    source<<"__local float4* imag0 = "<<variable<<"imag;\n";
    source<<"__local float4* real1 = &temp[blockStart];\n";
    source<<"__local float4* imag1 = &temp[blockStart+get_local_size(0)];\n";

    // Factor size, generating an appropriate block of code for each factor.

    while (unfactored > 1) {
        int input = stage%2;
        int output = 1-input;
        source<<"{\n";
        if (unfactored%5 == 0) {
            L = L/5;
            source<<"// Pass "<<(stage+1)<<" (radix 5)\n";
            source<<"if (indexInBlock < "<<(L*m)<<") {\n";
            source<<"int i = indexInBlock;\n";
            source<<"int j = i/"<<m<<";\n";
            source<<"float4 c0r = real"<<input<<"[i];\n";
            source<<"float4 c0i = imag"<<input<<"[i];\n";
            source<<"float4 c1r = real"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"float4 c1i = imag"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"float4 c2r = real"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"float4 c2i = imag"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"float4 c3r = real"<<input<<"[i+"<<(3*L*m)<<"];\n";
            source<<"float4 c3i = imag"<<input<<"[i+"<<(3*L*m)<<"];\n";
            source<<"float4 c4r = real"<<input<<"[i+"<<(4*L*m)<<"];\n";
            source<<"float4 c4i = imag"<<input<<"[i+"<<(4*L*m)<<"];\n";
            source<<"float4 d0r = c1r+c4r;\n";
            source<<"float4 d0i = c1i+c4i;\n";
            source<<"float4 d1r = c2r+c3r;\n";
            source<<"float4 d1i = c2i+c3i;\n";
            source<<"float4 d2r = "<<OpenCLExpressionUtilities::doubleToString(sin(0.4*M_PI))<<"*(c1r-c4r);\n";
            source<<"float4 d2i = "<<OpenCLExpressionUtilities::doubleToString(sin(0.4*M_PI))<<"*(c1i-c4i);\n";
            source<<"float4 d3r = "<<OpenCLExpressionUtilities::doubleToString(sin(0.4*M_PI))<<"*(c2r-c3r);\n";
            source<<"float4 d3i = "<<OpenCLExpressionUtilities::doubleToString(sin(0.4*M_PI))<<"*(c2i-c3i);\n";
            source<<"float4 d4r = d0r+d1r;\n";
            source<<"float4 d4i = d0i+d1i;\n";
            source<<"float4 d5r = "<<OpenCLExpressionUtilities::doubleToString(0.25*sqrt(5.0))<<"*(d0r-d1r);\n";
            source<<"float4 d5i = "<<OpenCLExpressionUtilities::doubleToString(0.25*sqrt(5.0))<<"*(d0i-d1i);\n";
            source<<"float4 d6r = c0r-0.25f*d4r;\n";
            source<<"float4 d6i = c0i-0.25f*d4i;\n";
            source<<"float4 d7r = d6r+d5r;\n";
            source<<"float4 d7i = d6i+d5i;\n";
            source<<"float4 d8r = d6r-d5r;\n";
            source<<"float4 d8i = d6i-d5i;\n";
            string coeff = OpenCLExpressionUtilities::doubleToString(sin(0.2*M_PI)/sin(0.4*M_PI));
            source<<"float4 d9r = "<<sign<<"*(d2i+"<<coeff<<"*d3i);\n";
            source<<"float4 d9i = "<<sign<<"*(-d2r-"<<coeff<<"*d3r);\n";
            source<<"float4 d10r = "<<sign<<"*("<<coeff<<"*d2i-d3i);\n";
            source<<"float4 d10i = "<<sign<<"*(d3r-"<<coeff<<"*d2r);\n";
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
            source<<"}\n";
            m = m*5;
            unfactored /= 5;
        }
        else if (unfactored%4 == 0) {
            L = L/4;
            source<<"// Pass "<<(stage+1)<<" (radix 4)\n";
            source<<"if (indexInBlock < "<<(L*m)<<") {\n";
            source<<"int i = indexInBlock;\n";
            source<<"int j = i/"<<m<<";\n";
            source<<"float4 c0r = real"<<input<<"[i];\n";
            source<<"float4 c0i = imag"<<input<<"[i];\n";
            source<<"float4 c1r = real"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"float4 c1i = imag"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"float4 c2r = real"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"float4 c2i = imag"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"float4 c3r = real"<<input<<"[i+"<<(3*L*m)<<"];\n";
            source<<"float4 c3i = imag"<<input<<"[i+"<<(3*L*m)<<"];\n";
            source<<"float4 d0r = c0r+c2r;\n";
            source<<"float4 d0i = c0i+c2i;\n";
            source<<"float4 d1r = c0r-c2r;\n";
            source<<"float4 d1i = c0i-c2i;\n";
            source<<"float4 d2r = c1r+c3r;\n";
            source<<"float4 d2i = c1i+c3i;\n";
            source<<"float4 d3r = "<<sign<<"*(c1i-c3i);\n";
            source<<"float4 d3i = "<<sign<<"*(c3r-c1r);\n";
            source<<"real"<<output<<"[i+3*j*"<<m<<"] = d0r+d2r;\n";
            source<<"imag"<<output<<"[i+3*j*"<<m<<"] = d0i+d2i;\n";
            source<<"real"<<output<<"[i+(3*j+1)*"<<m<<"] = "<<multReal<<"(w[j*"<<size<<"/"<<(4*L)<<"], d1r+d3r, d1i+d3i);\n";
            source<<"imag"<<output<<"[i+(3*j+1)*"<<m<<"] = "<<multImag<<"(w[j*"<<size<<"/"<<(4*L)<<"], d1r+d3r, d1i+d3i);\n";
            source<<"real"<<output<<"[i+(3*j+2)*"<<m<<"] = "<<multReal<<"(w[j*"<<(2*size)<<"/"<<(4*L)<<"], d0r-d2r, d0i-d2i);\n";
            source<<"imag"<<output<<"[i+(3*j+2)*"<<m<<"] = "<<multImag<<"(w[j*"<<(2*size)<<"/"<<(4*L)<<"], d0r-d2r, d0i-d2i);\n";
            source<<"real"<<output<<"[i+(3*j+3)*"<<m<<"] = "<<multReal<<"(w[j*"<<(3*size)<<"/"<<(4*L)<<"], d1r-d3r, d1i-d3i);\n";
            source<<"imag"<<output<<"[i+(3*j+3)*"<<m<<"] = "<<multImag<<"(w[j*"<<(3*size)<<"/"<<(4*L)<<"], d1r-d3r, d1i-d3i);\n";
            source<<"}\n";
            m = m*4;
            unfactored /= 4;
        }
        else if (unfactored%3 == 0) {
            L = L/3;
            source<<"// Pass "<<(stage+1)<<" (radix 3)\n";
            source<<"if (indexInBlock < "<<(L*m)<<") {\n";
            source<<"int i = indexInBlock;\n";
            source<<"int j = i/"<<m<<";\n";
            source<<"float4 c0r = real"<<input<<"[i];\n";
            source<<"float4 c0i = imag"<<input<<"[i];\n";
            source<<"float4 c1r = real"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"float4 c1i = imag"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"float4 c2r = real"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"float4 c2i = imag"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"float4 d0r = c1r+c2r;\n";
            source<<"float4 d0i = c1i+c2i;\n";
            source<<"float4 d1r = c0r-0.5f*d0r;\n";
            source<<"float4 d1i = c0i-0.5f*d0i;\n";
            source<<"float4 d2r = "<<sign<<"*"<<OpenCLExpressionUtilities::doubleToString(sin(M_PI/3.0))<<"*(c1i-c2i);\n";
            source<<"float4 d2i = "<<sign<<"*"<<OpenCLExpressionUtilities::doubleToString(sin(M_PI/3.0))<<"*(c2r-c1r);\n";
            source<<"real"<<output<<"[i+2*j*"<<m<<"] = c0r+d0r;\n";
            source<<"imag"<<output<<"[i+2*j*"<<m<<"] = c0i+d0i;\n";
            source<<"real"<<output<<"[i+(2*j+1)*"<<m<<"] = "<<multReal<<"(w[j*"<<size<<"/"<<(3*L)<<"], d1r+d2r, d1i+d2i);\n";
            source<<"imag"<<output<<"[i+(2*j+1)*"<<m<<"] = "<<multImag<<"(w[j*"<<size<<"/"<<(3*L)<<"], d1r+d2r, d1i+d2i);\n";
            source<<"real"<<output<<"[i+(2*j+2)*"<<m<<"] = "<<multReal<<"(w[j*"<<(2*size)<<"/"<<(3*L)<<"], d1r-d2r, d1i-d2i);\n";
            source<<"imag"<<output<<"[i+(2*j+2)*"<<m<<"] = "<<multImag<<"(w[j*"<<(2*size)<<"/"<<(3*L)<<"], d1r-d2r, d1i-d2i);\n";
            source<<"}\n";
            m = m*3;
            unfactored /= 3;
        }
        else if (unfactored%2 == 0) {
            L = L/2;
            source<<"// Pass "<<(stage+1)<<" (radix 2)\n";
            source<<"if (indexInBlock < "<<(L*m)<<") {\n";
            source<<"int i = indexInBlock;\n";
            source<<"int j = i/"<<m<<";\n";
            source<<"float4 c0r = real"<<input<<"[i];\n";
            source<<"float4 c0i = imag"<<input<<"[i];\n";
            source<<"float4 c1r = real"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"float4 c1i = imag"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"real"<<output<<"[i+j*"<<m<<"] = c0r+c1r;\n";
            source<<"imag"<<output<<"[i+j*"<<m<<"] = c0i+c1i;\n";
            source<<"real"<<output<<"[i+(j+1)*"<<m<<"] = "<<multReal<<"(w[j*"<<size<<"/"<<(2*L)<<"], c0r-c1r, c0i-c1i);\n";
            source<<"imag"<<output<<"[i+(j+1)*"<<m<<"] = "<<multImag<<"(w[j*"<<size<<"/"<<(2*L)<<"], c0r-c1r, c0i-c1i);\n";
            source<<"}\n";
            m = m*2;
            unfactored /= 2;
        }
        else
            throw OpenMMException("Illegal size for FFT: "+OpenCLExpressionUtilities::intToString(size));
        source<<"barrier(CLK_LOCAL_MEM_FENCE);\n";
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

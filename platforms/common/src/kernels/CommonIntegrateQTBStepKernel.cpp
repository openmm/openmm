/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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

#include "openmm/common/CommonIntegrateQTBStepKernel.h"
#include "openmm/common/CommonKernelUtilities.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/QTBIntegratorUtilities.h"
#include "CommonKernelSources.h"
#include "SimTKOpenMMUtilities.h"

using namespace OpenMM;
using namespace std;

void CommonIntegrateQTBStepKernel::initialize(const System& system, const QTBIntegrator& integrator) {
    cc.initializeContexts();
    ContextSelector selector(cc);
    cc.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    int numParticles = system.getNumParticles();
    dt = integrator.getStepSize();
    friction = integrator.getFriction();
    segmentLength = (int) ceil(integrator.getSegmentLength()/integrator.getStepSize());
    numFreq = (3*segmentLength+1)/2;
    noise.initialize<float>(cc, 3*3*segmentLength*numParticles, "noise");
    SimTKOpenMMUtilities::setRandomNumberSeed(integrator.getRandomNumberSeed());
    vector<float> noiseVec(noise.getSize());
    for (int i = 0; i < noiseVec.size(); i++)
        noiseVec[i] = (float) SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
    noise.upload(noiseVec);
    int elementSize = (cc.getUseMixedPrecision() || cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    randomForce.initialize(cc, 3*segmentLength*numParticles, elementSize, "randomForce");
    segmentVelocity.initialize(cc, 3*segmentLength*numParticles, elementSize, "segmentVelocity");
    oldDelta.initialize(cc, cc.getPaddedNumAtoms(), 4*elementSize, "oldDelta");
    thetad.initialize(cc, numFreq, elementSize, "thetad");
    cutoffFunction.initialize(cc, numFreq, elementSize, "cutoffFunction");
    vector<double> cf(numFreq);
    double cutoff = integrator.getCutoffFrequency();
    double cutoffWidth = cutoff/100;
    for (int i = 0; i < numFreq; i++) {
        double w = M_PI*i/(numFreq*dt);
        cf[i] = 1.0/(1.0+exp((w-cutoff)/cutoffWidth));
    }
    cutoffFunction.upload(cf, true);

    map<string, string> defines, replacements;
    defines["M_PI"] = cc.doubleToString(M_PI);
    defines["FFT_LENGTH"] = cc.intToString(3*segmentLength);
    int outputIndex;
    replacements["FFT_FORWARD"] = createFFT(3*segmentLength, 0, outputIndex, true);
    replacements["RECIP_DATA"] = (outputIndex == 0 ? "data0" : "data1");
    replacements["FFT_BACKWARD"] = createFFT(3*segmentLength, outputIndex, outputIndex, false);
    ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::qtb, replacements), defines);
    kernel1 = program->createKernel("integrateQTBPart1");
    kernel2 = program->createKernel("integrateQTBPart2");
    kernel3 = program->createKernel("integrateQTBPart3");
    noiseKernel = program->createKernel("generateNoise");
    forceKernel = program->createKernel("generateRandomForce");
    stepIndex = 0;
    prevTemp = -1.0;
}

void CommonIntegrateQTBStepKernel::execute(ContextImpl& context, const QTBIntegrator& integrator) {
    ContextSelector selector(cc);
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    int numAtoms = cc.getNumAtoms();
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    double temperature = integrator.getTemperature();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        bool useDouble = cc.getUseDoublePrecision() || cc.getUseMixedPrecision();
        kernel1->addArg(numAtoms);
        kernel1->addArg(paddedNumAtoms);
        if (useDouble)
            kernel1->addArg(dt);
        else
            kernel1->addArg((float) dt);
        kernel1->addArg();
        kernel1->addArg(cc.getVelm());
        kernel1->addArg(cc.getLongForceBuffer());
        kernel1->addArg(segmentVelocity);
        kernel2->addArg(numAtoms);
        if (useDouble) {
            kernel2->addArg(dt);
            kernel2->addArg(friction);
        }
        else {
            kernel2->addArg((float) dt);
            kernel2->addArg((float) friction);
        }
        kernel2->addArg();
        kernel2->addArg(cc.getVelm());
        kernel2->addArg(integration.getPosDelta());
        kernel2->addArg(oldDelta);
        kernel2->addArg(randomForce);
        kernel3->addArg(numAtoms);
        if (useDouble)
            kernel3->addArg(dt);
        else
            kernel3->addArg((float) dt);
        kernel3->addArg(cc.getPosq());
        kernel3->addArg(cc.getVelm());
        kernel3->addArg(integration.getPosDelta());
        kernel3->addArg(oldDelta);
        if (cc.getUseMixedPrecision())
            kernel3->addArg(cc.getPosqCorrection());
        noiseKernel->addArg(numAtoms);
        noiseKernel->addArg(segmentLength);
        noiseKernel->addArg(noise);
        noiseKernel->addArg(integration.getRandom());
        noiseKernel->addArg(); // Random index will be set just before it is executed.
        forceKernel->addArg(numAtoms);
        forceKernel->addArg(segmentLength);
        if (useDouble) {
            forceKernel->addArg(dt);
            forceKernel->addArg(friction);
        }
        else {
            forceKernel->addArg((float) dt);
            forceKernel->addArg((float) friction);
        }
        forceKernel->addArg(noise);
        forceKernel->addArg(randomForce);
        forceKernel->addArg(cc.getVelm());
        forceKernel->addArg(thetad);
        forceKernel->addArg(cutoffFunction);
        forceKernel->addArg();
    }
    cc.getIntegrationUtilities().setNextStepSize(dt);

    // Update the random force at the start of each segment.

    if (stepIndex%segmentLength == 0) {
        if (temperature != prevTemp) {
            vector<double> thetaVec, thetadVec;
            QTBIntegratorUtilities::calculateSpectrum(temperature, friction, dt, numFreq, thetaVec, thetadVec, cc.getThreadPool());
            thetad.upload(thetadVec, true);
            prevTemp = temperature;
        }
        noiseKernel->setArg(4, integration.prepareRandomNumbers(numAtoms*segmentLength));
        noiseKernel->execute(numAtoms*segmentLength);
        forceKernel->setArg(9, (float) (BOLTZ*integrator.getTemperature()));
        forceKernel->execute(3*numAtoms*64, 64);
        stepIndex = 0;
    }

    // Perform the integration.

    kernel1->setArg(3, stepIndex);
    kernel2->setArg(3, stepIndex);
    kernel1->execute(numAtoms);
    integration.applyVelocityConstraints(integrator.getConstraintTolerance());
    kernel2->execute(numAtoms);
    integration.applyConstraints(integrator.getConstraintTolerance());
    kernel3->execute(numAtoms);
    integration.computeVirtualSites();

    // Update the time and step count.

    cc.setTime(cc.getTime()+dt);
    cc.setStepCount(cc.getStepCount()+1);
    stepIndex++;
    cc.reorderAtoms();

    // Reduce UI lag.

    flushPeriodically(cc);
}

double CommonIntegrateQTBStepKernel::computeKineticEnergy(ContextImpl& context, const QTBIntegrator& integrator) {
    return cc.getIntegrationUtilities().computeKineticEnergy(0.0);
}

void CommonIntegrateQTBStepKernel::getAdaptedFriction(ContextImpl& context, int particle, std::vector<double>& friction) const {
}

void CommonIntegrateQTBStepKernel::setAdaptedFriction(ContextImpl& context, int particle, const std::vector<double>& friction) {
}

void CommonIntegrateQTBStepKernel::createCheckpoint(ContextImpl& context, ostream& stream) const {
}

void CommonIntegrateQTBStepKernel::loadCheckpoint(ContextImpl& context, istream& stream) {
}

string CommonIntegrateQTBStepKernel::createFFT(int size, int inputIndex, int& outputIndex, bool forward) {
    stringstream source;
    int stage = 0;
    int L = size;
    int m = 1;
    string sign = (forward ? "1" : "-1");

    // Factor size, generating an appropriate block of code for each factor.

    while (L > 1) {
        int input = (inputIndex+stage)%2;
        int output = 1-input;
        int radix;
        if (L%7 == 0)
            radix = 7;
        else if (L%5 == 0)
            radix = 5;
        else if (L%4 == 0)
            radix = 4;
        else if (L%3 == 0)
            radix = 3;
        else if (L%2 == 0)
            radix = 2;
        else
            throw OpenMMException("QTBIntegrator: Number of steps in a segment must be a multiple of powers of 2, 3, 5, and 7");
        source<<"{\n";
        L = L/radix;
        source<<"// Pass "<<(stage+1)<<" (radix "<<radix<<")\n";
        source<<"for (int i = LOCAL_ID; i < "<<(L*m)<<"; i += LOCAL_SIZE) {\n";
        source<<"int base = i;\n";
        source<<"int j = i/"<<m<<";\n";
        if (radix == 7) {
            source<<"real2 c0 = data"<<input<<"[base];\n";
            source<<"real2 c1 = data"<<input<<"[base+"<<(L*m)<<"];\n";
            source<<"real2 c2 = data"<<input<<"[base+"<<(2*L*m)<<"];\n";
            source<<"real2 c3 = data"<<input<<"[base+"<<(3*L*m)<<"];\n";
            source<<"real2 c4 = data"<<input<<"[base+"<<(4*L*m)<<"];\n";
            source<<"real2 c5 = data"<<input<<"[base+"<<(5*L*m)<<"];\n";
            source<<"real2 c6 = data"<<input<<"[base+"<<(6*L*m)<<"];\n";
            source<<"real2 d0 = c1+c6;\n";
            source<<"real2 d1 = c1-c6;\n";
            source<<"real2 d2 = c2+c5;\n";
            source<<"real2 d3 = c2-c5;\n";
            source<<"real2 d4 = c4+c3;\n";
            source<<"real2 d5 = c4-c3;\n";
            source<<"real2 d6 = d2+d0;\n";
            source<<"real2 d7 = d5+d3;\n";
            source<<"real2 b0 = c0+d6+d4;\n";
            source<<"real2 b1 = "<<cc.doubleToString((cos(2*M_PI/7)+cos(4*M_PI/7)+cos(6*M_PI/7))/3-1)<<"*(d6+d4);\n";
            source<<"real2 b2 = "<<cc.doubleToString((2*cos(2*M_PI/7)-cos(4*M_PI/7)-cos(6*M_PI/7))/3)<<"*(d0-d4);\n";
            source<<"real2 b3 = "<<cc.doubleToString((cos(2*M_PI/7)-2*cos(4*M_PI/7)+cos(6*M_PI/7))/3)<<"*(d4-d2);\n";
            source<<"real2 b4 = "<<cc.doubleToString((cos(2*M_PI/7)+cos(4*M_PI/7)-2*cos(6*M_PI/7))/3)<<"*(d2-d0);\n";
            source<<"real2 b5 = -("<<sign<<")*"<<cc.doubleToString((sin(2*M_PI/7)+sin(4*M_PI/7)-sin(6*M_PI/7))/3)<<"*(d7+d1);\n";
            source<<"real2 b6 = -("<<sign<<")*"<<cc.doubleToString((2*sin(2*M_PI/7)-sin(4*M_PI/7)+sin(6*M_PI/7))/3)<<"*(d1-d5);\n";
            source<<"real2 b7 = -("<<sign<<")*"<<cc.doubleToString((sin(2*M_PI/7)-2*sin(4*M_PI/7)-sin(6*M_PI/7))/3)<<"*(d5-d3);\n";
            source<<"real2 b8 = -("<<sign<<")*"<<cc.doubleToString((sin(2*M_PI/7)+sin(4*M_PI/7)+2*sin(6*M_PI/7))/3)<<"*(d3-d1);\n";
            source<<"real2 t0 = b0+b1;\n";
            source<<"real2 t1 = b2+b3;\n";
            source<<"real2 t2 = b4-b3;\n";
            source<<"real2 t3 = -b2-b4;\n";
            source<<"real2 t4 = b6+b7;\n";
            source<<"real2 t5 = b8-b7;\n";
            source<<"real2 t6 = -b8-b6;\n";
            source<<"real2 t7 = t0+t1;\n";
            source<<"real2 t8 = t0+t2;\n";
            source<<"real2 t9 = t0+t3;\n";
            source<<"real2 t10 = (real2) (t4.y+b5.y, -(t4.x+b5.x));\n";
            source<<"real2 t11 = (real2) (t5.y+b5.y, -(t5.x+b5.x));\n";
            source<<"real2 t12 = (real2) (t6.y+b5.y, -(t6.x+b5.x));\n";
            source<<"data"<<output<<"[base+6*j*"<<m<<"] = b0;\n";
            source<<"data"<<output<<"[base+(6*j+1)*"<<m<<"] = multiplyComplex(w[j*"<<size<<"/"<<(7*L)<<"], t7-t10);\n";
            source<<"data"<<output<<"[base+(6*j+2)*"<<m<<"] = multiplyComplex(w[j*"<<(2*size)<<"/"<<(7*L)<<"], t9-t12);\n";
            source<<"data"<<output<<"[base+(6*j+3)*"<<m<<"] = multiplyComplex(w[j*"<<(3*size)<<"/"<<(7*L)<<"], t8+t11);\n";
            source<<"data"<<output<<"[base+(6*j+4)*"<<m<<"] = multiplyComplex(w[j*"<<(4*size)<<"/"<<(7*L)<<"], t8-t11);\n";
            source<<"data"<<output<<"[base+(6*j+5)*"<<m<<"] = multiplyComplex(w[j*"<<(5*size)<<"/"<<(7*L)<<"], t9+t12);\n";
            source<<"data"<<output<<"[base+(6*j+6)*"<<m<<"] = multiplyComplex(w[j*"<<(6*size)<<"/"<<(7*L)<<"], t7+t10);\n";
        }
        else if (radix == 5) {
            source<<"real2 c0 = data"<<input<<"[base];\n";
            source<<"real2 c1 = data"<<input<<"[base+"<<(L*m)<<"];\n";
            source<<"real2 c2 = data"<<input<<"[base+"<<(2*L*m)<<"];\n";
            source<<"real2 c3 = data"<<input<<"[base+"<<(3*L*m)<<"];\n";
            source<<"real2 c4 = data"<<input<<"[base+"<<(4*L*m)<<"];\n";
            source<<"real2 d0 = c1+c4;\n";
            source<<"real2 d1 = c2+c3;\n";
            source<<"real2 d2 = "<<cc.doubleToString(sin(0.4*M_PI))<<"*(c1-c4);\n";
            source<<"real2 d3 = "<<cc.doubleToString(sin(0.4*M_PI))<<"*(c2-c3);\n";
            source<<"real2 d4 = d0+d1;\n";
            source<<"real2 d5 = "<<cc.doubleToString(0.25*sqrt(5.0))<<"*(d0-d1);\n";
            source<<"real2 d6 = c0-0.25f*d4;\n";
            source<<"real2 d7 = d6+d5;\n";
            source<<"real2 d8 = d6-d5;\n";
            string coeff = cc.doubleToString(sin(0.2*M_PI)/sin(0.4*M_PI));
            source<<"real2 d9 = "<<sign<<"*(real2) (d2.y+"<<coeff<<"*d3.y, -d2.x-"<<coeff<<"*d3.x);\n";
            source<<"real2 d10 = "<<sign<<"*(real2) ("<<coeff<<"*d2.y-d3.y, d3.x-"<<coeff<<"*d2.x);\n";
            source<<"data"<<output<<"[base+4*j*"<<m<<"] = c0+d4;\n";
            source<<"data"<<output<<"[base+(4*j+1)*"<<m<<"] = multiplyComplex(w[j*"<<size<<"/"<<(5*L)<<"], d7+d9);\n";
            source<<"data"<<output<<"[base+(4*j+2)*"<<m<<"] = multiplyComplex(w[j*"<<(2*size)<<"/"<<(5*L)<<"], d8+d10);\n";
            source<<"data"<<output<<"[base+(4*j+3)*"<<m<<"] = multiplyComplex(w[j*"<<(3*size)<<"/"<<(5*L)<<"], d8-d10);\n";
            source<<"data"<<output<<"[base+(4*j+4)*"<<m<<"] = multiplyComplex(w[j*"<<(4*size)<<"/"<<(5*L)<<"], d7-d9);\n";
        }
        else if (radix == 4) {
            source<<"real2 c0 = data"<<input<<"[base];\n";
            source<<"real2 c1 = data"<<input<<"[base+"<<(L*m)<<"];\n";
            source<<"real2 c2 = data"<<input<<"[base+"<<(2*L*m)<<"];\n";
            source<<"real2 c3 = data"<<input<<"[base+"<<(3*L*m)<<"];\n";
            source<<"real2 d0 = c0+c2;\n";
            source<<"real2 d1 = c0-c2;\n";
            source<<"real2 d2 = c1+c3;\n";
            source<<"real2 d3 = "<<sign<<"*(real2) (c1.y-c3.y, c3.x-c1.x);\n";
            source<<"data"<<output<<"[base+3*j*"<<m<<"] = d0+d2;\n";
            source<<"data"<<output<<"[base+(3*j+1)*"<<m<<"] = multiplyComplex(w[j*"<<size<<"/"<<(4*L)<<"], d1+d3);\n";
            source<<"data"<<output<<"[base+(3*j+2)*"<<m<<"] = multiplyComplex(w[j*"<<(2*size)<<"/"<<(4*L)<<"], d0-d2);\n";
            source<<"data"<<output<<"[base+(3*j+3)*"<<m<<"] = multiplyComplex(w[j*"<<(3*size)<<"/"<<(4*L)<<"], d1-d3);\n";
        }
        else if (radix == 3) {
            source<<"real2 c0 = data"<<input<<"[base];\n";
            source<<"real2 c1 = data"<<input<<"[base+"<<(L*m)<<"];\n";
            source<<"real2 c2 = data"<<input<<"[base+"<<(2*L*m)<<"];\n";
            source<<"real2 d0 = c1+c2;\n";
            source<<"real2 d1 = c0-0.5f*d0;\n";
            source<<"real2 d2 = "<<sign<<"*"<<cc.doubleToString(sin(M_PI/3.0))<<"*(real2) (c1.y-c2.y, c2.x-c1.x);\n";
            source<<"data"<<output<<"[base+2*j*"<<m<<"] = c0+d0;\n";
            source<<"data"<<output<<"[base+(2*j+1)*"<<m<<"] = multiplyComplex(w[j*"<<size<<"/"<<(3*L)<<"], d1+d2);\n";
            source<<"data"<<output<<"[base+(2*j+2)*"<<m<<"] = multiplyComplex(w[j*"<<(2*size)<<"/"<<(3*L)<<"], d1-d2);\n";
        }
        else if (radix == 2) {
            source<<"real2 c0 = data"<<input<<"[base];\n";
            source<<"real2 c1 = data"<<input<<"[base+"<<(L*m)<<"];\n";
            source<<"data"<<output<<"[base+j*"<<m<<"] = c0+c1;\n";
            source<<"data"<<output<<"[base+(j+1)*"<<m<<"] = multiplyComplex(w[j*"<<size<<"/"<<(2*L)<<"], c0-c1);\n";
        }
        source<<"}\n";
        m = m*radix;
        source<<"SYNC_THREADS\n";
        source<<"}\n";
        ++stage;
    }
    outputIndex = (inputIndex+stage)%2;
    return source.str();
}

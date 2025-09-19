/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
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
    workspace.initialize(cc, 18*segmentLength*cc.getNumThreadBlocks(), elementSize, "workspace");
    cutoffFunction.initialize(cc, numFreq, elementSize, "cutoffFunction");
    vector<double> cf(numFreq);
    double cutoff = integrator.getCutoffFrequency();
    double cutoffWidth = cutoff/100;
    for (int i = 0; i < numFreq; i++) {
        double w = M_PI*i/(numFreq*dt);
        cf[i] = 1.0/(1.0+exp((w-cutoff)/cutoffWidth));
    }
    cutoffFunction.upload(cf, true);
    vector<double> typeMassVec, typeAdaptationRateVec;
    vector<vector<int> > typeParticles;
    QTBIntegratorUtilities::findTypes(system, integrator, particleTypeVec, typeParticles, typeMassVec, typeAdaptationRateVec);
    particleType.initialize<int>(cc, particleTypeVec.size(), "particleType");
    particleType.upload(particleTypeVec);
    typeAdaptationRate.initialize<float>(cc, typeAdaptationRateVec.size(), "typeAdaptationRate");
    typeAdaptationRate.upload(typeAdaptationRateVec, true);
    vector<int> typeParticleCountVec(typeParticles.size());
    for (int i = 0; i < typeParticleCountVec.size(); i++)
        typeParticleCountVec[i] = typeParticles[i].size();
    typeParticleCount.initialize<int>(cc, typeParticleCountVec.size(), "typeParticleCount");
    typeParticleCount.upload(typeParticleCountVec);
    numTypes = typeParticles.size();
    adaptedFriction.initialize(cc, numFreq*numTypes, elementSize, "adaptedFriction");
    vector<double> adaptedFrictionVec(adaptedFriction.getSize(), friction);
    adaptedFriction.upload(adaptedFrictionVec, true);
    dfdt.initialize(cc, adaptedFriction.getSize(), sizeof(long long), "dfdt");

    // Create the kernels.

    map<string, string> defines, replacements;
    defines["M_PI"] = cc.doubleToString(M_PI);
    int outputIndex;
    replacements["FFT_FORWARD"] = createFFT(3*segmentLength, 0, outputIndex, true);
    replacements["RECIP_DATA"] = (outputIndex == 0 ? "data0" : "data1");
    replacements["FFT_BACKWARD"] = createFFT(3*segmentLength, outputIndex, outputIndex, false);
    replacements["ADAPTATION_FFT"] = createFFT(3*segmentLength, 0, outputIndex, true);
    replacements["ADAPTATION_RECIP"] = (outputIndex == 0 ? "data0" : "data1");
    ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::qtb, replacements), defines);
    kernel1 = program->createKernel("integrateQTBPart1");
    kernel2 = program->createKernel("integrateQTBPart2");
    kernel3 = program->createKernel("integrateQTBPart3");
    noiseKernel = program->createKernel("generateNoise");
    forceKernel = program->createKernel("generateRandomForce");
    adapt1Kernel = program->createKernel("adaptFrictionPart1");
    adapt2Kernel = program->createKernel("adaptFrictionPart2");
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
        kernel1->addArg(cc.getVelm());
        kernel1->addArg(cc.getLongForceBuffer());
        kernel1->addArg(cc.getAtomIndexArray());
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
        kernel2->addArg(segmentVelocity);
        kernel2->addArg(cc.getAtomIndexArray());
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
        forceKernel->addArg(particleType);
        forceKernel->addArg(adaptedFriction);
        forceKernel->addArg(workspace);
        adapt1Kernel->addArg(numAtoms);
        adapt1Kernel->addArg(segmentLength);
        adapt1Kernel->addArg(cc.getVelm());
        adapt1Kernel->addArg(particleType);
        adapt1Kernel->addArg(randomForce);
        adapt1Kernel->addArg(segmentVelocity);
        adapt1Kernel->addArg(adaptedFriction);
        adapt1Kernel->addArg(dfdt);
        adapt1Kernel->addArg(workspace);
        adapt2Kernel->addArg(numTypes);
        adapt2Kernel->addArg(segmentLength);
        if (useDouble)
            adapt2Kernel->addArg(dt);
        else
            adapt2Kernel->addArg((float) dt);
        adapt2Kernel->addArg(typeParticleCount);
        adapt2Kernel->addArg(typeAdaptationRate);
        adapt2Kernel->addArg(adaptedFriction);
        adapt2Kernel->addArg(dfdt);
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
        cc.clearBuffer(dfdt);
        adapt1Kernel->execute(3*numAtoms*128, 128);
        adapt2Kernel->execute(numTypes*128, 128);
        noiseKernel->setArg(4, integration.prepareRandomNumbers(3*numAtoms*segmentLength/4+1));
        noiseKernel->execute(3*numAtoms*128, 128);
        forceKernel->execute(3*numAtoms*128, 128);
        stepIndex = 0;
    }

    // Perform the integration.

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
    ASSERT_VALID_INDEX(particle, particleTypeVec);
    int type = particleTypeVec[particle];
    friction.resize(numFreq);
    ContextSelector selector(cc);
    if (cc.getUseMixedPrecision() || cc.getUseDoublePrecision()) {
        vector<double> adaptedFrictionVec;
        adaptedFriction.download(adaptedFrictionVec);
        for (int i = 0; i < numFreq; i++)
            friction[i] = adaptedFrictionVec[type*numFreq+i];
    }
    else {
        vector<float> adaptedFrictionVec;
        adaptedFriction.download(adaptedFrictionVec);
        for (int i = 0; i < numFreq; i++)
            friction[i] = adaptedFrictionVec[type*numFreq+i];
    }
}

void CommonIntegrateQTBStepKernel::setAdaptedFriction(ContextImpl& context, int particle, const std::vector<double>& friction) {
    ASSERT_VALID_INDEX(particle, particleTypeVec);
    int type = particleTypeVec[particle];
    ContextSelector selector(cc);
    if (cc.getUseMixedPrecision() || cc.getUseDoublePrecision()) {
        vector<double> adaptedFrictionVec;
        adaptedFriction.download(adaptedFrictionVec);
        for (int i = 0; i < numFreq; i++)
            adaptedFrictionVec[type*numFreq+i] = friction[i];
        adaptedFriction.upload(adaptedFrictionVec);
    }
    else {
        vector<float> adaptedFrictionVec;
        adaptedFriction.download(adaptedFrictionVec);
        for (int i = 0; i < numFreq; i++)
            adaptedFrictionVec[type*numFreq+i] = friction[i];
        adaptedFriction.upload(adaptedFrictionVec);
    }
}

void CommonIntegrateQTBStepKernel::createCheckpoint(ContextImpl& context, ostream& stream) const {
    ContextSelector selector(cc);
    stream.write((char*) &stepIndex, sizeof(int));
    vector<float> f;
    noise.download(f);
    stream.write((char*) f.data(), sizeof(float)*f.size());
    if (cc.getUseMixedPrecision() || cc.getUseDoublePrecision()) {
        vector<double> d;
        randomForce.download(d);
        stream.write((char*) d.data(), sizeof(double)*d.size());
        segmentVelocity.download(d);
        stream.write((char*) d.data(), sizeof(double)*d.size());
        adaptedFriction.download(d);
        stream.write((char*) d.data(), sizeof(double)*d.size());
    }
    else {
        randomForce.download(f);
        stream.write((char*) f.data(), sizeof(float)*f.size());
        segmentVelocity.download(f);
        stream.write((char*) f.data(), sizeof(float)*f.size());
        adaptedFriction.download(f);
        stream.write((char*) f.data(), sizeof(float)*f.size());
    }
}

void CommonIntegrateQTBStepKernel::loadCheckpoint(ContextImpl& context, istream& stream) {
    ContextSelector selector(cc);
    stream.read((char*) &stepIndex, sizeof(int));
    vector<float> f(noise.getSize());
    stream.read((char*) f.data(), sizeof(float)*noise.getSize());
    noise.upload(f);
    if (cc.getUseMixedPrecision() || cc.getUseDoublePrecision()) {
        vector<double> d;
        d.resize(randomForce.getSize());
        stream.read((char*) d.data(), sizeof(double)*d.size());
        randomForce.upload(d);
        d.resize(segmentVelocity.getSize());
        stream.read((char*) d.data(), sizeof(double)*d.size());
        segmentVelocity.upload(d);
        d.resize(adaptedFriction.getSize());
        stream.read((char*) d.data(), sizeof(double)*d.size());
        adaptedFriction.upload(d);
    }
    else {
        f.resize(randomForce.getSize());
        stream.read((char*) f.data(), sizeof(float)*f.size());
        randomForce.upload(f);
        f.resize(segmentVelocity.getSize());
        stream.read((char*) f.data(), sizeof(float)*f.size());
        segmentVelocity.upload(f);
        f.resize(adaptedFriction.getSize());
        stream.read((char*) f.data(), sizeof(float)*f.size());
        adaptedFriction.upload(f);
    }
}

string CommonIntegrateQTBStepKernel::createFFT(int size, int inputIndex, int& outputIndex, bool forward) {
    stringstream source;
    int stage = 0;
    int L = size;
    int m = 1;
    string sign = (forward ? "1" : "-1");
    string mult = (forward ? "multiplyComplex" : "multiplyComplexConj");

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
            source<<"mixed2 c0 = data"<<input<<"[base];\n";
            source<<"mixed2 c1 = data"<<input<<"[base+"<<(L*m)<<"];\n";
            source<<"mixed2 c2 = data"<<input<<"[base+"<<(2*L*m)<<"];\n";
            source<<"mixed2 c3 = data"<<input<<"[base+"<<(3*L*m)<<"];\n";
            source<<"mixed2 c4 = data"<<input<<"[base+"<<(4*L*m)<<"];\n";
            source<<"mixed2 c5 = data"<<input<<"[base+"<<(5*L*m)<<"];\n";
            source<<"mixed2 c6 = data"<<input<<"[base+"<<(6*L*m)<<"];\n";
            source<<"mixed2 d0 = c1+c6;\n";
            source<<"mixed2 d1 = c1-c6;\n";
            source<<"mixed2 d2 = c2+c5;\n";
            source<<"mixed2 d3 = c2-c5;\n";
            source<<"mixed2 d4 = c4+c3;\n";
            source<<"mixed2 d5 = c4-c3;\n";
            source<<"mixed2 d6 = d2+d0;\n";
            source<<"mixed2 d7 = d5+d3;\n";
            source<<"mixed2 b0 = c0+d6+d4;\n";
            source<<"mixed2 b1 = "<<cc.doubleToString((cos(2*M_PI/7)+cos(4*M_PI/7)+cos(6*M_PI/7))/3-1)<<"*(d6+d4);\n";
            source<<"mixed2 b2 = "<<cc.doubleToString((2*cos(2*M_PI/7)-cos(4*M_PI/7)-cos(6*M_PI/7))/3)<<"*(d0-d4);\n";
            source<<"mixed2 b3 = "<<cc.doubleToString((cos(2*M_PI/7)-2*cos(4*M_PI/7)+cos(6*M_PI/7))/3)<<"*(d4-d2);\n";
            source<<"mixed2 b4 = "<<cc.doubleToString((cos(2*M_PI/7)+cos(4*M_PI/7)-2*cos(6*M_PI/7))/3)<<"*(d2-d0);\n";
            source<<"mixed2 b5 = -("<<sign<<")*"<<cc.doubleToString((sin(2*M_PI/7)+sin(4*M_PI/7)-sin(6*M_PI/7))/3)<<"*(d7+d1);\n";
            source<<"mixed2 b6 = -("<<sign<<")*"<<cc.doubleToString((2*sin(2*M_PI/7)-sin(4*M_PI/7)+sin(6*M_PI/7))/3)<<"*(d1-d5);\n";
            source<<"mixed2 b7 = -("<<sign<<")*"<<cc.doubleToString((sin(2*M_PI/7)-2*sin(4*M_PI/7)-sin(6*M_PI/7))/3)<<"*(d5-d3);\n";
            source<<"mixed2 b8 = -("<<sign<<")*"<<cc.doubleToString((sin(2*M_PI/7)+sin(4*M_PI/7)+2*sin(6*M_PI/7))/3)<<"*(d3-d1);\n";
            source<<"mixed2 t0 = b0+b1;\n";
            source<<"mixed2 t1 = b2+b3;\n";
            source<<"mixed2 t2 = b4-b3;\n";
            source<<"mixed2 t3 = -b2-b4;\n";
            source<<"mixed2 t4 = b6+b7;\n";
            source<<"mixed2 t5 = b8-b7;\n";
            source<<"mixed2 t6 = -b8-b6;\n";
            source<<"mixed2 t7 = t0+t1;\n";
            source<<"mixed2 t8 = t0+t2;\n";
            source<<"mixed2 t9 = t0+t3;\n";
            source<<"mixed2 t10 = make_mixed2(t4.y+b5.y, -(t4.x+b5.x));\n";
            source<<"mixed2 t11 = make_mixed2(t5.y+b5.y, -(t5.x+b5.x));\n";
            source<<"mixed2 t12 = make_mixed2(t6.y+b5.y, -(t6.x+b5.x));\n";
            source<<"data"<<output<<"[base+6*j*"<<m<<"] = b0;\n";
            source<<"data"<<output<<"[base+(6*j+1)*"<<m<<"] = "<<mult<<"(w[j*"<<size<<"/"<<(7*L)<<"], t7-t10);\n";
            source<<"data"<<output<<"[base+(6*j+2)*"<<m<<"] = "<<mult<<"(w[j*"<<(2*size)<<"/"<<(7*L)<<"], t9-t12);\n";
            source<<"data"<<output<<"[base+(6*j+3)*"<<m<<"] = "<<mult<<"(w[j*"<<(3*size)<<"/"<<(7*L)<<"], t8+t11);\n";
            source<<"data"<<output<<"[base+(6*j+4)*"<<m<<"] = "<<mult<<"(w[j*"<<(4*size)<<"/"<<(7*L)<<"], t8-t11);\n";
            source<<"data"<<output<<"[base+(6*j+5)*"<<m<<"] = "<<mult<<"(w[j*"<<(5*size)<<"/"<<(7*L)<<"], t9+t12);\n";
            source<<"data"<<output<<"[base+(6*j+6)*"<<m<<"] = "<<mult<<"(w[j*"<<(6*size)<<"/"<<(7*L)<<"], t7+t10);\n";
        }
        else if (radix == 5) {
            source<<"mixed2 c0 = data"<<input<<"[base];\n";
            source<<"mixed2 c1 = data"<<input<<"[base+"<<(L*m)<<"];\n";
            source<<"mixed2 c2 = data"<<input<<"[base+"<<(2*L*m)<<"];\n";
            source<<"mixed2 c3 = data"<<input<<"[base+"<<(3*L*m)<<"];\n";
            source<<"mixed2 c4 = data"<<input<<"[base+"<<(4*L*m)<<"];\n";
            source<<"mixed2 d0 = c1+c4;\n";
            source<<"mixed2 d1 = c2+c3;\n";
            source<<"mixed2 d2 = "<<cc.doubleToString(sin(0.4*M_PI))<<"*(c1-c4);\n";
            source<<"mixed2 d3 = "<<cc.doubleToString(sin(0.4*M_PI))<<"*(c2-c3);\n";
            source<<"mixed2 d4 = d0+d1;\n";
            source<<"mixed2 d5 = "<<cc.doubleToString(0.25*sqrt(5.0))<<"*(d0-d1);\n";
            source<<"mixed2 d6 = c0-0.25f*d4;\n";
            source<<"mixed2 d7 = d6+d5;\n";
            source<<"mixed2 d8 = d6-d5;\n";
            string coeff = cc.doubleToString(sin(0.2*M_PI)/sin(0.4*M_PI));
            source<<"mixed2 d9 = "<<sign<<"*make_mixed2(d2.y+"<<coeff<<"*d3.y, -d2.x-"<<coeff<<"*d3.x);\n";
            source<<"mixed2 d10 = "<<sign<<"*make_mixed2("<<coeff<<"*d2.y-d3.y, d3.x-"<<coeff<<"*d2.x);\n";
            source<<"data"<<output<<"[base+4*j*"<<m<<"] = c0+d4;\n";
            source<<"data"<<output<<"[base+(4*j+1)*"<<m<<"] = "<<mult<<"(w[j*"<<size<<"/"<<(5*L)<<"], d7+d9);\n";
            source<<"data"<<output<<"[base+(4*j+2)*"<<m<<"] = "<<mult<<"(w[j*"<<(2*size)<<"/"<<(5*L)<<"], d8+d10);\n";
            source<<"data"<<output<<"[base+(4*j+3)*"<<m<<"] = "<<mult<<"(w[j*"<<(3*size)<<"/"<<(5*L)<<"], d8-d10);\n";
            source<<"data"<<output<<"[base+(4*j+4)*"<<m<<"] = "<<mult<<"(w[j*"<<(4*size)<<"/"<<(5*L)<<"], d7-d9);\n";
        }
        else if (radix == 4) {
            source<<"mixed2 c0 = data"<<input<<"[base];\n";
            source<<"mixed2 c1 = data"<<input<<"[base+"<<(L*m)<<"];\n";
            source<<"mixed2 c2 = data"<<input<<"[base+"<<(2*L*m)<<"];\n";
            source<<"mixed2 c3 = data"<<input<<"[base+"<<(3*L*m)<<"];\n";
            source<<"mixed2 d0 = c0+c2;\n";
            source<<"mixed2 d1 = c0-c2;\n";
            source<<"mixed2 d2 = c1+c3;\n";
            source<<"mixed2 d3 = "<<sign<<"*make_mixed2(c1.y-c3.y, c3.x-c1.x);\n";
            source<<"data"<<output<<"[base+3*j*"<<m<<"] = d0+d2;\n";
            source<<"data"<<output<<"[base+(3*j+1)*"<<m<<"] = "<<mult<<"(w[j*"<<size<<"/"<<(4*L)<<"], d1+d3);\n";
            source<<"data"<<output<<"[base+(3*j+2)*"<<m<<"] = "<<mult<<"(w[j*"<<(2*size)<<"/"<<(4*L)<<"], d0-d2);\n";
            source<<"data"<<output<<"[base+(3*j+3)*"<<m<<"] = "<<mult<<"(w[j*"<<(3*size)<<"/"<<(4*L)<<"], d1-d3);\n";
        }
        else if (radix == 3) {
            source<<"mixed2 c0 = data"<<input<<"[base];\n";
            source<<"mixed2 c1 = data"<<input<<"[base+"<<(L*m)<<"];\n";
            source<<"mixed2 c2 = data"<<input<<"[base+"<<(2*L*m)<<"];\n";
            source<<"mixed2 d0 = c1+c2;\n";
            source<<"mixed2 d1 = c0-0.5f*d0;\n";
            source<<"mixed2 d2 = "<<sign<<"*"<<cc.doubleToString(sin(M_PI/3.0))<<"*make_mixed2(c1.y-c2.y, c2.x-c1.x);\n";
            source<<"data"<<output<<"[base+2*j*"<<m<<"] = c0+d0;\n";
            source<<"data"<<output<<"[base+(2*j+1)*"<<m<<"] = "<<mult<<"(w[j*"<<size<<"/"<<(3*L)<<"], d1+d2);\n";
            source<<"data"<<output<<"[base+(2*j+2)*"<<m<<"] = "<<mult<<"(w[j*"<<(2*size)<<"/"<<(3*L)<<"], d1-d2);\n";
        }
        else if (radix == 2) {
            source<<"mixed2 c0 = data"<<input<<"[base];\n";
            source<<"mixed2 c1 = data"<<input<<"[base+"<<(L*m)<<"];\n";
            source<<"data"<<output<<"[base+j*"<<m<<"] = c0+c1;\n";
            source<<"data"<<output<<"[base+(j+1)*"<<m<<"] = "<<mult<<"(w[j*"<<size<<"/"<<(2*L)<<"], c0-c1);\n";
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

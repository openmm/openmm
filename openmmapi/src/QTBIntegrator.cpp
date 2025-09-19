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

#include "openmm/QTBIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <set>
#include <string>

using namespace OpenMM;
using std::map;
using std::set;
using std::string;
using std::vector;

QTBIntegrator::QTBIntegrator(double temperature, double frictionCoeff, double stepSize) {
    setTemperature(temperature);
    setFriction(frictionCoeff);
    setSegmentLength(1.0);
    setCutoffFrequency(500.0);
    setDefaultAdaptationRate(0.01);
    setStepSize(stepSize);
    setConstraintTolerance(1e-5);
    setRandomNumberSeed(0);
}

void QTBIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    context = &contextRef;
    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateQTBStepKernel::Name(), contextRef);
    kernel.getAs<IntegrateQTBStepKernel>().initialize(contextRef.getSystem(), *this);
}

void QTBIntegrator::setTemperature(double temp) {
    if (temp < 0)
        throw OpenMMException("Temperature cannot be negative");
    temperature = temp;
}

void QTBIntegrator::setFriction(double coeff) {
    if (coeff < 0)
        throw OpenMMException("Friction cannot be negative");
    friction = coeff;
}

double QTBIntegrator::getSegmentLength() const {
    return segmentLength;
}

void QTBIntegrator::setSegmentLength(double length) {
    segmentLength = length;
}

double QTBIntegrator::getCutoffFrequency() const {
    return cutoffFrequency;
}

void QTBIntegrator::setCutoffFrequency(double cutoff) {
    cutoffFrequency = cutoff;
}

const map<int, int>& QTBIntegrator::getParticleTypes() const {
    return particleType;
}

void QTBIntegrator::setParticleType(int index, int type) {
    particleType[index] = type;
}

double QTBIntegrator::getDefaultAdaptationRate() const {
    return defaultAdaptationRate;
}

void QTBIntegrator::setDefaultAdaptationRate(double rate) {
    defaultAdaptationRate = rate;
}

const map<int, double>& QTBIntegrator::getTypeAdaptationRates() const {
    return typeAdaptationRates;
}

void QTBIntegrator::setTypeAdaptationRate(int type, double rate) {
    typeAdaptationRates[type] = rate;
}

void QTBIntegrator::getAdaptedFriction(int particle, vector<double>& friction) const {
    kernel.getAs<IntegrateQTBStepKernel>().getAdaptedFriction(*context, particle, friction);
}

void QTBIntegrator::cleanup() {
    kernel = Kernel();
}

vector<string> QTBIntegrator::getKernelNames() {
    vector<std::string> names;
    names.push_back(IntegrateQTBStepKernel::Name());
    return names;
}

double QTBIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateQTBStepKernel>().computeKineticEnergy(*context, *this);
}

bool QTBIntegrator::kineticEnergyRequiresForce() const {
    return false;
}

void QTBIntegrator::step(int steps) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");
    for (int i = 0; i < steps; ++i) {
        context->updateContextState();
        context->calcForcesAndEnergy(true, false, getIntegrationForceGroups());
        kernel.getAs<IntegrateQTBStepKernel>().execute(*context, *this);
    }
}

void QTBIntegrator::createCheckpoint(std::ostream& stream) const {
    kernel.getAs<IntegrateQTBStepKernel>().createCheckpoint(*context, stream);
}

void QTBIntegrator::loadCheckpoint(std::istream& stream) {
    kernel.getAs<IntegrateQTBStepKernel>().loadCheckpoint(*context, stream);
}

void QTBIntegrator::serializeParameters(SerializationNode& node) const {
    node.setIntProperty("version", 1);
    vector<int> particles;
    set<int> types;
    for (int i = 0; i < context->getSystem().getNumParticles(); i++) {
        if (particleType.find(i) == particleType.end())
            particles.push_back(i);
        else if (types.find(particleType.at(i)) == types.end()) {
            particles.push_back(i);
            types.insert(particleType.at(i));
        }
    }
    vector<double> friction;
    for (int i : particles) {
        SerializationNode& particle = node.createChildNode("Friction").setIntProperty("particle", i);
        getAdaptedFriction(i, friction);
        for (double f : friction)
            particle.createChildNode("F").setDoubleProperty("v", f);
    }
}

void QTBIntegrator::deserializeParameters(const SerializationNode& node) {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    for (const SerializationNode& child : node.getChildren()) {
        int particle = child.getIntProperty("particle");
        vector<double> friction;
        for (const SerializationNode& f : child.getChildren())
            friction.push_back(f.getDoubleProperty("v"));
        kernel.getAs<IntegrateQTBStepKernel>().setAdaptedFriction(*context, particle, friction);
    }
}

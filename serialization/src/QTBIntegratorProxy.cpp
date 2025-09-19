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

#include "openmm/serialization/QTBIntegratorProxy.h"
#include <OpenMM.h>

using namespace std;
using namespace OpenMM;

QTBIntegratorProxy::QTBIntegratorProxy() : SerializationProxy("QTBIntegrator") {
}

void QTBIntegratorProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const QTBIntegrator& integrator = *reinterpret_cast<const QTBIntegrator*>(object);
    node.setDoubleProperty("stepSize", integrator.getStepSize());
    node.setDoubleProperty("constraintTolerance", integrator.getConstraintTolerance());
    node.setDoubleProperty("temperature", integrator.getTemperature());
    node.setDoubleProperty("friction", integrator.getFriction());
    node.setDoubleProperty("segmentLength", integrator.getSegmentLength());
    node.setDoubleProperty("cutoffFrequency", integrator.getCutoffFrequency());
    node.setDoubleProperty("defaultAdaptationRate", integrator.getDefaultAdaptationRate());
    node.setIntProperty("randomSeed", integrator.getRandomNumberSeed());
    node.setIntProperty("integrationForceGroups", integrator.getIntegrationForceGroups());
    SerializationNode& particleTypes = node.createChildNode("ParticleTypes");
    for (auto type : integrator.getParticleTypes())
        particleTypes.createChildNode("Particle").setIntProperty("index", type.first).setIntProperty("type", type.second);
    SerializationNode& adaptationRates = node.createChildNode("TypeAdaptationRates");
    for (auto type : integrator.getTypeAdaptationRates())
        adaptationRates.createChildNode("Type").setIntProperty("index", type.first).setDoubleProperty("rate", type.second);
}

void* QTBIntegratorProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    QTBIntegrator *integrator = new QTBIntegrator(node.getDoubleProperty("temperature"),
                                                  node.getDoubleProperty("friction"),
                                                  node.getDoubleProperty("stepSize"));
    integrator->setConstraintTolerance(node.getDoubleProperty("constraintTolerance"));
    integrator->setSegmentLength(node.getDoubleProperty("segmentLength"));
    integrator->setCutoffFrequency(node.getDoubleProperty("cutoffFrequency"));
    integrator->setDefaultAdaptationRate(node.getDoubleProperty("defaultAdaptationRate"));
    integrator->setRandomNumberSeed(node.getIntProperty("randomSeed"));
    integrator->setIntegrationForceGroups(node.getIntProperty("integrationForceGroups", 0xFFFFFFFF));
    const SerializationNode& particleTypes = node.getChildNode("ParticleTypes");
    for (auto& particle : particleTypes.getChildren())
        integrator->setParticleType(particle.getIntProperty("index"), particle.getIntProperty("type"));
    const SerializationNode& adaptationRates = node.getChildNode("TypeAdaptationRates");
    for (auto& type : adaptationRates.getChildren())
        integrator->setTypeAdaptationRate(type.getIntProperty("index"), type.getDoubleProperty("rate"));
    return integrator;
}
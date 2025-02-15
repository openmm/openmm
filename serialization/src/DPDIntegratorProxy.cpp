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

#include "openmm/serialization/DPDIntegratorProxy.h"
#include <OpenMM.h>

using namespace std;
using namespace OpenMM;

DPDIntegratorProxy::DPDIntegratorProxy() : SerializationProxy("DPDIntegrator") {

}

void DPDIntegratorProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const DPDIntegrator& integrator = *reinterpret_cast<const DPDIntegrator*>(object);
    node.setDoubleProperty("stepSize", integrator.getStepSize());
    node.setDoubleProperty("constraintTolerance", integrator.getConstraintTolerance());
    node.setDoubleProperty("temperature", integrator.getTemperature());
    node.setDoubleProperty("defaultFriction", integrator.getDefaultFriction());
    node.setDoubleProperty("defaultCutoff", integrator.getDefaultCutoff());
    node.setIntProperty("randomSeed", integrator.getRandomNumberSeed());
    node.setIntProperty("integrationForceGroups", integrator.getIntegrationForceGroups());
    SerializationNode& types = node.createChildNode("ParticleTypes");
    for (auto type: integrator.getParticleTypes())
        types.createChildNode("Type").setIntProperty("particle", type.first).setIntProperty("type", type.second);
    SerializationNode& pairs = node.createChildNode("TypePairs");
    for (int i = 0; i < integrator.getNumTypePairs(); i++) {
        int type1, type2;
        double friction, cutoff;
        integrator.getTypePairParameters(i, type1, type2, friction, cutoff);
        pairs.createChildNode("Pair").setIntProperty("type1", type1).setIntProperty("type2", type2).setDoubleProperty("friction", friction).setDoubleProperty("cutoff", cutoff);
    }
}

void* DPDIntegratorProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    DPDIntegrator *integrator = new DPDIntegrator(node.getDoubleProperty("temperature"),
                                                            node.getDoubleProperty("defaultFriction"),
                                                            node.getDoubleProperty("defaultCutoff"),
                                                            node.getDoubleProperty("stepSize"));
    integrator->setConstraintTolerance(node.getDoubleProperty("constraintTolerance"));
    integrator->setRandomNumberSeed(node.getIntProperty("randomSeed"));
    integrator->setIntegrationForceGroups(node.getIntProperty("integrationForceGroups", 0xFFFFFFFF));
    const SerializationNode& types = node.getChildNode("ParticleTypes");
    for (auto& type: types.getChildren())
        integrator->setParticleType(type.getIntProperty("particle"), type.getIntProperty("type"));
    const SerializationNode& pairs = node.getChildNode("TypePairs");
    for (auto& pair: pairs.getChildren())
        integrator->addTypePair(pair.getIntProperty("type1"), pair.getIntProperty("type2"), pair.getDoubleProperty("friction"), pair.getDoubleProperty("cutoff"));
    return integrator;
}

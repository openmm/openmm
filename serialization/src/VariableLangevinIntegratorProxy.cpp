/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2024 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Yutong Zhao                                        *
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

#include "openmm/serialization/VariableLangevinIntegratorProxy.h"
#include <OpenMM.h>

using namespace std;
using namespace OpenMM;

VariableLangevinIntegratorProxy::VariableLangevinIntegratorProxy() : SerializationProxy("VariableLangevinIntegrator") {
}

void VariableLangevinIntegratorProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 2);
    const VariableLangevinIntegrator& integrator = *reinterpret_cast<const VariableLangevinIntegrator*>(object);
    node.setDoubleProperty("stepSize", integrator.getStepSize());
    node.setDoubleProperty("constraintTolerance", integrator.getConstraintTolerance());
    node.setDoubleProperty("temperature", integrator.getTemperature());
    node.setDoubleProperty("friction", integrator.getFriction());
    node.setDoubleProperty("errorTol", integrator.getErrorTolerance());
    node.setDoubleProperty("maxStepSize", integrator.getMaximumStepSize());
    node.setIntProperty("randomSeed", integrator.getRandomNumberSeed());
    node.setIntProperty("integrationForceGroups", integrator.getIntegrationForceGroups());
}

void* VariableLangevinIntegratorProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version < 1 || version > 2)
        throw OpenMMException("Unsupported version number");
    VariableLangevinIntegrator *integrator = new VariableLangevinIntegrator(node.getDoubleProperty("temperature"),
                                                            node.getDoubleProperty("friction"),
                                                            node.getDoubleProperty("errorTol"));
    integrator->setStepSize(node.getDoubleProperty("stepSize"));
    integrator->setConstraintTolerance(node.getDoubleProperty("constraintTolerance"));
    integrator->setRandomNumberSeed(node.getIntProperty("randomSeed"));
    integrator->setIntegrationForceGroups(node.getIntProperty("integrationForceGroups", 0xFFFFFFFF));
    if (version > 1)
        integrator->setMaximumStepSize(node.getDoubleProperty("maxStepSize"));
    return integrator;
}
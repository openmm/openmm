/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2016 Stanford University and the Authors.      *
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

#include "openmm/serialization/MonteCarloBarostatProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/MonteCarloBarostat.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

MonteCarloBarostatProxy::MonteCarloBarostatProxy() : SerializationProxy("MonteCarloBarostat") {
}

void MonteCarloBarostatProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const MonteCarloBarostat& force = *reinterpret_cast<const MonteCarloBarostat*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setDoubleProperty("pressure", force.getDefaultPressure());
    node.setDoubleProperty("temperature", force.getDefaultTemperature());
    node.setIntProperty("frequency", force.getFrequency());
    node.setIntProperty("randomSeed", force.getRandomNumberSeed());
}

void* MonteCarloBarostatProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    MonteCarloBarostat* force = NULL;
    try {
        force = new MonteCarloBarostat(node.getDoubleProperty("pressure"), node.getDoubleProperty("temperature"), node.getIntProperty("frequency"));
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setRandomNumberSeed(node.getIntProperty("randomSeed"));
        return force;
    }
    catch (...) {
        if (force != NULL)
            delete force;
        throw;
    }
}

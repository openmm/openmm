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

#include "openmm/serialization/MonteCarloMembraneBarostatProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/MonteCarloMembraneBarostat.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

MonteCarloMembraneBarostatProxy::MonteCarloMembraneBarostatProxy() : SerializationProxy("MonteCarloMembraneBarostat") {
}

void MonteCarloMembraneBarostatProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const MonteCarloMembraneBarostat& force = *reinterpret_cast<const MonteCarloMembraneBarostat*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setDoubleProperty("pressure", force.getDefaultPressure());
    node.setDoubleProperty("surfaceTension", force.getDefaultSurfaceTension());
    node.setDoubleProperty("temperature", force.getDefaultTemperature());
    node.setIntProperty("xymode", force.getXYMode());
    node.setIntProperty("zmode", force.getZMode());
    node.setIntProperty("frequency", force.getFrequency());
    node.setIntProperty("randomSeed", force.getRandomNumberSeed());
}

void* MonteCarloMembraneBarostatProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    MonteCarloMembraneBarostat* force = NULL;
    try {
        MonteCarloMembraneBarostat::XYMode xymode = (MonteCarloMembraneBarostat::XYMode) node.getIntProperty("xymode");
        MonteCarloMembraneBarostat::ZMode zmode = (MonteCarloMembraneBarostat::ZMode) node.getIntProperty("zmode");
        force = new MonteCarloMembraneBarostat(node.getDoubleProperty("pressure"), node.getDoubleProperty("surfaceTension"),
                node.getDoubleProperty("temperature"), xymode, zmode, node.getIntProperty("frequency"));
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

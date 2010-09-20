/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

#include "openmm/serialization/HarmonicAngleForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/HarmonicAngleForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

HarmonicAngleForceProxy::HarmonicAngleForceProxy() : SerializationProxy("HarmonicAngleForce") {
}

void HarmonicAngleForceProxy::serialize(const void* object, SerializationNode& node) const {
    const HarmonicAngleForce& force = *reinterpret_cast<const HarmonicAngleForce*>(object);
    SerializationNode& bonds = node.createChildNode("Angles");
    for (int i = 0; i < force.getNumAngles(); i++) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(i, particle1, particle2, particle3, angle, k);
        bonds.createChildNode("Angle").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setIntProperty("p3", particle3).setDoubleProperty("a", angle).setDoubleProperty("k", k);
    }
}

void* HarmonicAngleForceProxy::deserialize(const SerializationNode& node) const {
    HarmonicAngleForce* force = new HarmonicAngleForce();
    try {
        const SerializationNode& bonds = node.getChildNode("Angles");
        for (int i = 0; i < (int) bonds.getChildren().size(); i++) {
            const SerializationNode& constraint = bonds.getChildren()[i];
            force->addAngle(constraint.getDoubleProperty("p1"), constraint.getDoubleProperty("p2"), constraint.getDoubleProperty("p3"), constraint.getDoubleProperty("a"), constraint.getDoubleProperty("k"));
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}


/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
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

#include "openmm/serialization/AmoebaStretchBendForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/AmoebaStretchBendForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

AmoebaStretchBendForceProxy::AmoebaStretchBendForceProxy() : SerializationProxy("AmoebaStretchBendForce") {
}

void AmoebaStretchBendForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 4);
    const AmoebaStretchBendForce& force = *reinterpret_cast<const AmoebaStretchBendForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setBoolProperty("usesPeriodic", force.usesPeriodicBoundaryConditions());
    SerializationNode& bonds = node.createChildNode("StretchBendAngles");
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force.getNumStretchBends()); ii++) {
        int particle1, particle2, particle3;
        double distanceAB, distanceCB, angle, k1, k2;
        force.getStretchBendParameters(ii, particle1, particle2, particle3, distanceAB, distanceCB, angle, k1, k2);
        bonds.createChildNode("StretchBendAngle").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setIntProperty("p3", particle3).setDoubleProperty("dAB", distanceAB).setDoubleProperty("dCB", distanceCB).setDoubleProperty("angle", angle).setDoubleProperty("k1", k1).setDoubleProperty("k2", k2);
    }

}

void* AmoebaStretchBendForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version < 1 || version > 4)
        throw OpenMMException("Unsupported version number");
    AmoebaStretchBendForce* force = new AmoebaStretchBendForce();
    try {
        if (version > 2)
            force->setForceGroup(node.getIntProperty("forceGroup", 0));
        if (version > 3)
            force->setUsesPeriodicBoundaryConditions(node.getBoolProperty("usesPeriodic"));
        const SerializationNode& bonds = node.getChildNode("StretchBendAngles");
        for (auto& bond : bonds.getChildren()) {
            double k1, k2;
            if (version == 1)
                k1 = k2 = bond.getDoubleProperty("k");
            else {
                k1 = bond.getDoubleProperty("k1");
                k2 = bond.getDoubleProperty("k2");
            }
            force->addStretchBend(bond.getIntProperty("p1"), bond.getIntProperty("p2"), bond.getIntProperty("p3"),  bond.getDoubleProperty("dAB"), bond.getDoubleProperty("dCB"), bond.getDoubleProperty("angle"), k1, k2);

        }
    }
    catch (...) {
        delete force;
        throw;
    }

    return force;
}

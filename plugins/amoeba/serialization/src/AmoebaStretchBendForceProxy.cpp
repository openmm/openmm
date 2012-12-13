/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
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
    node.setIntProperty("version", 1);
    const AmoebaStretchBendForce& force = *reinterpret_cast<const AmoebaStretchBendForce*>(object);
    SerializationNode& bonds = node.createChildNode("StretchBendAngles");
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force.getNumStretchBends()); ii++) {
        int particle1, particle2, particle3;
        double distanceAB, distanceCB, angle, k;
        force.getStretchBendParameters(ii, particle1, particle2, particle3, distanceAB, distanceCB, angle, k);
        bonds.createChildNode("StretchBendAngle").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setIntProperty("p3", particle3).setDoubleProperty("dAB", distanceAB).setDoubleProperty("dCB", distanceCB).setDoubleProperty("angle", angle).setDoubleProperty("k", k);
    }

}

void* AmoebaStretchBendForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    AmoebaStretchBendForce* force = new AmoebaStretchBendForce();
    try {
        const SerializationNode& bonds = node.getChildNode("StretchBendAngles");
        for ( unsigned int ii = 0; ii < (int) bonds.getChildren().size(); ii++) {
            const SerializationNode& bond = bonds.getChildren()[ii];
            force->addStretchBend(bond.getIntProperty("p1"), bond.getIntProperty("p2"), bond.getIntProperty("p3"),  bond.getDoubleProperty("dAB"), bond.getDoubleProperty("dCB"), bond.getDoubleProperty("angle"), bond.getDoubleProperty("k"));

        }
    }
    catch (...) {
        delete force;
        throw;
    }

    return force;
}

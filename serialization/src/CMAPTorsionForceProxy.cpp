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

#include "openmm/serialization/CMAPTorsionForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/CMAPTorsionForce.h"
#include <vector>

using namespace OpenMM;
using namespace std;

CMAPTorsionForceProxy::CMAPTorsionForceProxy() : SerializationProxy("CMAPTorsionForce") {
}

void CMAPTorsionForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 2);
    const CMAPTorsionForce& force = *reinterpret_cast<const CMAPTorsionForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setBoolProperty("usesPeriodic", force.usesPeriodicBoundaryConditions());
    SerializationNode& maps = node.createChildNode("Maps");
    for (int i = 0; i < force.getNumMaps(); i++) {
        int size;
        vector<double> energy;
        force.getMapParameters(i, size, energy);
        SerializationNode& map = maps.createChildNode("Map").setIntProperty("size", size);
        for (auto e : energy)
            map.createChildNode("Energy").setDoubleProperty("e", e);
    }
    SerializationNode& torsions = node.createChildNode("Torsions");
    for (int i = 0; i < force.getNumTorsions(); i++) {
        int map, a1, a2, a3, a4, b1, b2, b3, b4;
        force.getTorsionParameters(i, map, a1, a2, a3, a4, b1, b2, b3, b4);
        torsions.createChildNode("Torsion").setIntProperty("a1", a1).setIntProperty("a2", a2).setIntProperty("a3", a3).setIntProperty("a4", a4).setIntProperty("b1", b1).setIntProperty("b2", b2).setIntProperty("b3", b3).setIntProperty("b4", b4).setIntProperty("map", map);
    }
}

void* CMAPTorsionForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version < 1 || version > 2)
        throw OpenMMException("Unsupported version number");
    CMAPTorsionForce* force = new CMAPTorsionForce();
    try {
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        if (version > 1)
            force->setUsesPeriodicBoundaryConditions(node.getBoolProperty("usesPeriodic"));
        const SerializationNode& maps = node.getChildNode("Maps");
        for (auto& map : maps.getChildren()) {
            int size = map.getIntProperty("size");
            if (size*size != map.getChildren().size())
                throw OpenMMException("Wrong number of values specified for CMAP");
            vector<double> energy(size*size);
            for (int j = 0; j < (int) energy.size(); j++)
                energy[j] = map.getChildren()[j].getDoubleProperty("e");
            force->addMap(size, energy);
        }
        const SerializationNode& torsions = node.getChildNode("Torsions");
        for (auto& torsion : torsions.getChildren())
            force->addTorsion(torsion.getIntProperty("map"), torsion.getIntProperty("a1"), torsion.getIntProperty("a2"), torsion.getIntProperty("a3"), torsion.getIntProperty("a4"),
                    torsion.getIntProperty("b1"), torsion.getIntProperty("b2"), torsion.getIntProperty("b3"), torsion.getIntProperty("b4"));
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}

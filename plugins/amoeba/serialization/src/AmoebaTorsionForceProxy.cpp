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

#include "openmm/serialization/AmoebaTorsionForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/AmoebaTorsionForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

AmoebaTorsionForceProxy::AmoebaTorsionForceProxy() : SerializationProxy("AmoebaTorsionForce") {
}

static void addTorsionValues( SerializationNode& torsion, std::vector<double>& torsionValues ){
    for (int j = 0; j < torsionValues.size(); j++) {
        torsion.createChildNode("Value").setDoubleProperty("v", torsionValues[j]);
    }
}

static void loadTorsionValues( SerializationNode& torsion, std::vector<double>& torsionValues ){
    for (int j = 0; j < torsionValues.size(); j++) {
        torsion.createChildNode("Value").setDoubleProperty("v", torsionValues[j]);
    }
}

void AmoebaTorsionForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const AmoebaTorsionForce& force = *reinterpret_cast<const AmoebaTorsionForce*>(object);
    SerializationNode& bonds = node.createChildNode("Torsion");
    for (int i = 0; i < force.getNumTorsions(); i++) {
        int particle1, particle2, particle3, particle4;
        std::vector<double> torsion1;
        std::vector<double> torsion2;
        std::vector<double> torsion3;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, torsion1, torsion2, torsion3);
        bonds.createChildNode("Torsion").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setIntProperty("p3", particle3).setIntProperty("p4", particle4);
        SerializationNode& torsion = bonds.createChildNode("Torsion1").setIntProperty("size", torsion1.size());
        addTorsionValues( torsion, torsion1 );
        torsion = bonds.createChildNode("Torsion2").setIntProperty("size", torsion2.size());
        addTorsionValues( torsion, torsion2 );
        torsion = bonds.createChildNode("Torsion3").setIntProperty("size", torsion3.size());
        addTorsionValues( torsion, torsion3 );
    }
}

void* AmoebaTorsionForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    AmoebaTorsionForce* force = new AmoebaTorsionForce();
/*
    try {
        const SerializationNode& bonds = node.getChildNode("Torsion");
        for (int i = 0; i < (int) bonds.getChildren().size(); i++) {
            const SerializationNode& bond = bonds.getChildren()[i];
            //force->addTorsion(bond.getIntProperty("p1"), bond.getIntProperty("p2"), bond.getIntProperty("p3"),  bond.getIntProperty("p4") );
            std::vector<double> torsion1, std::vector<double> torsion2, std::vector<double> torsion3;
            const SerializationNode& torsion = bond.getChildNode("Torsion1");
        }
    }
    catch (...) {
        delete force;
        throw;
    }
*/
    return force;
}

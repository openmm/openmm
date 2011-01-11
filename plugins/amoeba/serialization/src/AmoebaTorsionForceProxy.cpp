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

static void addTorsionValues( SerializationNode& torsion, const std::vector<double>& torsionValues ){
    for (unsigned int jj = 0; jj < torsionValues.size(); jj++) {
        torsion.createChildNode("Value").setDoubleProperty("v", torsionValues[jj]);
    }
}

static void loadTorsionValues( SerializationNode& torsion, std::vector<double>& torsionValues ){
    int size = torsion.getIntProperty("size");
    torsionValues.resize(size);
    for (unsigned int jj = 0; jj < size; jj++) {
        torsionValues[jj] = ( torsion.getChildren()[jj].getDoubleProperty("v") );
    }
}

void AmoebaTorsionForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const AmoebaTorsionForce& force = *reinterpret_cast<const AmoebaTorsionForce*>(object);

    SerializationNode& bonds = node.createChildNode("Torsion").setIntProperty("size",  force.getNumTorsions() );
    for (unsigned int ii = 0; ii < force.getNumTorsions(); ii++) {

        int particle1, particle2, particle3, particle4;
        std::vector<double> torsion1;
        std::vector<double> torsion2;
        std::vector<double> torsion3;
        force.getTorsionParameters(ii, particle1, particle2, particle3, particle4, torsion1, torsion2, torsion3);

        SerializationNode& torsionBond        = bonds.createChildNode("Torsion");
        torsionBond.setIntProperty("p1", particle1).setIntProperty("p2", particle2).setIntProperty("p3", particle3).setIntProperty("p4", particle4);

        SerializationNode& torsionVector1     = torsionBond.createChildNode("Torsion1").setIntProperty("size", torsion1.size());
        addTorsionValues( torsionVector1, torsion1 );

        SerializationNode& torsionVector2     = torsionBond.createChildNode("Torsion2").setIntProperty("size", torsion2.size());
        addTorsionValues( torsionVector2, torsion2 );

        SerializationNode& torsionVector3     = torsionBond.createChildNode("Torsion3").setIntProperty("size", torsion3.size());
        addTorsionValues( torsionVector3, torsion3 );
    }
}

void* AmoebaTorsionForceProxy::deserialize(const SerializationNode& node) const {

    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");

    AmoebaTorsionForce* force = new AmoebaTorsionForce();
    try {
        const SerializationNode& bonds     = node.getChildNode("Torsion");
        vector<SerializationNode> children = bonds.getChildren();
        for (unsigned int i = 0; i < children.size(); i++) {
            SerializationNode& bond = children[i];

            std::vector<double> torsion1;
            std::vector<double> torsion2;
            std::vector<double> torsion3;

            SerializationNode& torsionNode1 = bond.getChildNode("Torsion1");
            loadTorsionValues( torsionNode1, torsion1 );

            SerializationNode& torsionNode2 = bond.getChildNode("Torsion2");
            loadTorsionValues( torsionNode2, torsion2 );

            SerializationNode& torsionNode3 = bond.getChildNode("Torsion3");
            loadTorsionValues( torsionNode3, torsion3 );

            force->addTorsion(bond.getIntProperty("p1"), bond.getIntProperty("p2"), bond.getIntProperty("p3"),  bond.getIntProperty("p4"), torsion1, torsion2, torsion3 );
        }
    }
    catch (...) {
        delete force;
        throw;
    }

    return force;
}

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

#include "openmm/serialization/AmoebaTorsionTorsionForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/AmoebaTorsionTorsionForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

AmoebaTorsionTorsionForceProxy::AmoebaTorsionTorsionForceProxy() : SerializationProxy("AmoebaTorsionTorsionForce") {
}

static void loadGrid(const SerializationNode& grid, std::vector< std::vector< std::vector<double> > >& gridVector) {

    const std::vector<SerializationNode>& gridSerializationRows  = grid.getChildren();
    gridVector.resize(gridSerializationRows.size());

    for (unsigned int ii = 0; ii < gridSerializationRows.size(); ii++) {
        const std::vector<SerializationNode>& gridSerializationColumns  = gridSerializationRows[ii].getChildren();
        gridVector[ii].resize(gridSerializationColumns.size());
        for (unsigned int jj = 0; jj < gridSerializationColumns.size(); jj++) {
            const SerializationNode& gridSerializationColumnNode = gridSerializationColumns[jj];
            gridVector[ii][jj].resize(6);
            gridVector[ii][jj][0] = gridSerializationColumnNode.getDoubleProperty("x");
            gridVector[ii][jj][1] = gridSerializationColumnNode.getDoubleProperty("y");
            gridVector[ii][jj][2] = gridSerializationColumnNode.getDoubleProperty("f");
            gridVector[ii][jj][3] = gridSerializationColumnNode.getDoubleProperty("fx");
            gridVector[ii][jj][4] = gridSerializationColumnNode.getDoubleProperty("fy");
            gridVector[ii][jj][5] = gridSerializationColumnNode.getDoubleProperty("fxy");
        }
    }
}

void AmoebaTorsionTorsionForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 3);
    const AmoebaTorsionTorsionForce& force = *reinterpret_cast<const AmoebaTorsionTorsionForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setBoolProperty("usesPeriodic", force.usesPeriodicBoundaryConditions());

    // grid[xIdx][yIdx][6 values]

    // value0 = x-Grid value
    // value1 = y-Grid value
    // value2 = F          function value
    // value3 = F_x        partial f wrt x
    // value4 = F_y        partial f wrt y
    // value5 = F_xy       partial f wrt x,y

    SerializationNode& grids = node.createChildNode("TorsionTorsionGrids");
    for (unsigned int kk = 0; kk < static_cast<unsigned int>(force.getNumTorsionTorsionGrids()); kk++) {

        const std::vector< std::vector< std::vector<double> > > grid = force.getTorsionTorsionGrid(kk);

        unsigned int gridCount = 0;
        unsigned int gridYsize =  grid[0].size();
        for (unsigned int ii = 0; ii < grid.size(); ii++) {
            gridCount += grid[ii].size();
        }

        SerializationNode& gridNode = grids.createChildNode("TorsionTorsionGrid");
        for (unsigned int ii = 0; ii < grid.size(); ii++) {
            SerializationNode& gridSerializationRow = gridNode.createChildNode("RowNode");
            gridSerializationRow.setIntProperty("dim", ii);
            for (unsigned int jj = 0; jj < grid[ii].size(); jj++) {
                SerializationNode& gridSerializationColumnNode = gridSerializationRow.createChildNode("ColumnNode");
                gridSerializationColumnNode.setIntProperty("dim", jj);
                unsigned int index = 0;
                gridSerializationColumnNode.setDoubleProperty("x",   grid[ii][jj][index++]);
                gridSerializationColumnNode.setDoubleProperty("y",   grid[ii][jj][index++]);
                gridSerializationColumnNode.setDoubleProperty("f",   grid[ii][jj][index++]);
                gridSerializationColumnNode.setDoubleProperty("fx",  grid[ii][jj][index++]);
                gridSerializationColumnNode.setDoubleProperty("fy",  grid[ii][jj][index++]);
                gridSerializationColumnNode.setDoubleProperty("fxy", grid[ii][jj][index++]);
            }
        }
    }

    SerializationNode& bonds = node.createChildNode("TorsionTorsion");
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force.getNumTorsionTorsions()); ii++) {
        int particle1, particle2, particle3, particle4, particle5;
        int chiralCheckAtomIndex, gridIndex;
        force.getTorsionTorsionParameters(ii, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
        bonds.createChildNode("TorsionTorsion").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setIntProperty("p3", particle3).setIntProperty("p4", particle4).setIntProperty("p5", particle5).setIntProperty("chiralCheckAtomIndex", chiralCheckAtomIndex).setIntProperty("gridIndex", gridIndex);

    }
}

void* AmoebaTorsionTorsionForceProxy::deserialize(const SerializationNode& node) const {

    int version = node.getIntProperty("version");
    if (version < 1 || version > 3)
        throw OpenMMException("Unsupported version number");

    AmoebaTorsionTorsionForce* force = new AmoebaTorsionTorsionForce();
    try {
        if (version > 1)
            force->setForceGroup(node.getIntProperty("forceGroup", 0));
        if (version > 2)
            force->setUsesPeriodicBoundaryConditions(node.getBoolProperty("usesPeriodic"));
        const SerializationNode& grids                    = node.getChildNode("TorsionTorsionGrids");
        const std::vector<SerializationNode>& gridList    = grids.getChildren();
        for (unsigned int ii = 0; ii < gridList.size(); ii++) {
            std::vector< std::vector< std::vector<double> > > gridVector;
            loadGrid(gridList[ii], gridVector);
            force->setTorsionTorsionGrid(ii, gridVector);
        }

        const SerializationNode& bonds     = node.getChildNode("TorsionTorsion");
        vector<SerializationNode> children = bonds.getChildren();
        for (unsigned int i = 0; i < children.size(); i++) {
            SerializationNode& bond = children[i];
            force->addTorsionTorsion(bond.getIntProperty("p1"), bond.getIntProperty("p2"), bond.getIntProperty("p3"),  bond.getIntProperty("p4"), bond.getIntProperty("p5"), bond.getIntProperty("chiralCheckAtomIndex"), bond.getIntProperty("gridIndex"));
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}

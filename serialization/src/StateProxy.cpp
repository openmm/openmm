/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2015 Stanford University and the Authors.      *
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

#include "openmm/serialization/StateProxy.h"
#include "openmm/Platform.h"
#include "openmm/State.h"
#include "openmm/Vec3.h"
#include <map>

using namespace std;
using namespace OpenMM;

StateProxy::StateProxy() : SerializationProxy("State") {

}

void StateProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    node.setStringProperty("openmmVersion", Platform::getOpenMMVersion());
    const State& s = *reinterpret_cast<const State*>(object);
    node.setDoubleProperty("time", s.getTime());
    Vec3 a,b,c;
    s.getPeriodicBoxVectors(a,b,c);
    SerializationNode& boxVectorsNode = node.createChildNode("PeriodicBoxVectors");
    boxVectorsNode.createChildNode("A").setDoubleProperty("x", a[0]).setDoubleProperty("y", a[1]).setDoubleProperty("z", a[2]);
    boxVectorsNode.createChildNode("B").setDoubleProperty("x", b[0]).setDoubleProperty("y", b[1]).setDoubleProperty("z", b[2]);
    boxVectorsNode.createChildNode("C").setDoubleProperty("x", c[0]).setDoubleProperty("y", c[1]).setDoubleProperty("z", c[2]);
    if ((s.getDataTypes()&State::Parameters) != 0) {
        s.getParameters();
        SerializationNode& parametersNode = node.createChildNode("Parameters");
        for (auto& param : s.getParameters())
            parametersNode.setDoubleProperty(param.first, param.second);
    }
    if ((s.getDataTypes()&State::Energy) != 0) {
        s.getPotentialEnergy();
        SerializationNode& energiesNode = node.createChildNode("Energies");
        energiesNode.setDoubleProperty("PotentialEnergy", s.getPotentialEnergy());
        energiesNode.setDoubleProperty("KineticEnergy", s.getKineticEnergy());
    }
    if ((s.getDataTypes()&State::Positions) != 0) {
        s.getPositions();
        SerializationNode& positionsNode = node.createChildNode("Positions");
        vector<Vec3> statePositions = s.getPositions();
        for (int i=0; i<statePositions.size();i++) {
           positionsNode.createChildNode("Position").setDoubleProperty("x", statePositions[i][0]).setDoubleProperty("y", statePositions[i][1]).setDoubleProperty("z", statePositions[i][2]);
        }
    }
    if ((s.getDataTypes()&State::Velocities) != 0) {
        s.getVelocities();
        SerializationNode& velocitiesNode = node.createChildNode("Velocities");
        vector<Vec3> stateVelocities = s.getVelocities();
        for (int i=0; i<stateVelocities.size();i++) {
           velocitiesNode.createChildNode("Velocity").setDoubleProperty("x", stateVelocities[i][0]).setDoubleProperty("y", stateVelocities[i][1]).setDoubleProperty("z", stateVelocities[i][2]);
        }
    }
    if ((s.getDataTypes()&State::Forces) != 0) {
        s.getForces();
        SerializationNode& forcesNode = node.createChildNode("Forces");
        vector<Vec3> stateForces = s.getForces();
        for (int i=0; i<stateForces.size();i++) {
            forcesNode.createChildNode("Force").setDoubleProperty("x", stateForces[i][0]).setDoubleProperty("y", stateForces[i][1]).setDoubleProperty("z", stateForces[i][2]);
        }
    }
}

void* StateProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    double outTime = node.getDoubleProperty("time");
    const SerializationNode& boxVectorsNode = node.getChildNode("PeriodicBoxVectors");
    const SerializationNode& AVec = boxVectorsNode.getChildNode("A");
    Vec3 outAVec(AVec.getDoubleProperty("x"),AVec.getDoubleProperty("y"),AVec.getDoubleProperty("z"));
    const SerializationNode& BVec = boxVectorsNode.getChildNode("B");
    Vec3 outBVec(BVec.getDoubleProperty("x"),BVec.getDoubleProperty("y"),BVec.getDoubleProperty("z"));
    const SerializationNode& CVec = boxVectorsNode.getChildNode("C");
    Vec3 outCVec(CVec.getDoubleProperty("x"),CVec.getDoubleProperty("y"),CVec.getDoubleProperty("z"));
    int types = 0;
    vector<int> arraySizes;
    State::StateBuilder builder(outTime);
    for (auto& child : node.getChildren()) {
        if (child.getName() == "Parameters") {
            map<string, double> outStateParams;
            for (auto& param : child.getProperties())
                outStateParams[param.first] = child.getDoubleProperty(param.first);
            builder.setParameters(outStateParams);
        }
        else if (child.getName() == "Energies") {
            double potentialEnergy = child.getDoubleProperty("PotentialEnergy");
            double kineticEnergy = child.getDoubleProperty("KineticEnergy");
            builder.setEnergy(kineticEnergy, potentialEnergy);
        }
        else if (child.getName() == "Positions") {
            vector<Vec3> outPositions;
            for (auto& particle : child.getChildren())
                outPositions.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
            builder.setPositions(outPositions);
            arraySizes.push_back(outPositions.size());
        }
        else if (child.getName() == "Velocities") {
            vector<Vec3> outVelocities;
            for (auto& particle : child.getChildren())
                outVelocities.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
            builder.setVelocities(outVelocities);
            arraySizes.push_back(outVelocities.size());
        }
        else if (child.getName() == "Forces") {
            vector<Vec3> outForces;
            for (auto& particle : child.getChildren())
                outForces.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
            builder.setForces(outForces);
            arraySizes.push_back(outForces.size());
        }
    }
    for (int i = 1; i < arraySizes.size(); i++) {
        if (arraySizes[i] != arraySizes[i-1]) {
            throw(OpenMMException("State Deserialization Particle Size Mismatch, check number of particles in Forces, Velocities, Positions!"));
        }
    }
    builder.setPeriodicBoxVectors(outAVec, outBVec, outCVec);
    State *s = new State();
    *s = builder.getState();
    return s;
}
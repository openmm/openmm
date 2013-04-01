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

#include "openmm/serialization/StateProxy.h"
#include <OpenMM.h>
#include <map>

using namespace std;
using namespace OpenMM;

StateProxy::StateProxy() : SerializationProxy("State") {

}

void StateProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const State& s = *reinterpret_cast<const State*>(object);
    node.setDoubleProperty("time", s.getTime());
    Vec3 a,b,c;
    s.getPeriodicBoxVectors(a,b,c);
    SerializationNode& boxVectorsNode = node.createChildNode("PeriodicBoxVectors");
    boxVectorsNode.createChildNode("A").setDoubleProperty("x", a[0]).setDoubleProperty("y", a[1]).setDoubleProperty("z", a[2]);
    boxVectorsNode.createChildNode("B").setDoubleProperty("x", b[0]).setDoubleProperty("y", b[1]).setDoubleProperty("z", b[2]);
    boxVectorsNode.createChildNode("C").setDoubleProperty("x", c[0]).setDoubleProperty("y", c[1]).setDoubleProperty("z", c[2]);
    try {
        s.getParameters();
        SerializationNode& parametersNode = node.createChildNode("Parameters");
        map<string, double> stateParams = s.getParameters();
        map<string, double>::const_iterator it;
        for (it = stateParams.begin(); it!=stateParams.end();it++) {
            parametersNode.setDoubleProperty(it->first, it->second);
        }
    } catch (const OpenMMException &) {
        // do nothing
    }
    try {
        s.getPotentialEnergy();
        SerializationNode& energiesNode = node.createChildNode("Energies");
        energiesNode.setDoubleProperty("PotentialEnergy", s.getPotentialEnergy());
        energiesNode.setDoubleProperty("KineticEnergy", s.getKineticEnergy());
    } catch (const OpenMMException &) {
        // do nothing
    }
    try {
        s.getPositions();
        SerializationNode& positionsNode = node.createChildNode("Positions");
        vector<Vec3> statePositions = s.getPositions();
        for (int i=0; i<statePositions.size();i++) {
           positionsNode.createChildNode("Position").setDoubleProperty("x", statePositions[i][0]).setDoubleProperty("y", statePositions[i][1]).setDoubleProperty("z", statePositions[i][2]);
        }
    } catch (const OpenMMException &) {
        // do nothing
    }
    try {
        s.getVelocities();
        SerializationNode& velocitiesNode = node.createChildNode("Velocities");
        vector<Vec3> stateVelocities = s.getVelocities();
        for (int i=0; i<stateVelocities.size();i++) {
           velocitiesNode.createChildNode("Velocity").setDoubleProperty("x", stateVelocities[i][0]).setDoubleProperty("y", stateVelocities[i][1]).setDoubleProperty("z", stateVelocities[i][2]);
        }
    } catch (const OpenMMException &) {
        // do nothing
    }
    try {
        s.getForces();
        SerializationNode& forcesNode = node.createChildNode("Forces");
        vector<Vec3> stateForces = s.getForces();
        for (int i=0; i<stateForces.size();i++) {
            forcesNode.createChildNode("Force").setDoubleProperty("x", stateForces[i][0]).setDoubleProperty("y", stateForces[i][1]).setDoubleProperty("z", stateForces[i][2]);
        }
    } catch (const OpenMMException &) {
        // do nothing
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
    map<string, double> outStateParams;
    try {
        const SerializationNode& parametersNode = node.getChildNode("Parameters");
        // inStateParams is really a <string,double> pair, where string is the name and double is the value
        // but we want to avoid casting a string to a double and instead use the built in routines,
        map<string, string> inStateParams = parametersNode.getProperties();
        for (map<string, string>::const_iterator pit = inStateParams.begin(); pit != inStateParams.end(); pit++) {
            outStateParams[pit->first] = parametersNode.getDoubleProperty(pit->first);
        }
        types = types | State::Parameters;
    } catch (const OpenMMException &) {
        // do nothing
    }
    double potentialEnergy;
    double kineticEnergy;
    try {
        const SerializationNode& energiesNode = node.getChildNode("Energies");
        potentialEnergy = energiesNode.getDoubleProperty("PotentialEnergy");
        kineticEnergy = energiesNode.getDoubleProperty("KineticEnergy");
        types = types | State::Energy;
    } catch (const OpenMMException &) {
        // do nothing    
    }
    vector<Vec3> outPositions;
    vector<Vec3> outVelocities;
    vector<Vec3> outForces;
    try {
        const SerializationNode& positionsNode = node.getChildNode("Positions");
        for (int i = 0; i < (int) positionsNode.getChildren().size(); i++) {
            const SerializationNode& particle = positionsNode.getChildren()[i];
            outPositions.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
        }
        types = types | State::Positions;
    } catch (const OpenMMException &) {
        // do nothing    
    }
    try {
        const SerializationNode& velocitiesNode = node.getChildNode("Velocities");
        for (int i = 0; i < (int) velocitiesNode.getChildren().size(); i++) {
            const SerializationNode& particle = velocitiesNode.getChildren()[i];
            outVelocities.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
        }
        types = types | State::Velocities;
    } catch (const OpenMMException &) {
        // do nothing    
    }
    try {
        const SerializationNode& forcesNode = node.getChildNode("Forces");
        for (int i = 0; i < (int) forcesNode.getChildren().size(); i++) {
            const SerializationNode& particle = forcesNode.getChildren()[i];
            outForces.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
        }
        types = types | State::Forces;
    } catch (const OpenMMException &) {
        // do nothing
    }
    vector<int> arraySizes;
    State::StateBuilder builder(outTime);
    if (types & State::Positions) {
        builder.setPositions(outPositions);
        arraySizes.push_back(outPositions.size());
    }  
    if (types & State::Velocities) {
        builder.setVelocities(outVelocities);
        arraySizes.push_back(outVelocities.size());
    }
    if (types & State::Forces) {
        builder.setForces(outForces);
        arraySizes.push_back(outForces.size());
    }
    if (types & State::Energy) {
        builder.setEnergy(kineticEnergy, potentialEnergy);
    }
    if (types & State::Parameters) {
        builder.setParameters(outStateParams);
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
/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2017 Stanford University and the Authors.      *
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

#include "openmm/serialization/SystemProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

SystemProxy::SystemProxy() : SerializationProxy("System") {
}

void SystemProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    node.setStringProperty("openmmVersion", Platform::getOpenMMVersion());
    const System& system = *reinterpret_cast<const System*>(object);
    Vec3 a, b, c;
    system.getDefaultPeriodicBoxVectors(a, b, c);
    SerializationNode& box = node.createChildNode("PeriodicBoxVectors");
    box.createChildNode("A").setDoubleProperty("x", a[0]).setDoubleProperty("y", a[1]).setDoubleProperty("z", a[2]);
    box.createChildNode("B").setDoubleProperty("x", b[0]).setDoubleProperty("y", b[1]).setDoubleProperty("z", b[2]);
    box.createChildNode("C").setDoubleProperty("x", c[0]).setDoubleProperty("y", c[1]).setDoubleProperty("z", c[2]);
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < system.getNumParticles(); i++) {
        SerializationNode& particle = particles.createChildNode("Particle").setDoubleProperty("mass", system.getParticleMass(i));
        if (system.isVirtualSite(i)) {
            const VirtualSite& vsite = system.getVirtualSite(i);
            if (typeid(vsite) == typeid(TwoParticleAverageSite)) {
                const TwoParticleAverageSite& site = dynamic_cast<const TwoParticleAverageSite&>(vsite);
                particle.createChildNode("TwoParticleAverageSite").setIntProperty("p1", site.getParticle(0)).setIntProperty("p2", site.getParticle(1)).setDoubleProperty("w1", site.getWeight(0)).setDoubleProperty("w2", site.getWeight(1));
            }
            else if (typeid(vsite) == typeid(ThreeParticleAverageSite)) {
                const ThreeParticleAverageSite& site = dynamic_cast<const ThreeParticleAverageSite&>(vsite);
                particle.createChildNode("ThreeParticleAverageSite").setIntProperty("p1", site.getParticle(0)).setIntProperty("p2", site.getParticle(1)).setIntProperty("p3", site.getParticle(2)).setDoubleProperty("w1", site.getWeight(0)).setDoubleProperty("w2", site.getWeight(1)).setDoubleProperty("w3", site.getWeight(2));
            }
            else if (typeid(vsite) == typeid(OutOfPlaneSite)) {
                const OutOfPlaneSite& site = dynamic_cast<const OutOfPlaneSite&>(vsite);
                particle.createChildNode("OutOfPlaneSite").setIntProperty("p1", site.getParticle(0)).setIntProperty("p2", site.getParticle(1)).setIntProperty("p3", site.getParticle(2)).setDoubleProperty("w12", site.getWeight12()).setDoubleProperty("w13", site.getWeight13()).setDoubleProperty("wc", site.getWeightCross());
            }
            else if (typeid(vsite) == typeid(LocalCoordinatesSite)) {
                const LocalCoordinatesSite& site = dynamic_cast<const LocalCoordinatesSite&>(vsite);
                int numParticles = site.getNumParticles();
                vector<double> wo, wx, wy;
                site.getOriginWeights(wo);
                site.getXWeights(wx);
                site.getYWeights(wy);
                Vec3 p = site.getLocalPosition();
                SerializationNode& siteNode = particle.createChildNode("LocalCoordinatesSite");
                siteNode.setDoubleProperty("pos1", p[0]).setDoubleProperty("pos2", p[1]).setDoubleProperty("pos3", p[2]);
                for (int j = 0; j < numParticles; j++) {
                    stringstream ss;
                    ss << (j+1);
                    string index = ss.str();
                    siteNode.setIntProperty("p"+index, site.getParticle(j));
                    siteNode.setDoubleProperty("wo"+index, wo[j]);
                    siteNode.setDoubleProperty("wx"+index, wx[j]);
                    siteNode.setDoubleProperty("wy"+index, wy[j]);
                }
            }
        }
    }
    SerializationNode& constraints = node.createChildNode("Constraints");
    for (int i = 0; i < system.getNumConstraints(); i++) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(i, particle1, particle2, distance);
        constraints.createChildNode("Constraint").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setDoubleProperty("d", distance);
    }
    SerializationNode& forces = node.createChildNode("Forces");
    for (int i = 0; i < system.getNumForces(); i++)
        forces.createChildNode("Force", &system.getForce(i));
}

void* SystemProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    System* system = new System();
    try {
        const SerializationNode& box = node.getChildNode("PeriodicBoxVectors");
        const SerializationNode& boxa = box.getChildNode("A");
        const SerializationNode& boxb = box.getChildNode("B");
        const SerializationNode& boxc = box.getChildNode("C");
        Vec3 a(boxa.getDoubleProperty("x"), boxa.getDoubleProperty("y"), boxa.getDoubleProperty("z"));
        Vec3 b(boxb.getDoubleProperty("x"), boxb.getDoubleProperty("y"), boxb.getDoubleProperty("z"));
        Vec3 c(boxc.getDoubleProperty("x"), boxc.getDoubleProperty("y"), boxc.getDoubleProperty("z"));
        system->setDefaultPeriodicBoxVectors(a, b, c);
        const SerializationNode& particles = node.getChildNode("Particles");
        for (int i = 0; i < (int) particles.getChildren().size(); i++) {
            system->addParticle(particles.getChildren()[i].getDoubleProperty("mass"));
            if (particles.getChildren()[i].getChildren().size() > 0) {
                const SerializationNode& vsite = particles.getChildren()[i].getChildren()[0];
                if (vsite.getName() == "TwoParticleAverageSite")
                    system->setVirtualSite(i, new TwoParticleAverageSite(vsite.getIntProperty("p1"), vsite.getIntProperty("p2"), vsite.getDoubleProperty("w1"), vsite.getDoubleProperty("w2")));
                else if (vsite.getName() == "ThreeParticleAverageSite")
                    system->setVirtualSite(i, new ThreeParticleAverageSite(vsite.getIntProperty("p1"), vsite.getIntProperty("p2"), vsite.getIntProperty("p3"), vsite.getDoubleProperty("w1"), vsite.getDoubleProperty("w2"), vsite.getDoubleProperty("w3")));
                else if (vsite.getName() == "OutOfPlaneSite")
                    system->setVirtualSite(i, new OutOfPlaneSite(vsite.getIntProperty("p1"), vsite.getIntProperty("p2"), vsite.getIntProperty("p3"), vsite.getDoubleProperty("w12"), vsite.getDoubleProperty("w13"), vsite.getDoubleProperty("wc")));
                else if (vsite.getName() == "LocalCoordinatesSite") {
                    vector<int> particles;
                    vector<double> wo, wx, wy;
                    for (int j = 0; ; j++) {
                        stringstream ss;
                        ss << (j+1);
                        string index = ss.str();
                        if (!vsite.hasProperty("p"+index))
                            break;
                        particles.push_back(vsite.getIntProperty("p"+index));
                        wo.push_back(vsite.getDoubleProperty("wo"+index));
                        wx.push_back(vsite.getDoubleProperty("wx"+index));
                        wy.push_back(vsite.getDoubleProperty("wy"+index));
                    }
                    Vec3 p(vsite.getDoubleProperty("pos1"), vsite.getDoubleProperty("pos2"), vsite.getDoubleProperty("pos3"));
                    system->setVirtualSite(i, new LocalCoordinatesSite(particles, wo, wx, wy, p));
                }
            }
        }
        const SerializationNode& constraints = node.getChildNode("Constraints");
        for (auto& constraint : constraints.getChildren())
            system->addConstraint(constraint.getIntProperty("p1"), constraint.getIntProperty("p2"), constraint.getDoubleProperty("d"));
        const SerializationNode& forces = node.getChildNode("Forces");
        for (auto& force : forces.getChildren())
            system->addForce(force.decodeObject<Force>());
    }
    catch (...) {
        delete system;
        throw;
    }
    return system;
}
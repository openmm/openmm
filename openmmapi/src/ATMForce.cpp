
#include "openmm/ATMForce.h"
#include "openmm/Force.h"
#include "openmm/serialization/XmlSerializer.h"
#include "openmm/internal/ATMForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

int ATMForce::addParticle(int particle, double dx, double dy, double dz) {
    particles.push_back(ParticleInfo(particle, dx, dy, dz));
    return particles.size()-1;
}

void ATMForce::getParticleParameters(int index, int& particle, double& dx, double &dy, double &dz) const {
    ASSERT_VALID_INDEX(index, particles);
    particle = particles[index].particle;
    dx = particles[index].dx;
    dy = particles[index].dy;
    dz = particles[index].dz;

}

void ATMForce::setParticleParameters(int index, int particle, double dx, double dy, double dz){
    ASSERT_VALID_INDEX(index, particles);
    particles[index].particle = particle;
    particles[index].dx = dx;
    particles[index].dy = dy;
    particles[index].dz = dz;
}

/* TO DO
int ATMForce::addForce(Force* force){
  Force* newforce = XmlSerializer::clone<Force>(*force);
  forces.push_back(newforce);
  return forces.size()-1;
}
*/

ForceImpl* ATMForce::createImpl() const {
    return new ATMForceImpl(*this);
}

void ATMForce::updateParametersInContext(OpenMM::Context& context) {
  dynamic_cast<ATMForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

double ATMForce::getPerturbationEnergy(const OpenMM::Context& context)  const {
  return dynamic_cast<const ATMForceImpl&>(getImplInContext(context)).getPerturbationEnergy(); 
}


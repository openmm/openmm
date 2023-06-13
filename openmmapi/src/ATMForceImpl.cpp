
#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/internal/ATMForceImpl.h"

#include "openmm/NonbondedForce.h"
#include "openmm/kernels.h"
#include "openmm/serialization/XmlSerializer.h"
#include "openmm/Vec3.h"

#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include <cmath>
#include <map>
#include <set>
#include <sstream>

#include <iostream>

using namespace OpenMM;
using namespace std;

ATMForceImpl::ATMForceImpl(const ATMForce& owner) : owner(owner), innerIntegrator1(1.0), innerIntegrator2(1.0), hasInitializedInnerContexts(false){
}

ATMForceImpl::~ATMForceImpl() {
}

void ATMForceImpl::copysystem(const OpenMM::System& system, OpenMM::System& innerSystem){
  
  //copy particles
  for (int i = 0; i < system.getNumParticles(); i++) 
    innerSystem.addParticle(system.getParticleMass(i));

  //copy constraints
  for (int i = 0; i < system.getNumConstraints(); i++) {
    int particle1, particle2;
    double distance;
    system.getConstraintParameters(i, particle1, particle2, distance);
    innerSystem.addConstraint(particle1, particle2, distance);
  }

  //copy periodic box dimensions
  Vec3 a, b, c;
  system.getDefaultPeriodicBoxVectors(a, b, c);
  innerSystem.setDefaultPeriodicBoxVectors(a, b, c);

  //add system forces other than those belonging to the ATMMetaForce group to the inner contexts
  int atmforcegroup = owner.getForceGroup();
  int numforces = system.getNumForces();
  for (int i=0; i<numforces; i++){
    const Force &force = system.getForce(i);
    int group = force.getForceGroup();
    if (group != atmforcegroup){
      Force* newforce = XmlSerializer::clone<Force>(force);
      newforce->setForceGroup(group);
      NonbondedForce* nonbonded = dynamic_cast<NonbondedForce*>(newforce);
      if (nonbonded != NULL)
	nonbonded->setReciprocalSpaceForceGroup(-1);
      innerSystem.addForce(newforce);
    }
  }

}
  
void ATMForceImpl::initialize(ContextImpl& context) {
  const OpenMM::System& system = context.getSystem();

  int atmforcegroup = owner.getForceGroup();

  // only forces in designated force groups are evaluated in the inner contexts
  variable_force_groups_mask = 0;
  vector<int> varforcegroups = owner.getVariableForceGroups();
  for (int i=0; i<varforcegroups.size() ; i++){
    if(varforcegroups[i] == atmforcegroup)
      throw OpenMMException("The ATM Meta Force group cannot be one of the variable force groups.");
    variable_force_groups_mask += 1<<varforcegroups[i];
  }
  
  copysystem(system, innerSystem1);
  copysystem(system, innerSystem2);

  kernel = context.getPlatform().createKernel(CalcATMForceKernel::Name(), context);
  kernel.getAs<CalcATMForceKernel>().initialize(context.getSystem(), owner);
}

double ATMForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
  if(! hasInitializedInnerContexts) {
    hasInitializedInnerContexts = true;
    
    innerContext1 = context.createLinkedContext(innerSystem1, innerIntegrator1);
    innerContext2 = context.createLinkedContext(innerSystem2, innerIntegrator2);

    vector<Vec3> pos;
    context.getPositions(pos);
    innerContext1->setPositions(pos);
    innerContext2->setPositions(pos);
  }

  if ((groups&(1<<owner.getForceGroup())) == 0) return 0.0;
  
  bool do_energy = true;
  ContextImpl& innercontextimpl1 = getContextImpl(*innerContext1);
  ContextImpl& innercontextimpl2 = getContextImpl(*innerContext2);

  //copies the coordinates etc. from the context to the inner contexts
  kernel.getAs<CalcATMForceKernel>().copyState(context, innercontextimpl1, innercontextimpl2);

  //evaluate variable energy and forces for original system
  double State1Energy = innercontextimpl1.calcForcesAndEnergy(true, do_energy, variable_force_groups_mask);
  
  //evaluate variable energy and force for the displaced system
  double State2Energy = innercontextimpl2.calcForcesAndEnergy(true, do_energy, variable_force_groups_mask);
  
  //evaluate the alchemical energy
  double energy = kernel.getAs<CalcATMForceKernel>().execute(context,
  								 innercontextimpl1, innercontextimpl2,
  								 State1Energy, State2Energy,
  								 includeForces, includeEnergy);

  //retrieve the perturbation energy
  PerturbationEnergy = kernel.getAs<CalcATMForceKernel>().getPerturbationEnergy();

  return (includeEnergy ? energy : 0.0);
}

std::map<std::string, double> ATMForceImpl::getDefaultParameters(){
  std::map<std::string, double> parameters;
  parameters[ATMForce::Lambda1()] = getOwner().getDefaultLambda1();
  parameters[ATMForce::Lambda2()] = getOwner().getDefaultLambda2();
  parameters[ATMForce::Alpha()]   = getOwner().getDefaultAlpha();
  parameters[ATMForce::U0()]      = getOwner().getDefaultU0();
  parameters[ATMForce::W0()]      = getOwner().getDefaultW0();
  parameters[ATMForce::Umax()]    = getOwner().getDefaultUmax();
  parameters[ATMForce::Ubcore()]  = getOwner().getDefaultUbcore();
  parameters[ATMForce::Acore()]   = getOwner().getDefaultAcore();
  parameters[ATMForce::Direction()]   = getOwner().getDefaultDirection();
  return parameters;
}


std::vector<std::string> ATMForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcATMForceKernel::Name());
    return names;
}

void ATMForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcATMForceKernel>().copyParametersToContext(context, owner);
}

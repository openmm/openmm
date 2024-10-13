//
// Created by babaid on 05.10.24.
//

#include "openmm/internal/ExternalPuremdForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <sstream>
#include <algorithm>
#include "openmm/PuremdInterface.h"


using namespace OpenMM;
using namespace std;
ExternalPuremdForceImpl::ExternalPuremdForceImpl(const ExternalPuremdForce &owner): CustomCPPForceImpl(owner), owner(owner)
{
    for(int i = 0; i<owner.getNumAtoms(); ++i)
    {
      int particle;
      char* symbol;
      int isqm;
      owner.getParticleParameters(i, particle, symbol, isqm);
      if(isqm)
      {
        qmParticles.emplace_back(particle);
        qmSymbols.emplace_back(symbol[0]);
        qmSymbols.emplace_back(symbol[1]);
      }
      else
      {
        mmParticles.emplace_back(particle);
        mmSymbols.emplace_back(symbol[0]);
        mmSymbols.emplace_back(symbol[1]);
      }
    }
}
void ExternalPuremdForceImpl::getBoxInfo(ContextImpl& context, std::vector<double>& simBoxInfo)
{
  std::vector<Vec3> PeriodicBoxVectors(3);
  context.getPeriodicBoxVectors(PeriodicBoxVectors[0], PeriodicBoxVectors[1], PeriodicBoxVectors[2]);
  for (int i=0; i<3; i++)
  {
    //AA
    simBoxInfo[i] = std::sqrt(PeriodicBoxVectors[i].dot(PeriodicBoxVectors[i]))*10;
  }
  simBoxInfo[3] = std::acos(PeriodicBoxVectors[1].dot(PeriodicBoxVectors[2])/(simBoxInfo[1]*simBoxInfo[2]))*180.0/M_PI;
  simBoxInfo[4] =  std::acos(PeriodicBoxVectors[0].dot(PeriodicBoxVectors[2])/(simBoxInfo[0]*simBoxInfo[2]))*180.0/M_PI;
  simBoxInfo[5] =  std::acos(PeriodicBoxVectors[0].dot(PeriodicBoxVectors[1])/(simBoxInfo[0]*simBoxInfo[1]))*180.0/M_PI;
}


double ExternalPuremdForceImpl::computeForce(ContextImpl& context, const std::vector<Vec3> &positions, std::vector<Vec3>& forces)
{
  // need to seperate positions
  std::vector<Vec3> transoformedPositions;
  //next we need to seperate and flatten the QM/MM positions and convert to AA
  int numQm = qmParticles.size(), numMm = mmParticles.size();
  std::vector<double> qmPos, mmPos_q;

  std::for_each(qmParticles.begin(), qmParticles.end(), [&](int Index){
    qmPos.emplace_back(positions[Index][0]*10);
    qmPos.emplace_back(positions[Index][1]*10);
    qmPos.emplace_back(positions[Index][2]*10);
  });



  //retrieve charges from the context. Had to introduce some changes to classes Context, ContextImpl,
  // UpdateStateDataKernel, CommonUpdateStateDataKernel
  std::vector<double> charges;
  context.getCharges(charges);

  std::for_each(mmParticles.begin(), mmParticles.end(), [&](int Index){
    mmPos_q.emplace_back(positions[Index][0]*10);
    mmPos_q.emplace_back(positions[Index][1]*10);
    mmPos_q.emplace_back(positions[Index][2]*10);
    mmPos_q.emplace_back(charges[Index]);
  });


  //get the box size. move this into a function
  std::vector<double> simBoxInfo(6);
  getBoxInfo(context, simBoxInfo);

  std::vector<double> qmForces, mmForces;
  std::vector<double> qmQ;
  double energy;
  Interface.getReaxffPuremdForces(numQm, qmSymbols, qmPos,
                                  numMm, mmSymbols, mmPos_q,
                                  simBoxInfo, qmForces, mmForces, qmQ,
                                  energy);

  // merge the qm and mm forces, additionally transform the scale
  std::vector<Vec3> transformedForces(owner.getNumAtoms());
  int Index;

  for (size_t i=0; i<qmParticles.size(); ++i)
  {
      Index = qmParticles[i];

      transformedForces[Index][0] = qmForces[i*3];
      transformedForces[Index][1] = qmForces[i*3 + 1];
      transformedForces[Index][2] = qmForces[i*3 + 2];
      charges[Index] = qmQ[i];
  }

  for (size_t i=0; i<mmParticles.size(); ++i)
  {
      Index = mmParticles[i];
      transformedForces[Index][0] = mmForces[i*3];
      transformedForces[Index][1] = mmForces[i*3 + 1];
      transformedForces[Index][2] = mmForces[i*3 + 2];
  }

  //update charges
  context.setCharges(charges);
  //copy forces and transform from Angstroms * Daltons / ps^2 to kJ/mol/nm
  double conversionFactor = 2.76E-7;
  for(size_t i =0;i<forces.size();++i) forces[i] = transformedForces[i]*conversionFactor;
  //done
  return energy;
}
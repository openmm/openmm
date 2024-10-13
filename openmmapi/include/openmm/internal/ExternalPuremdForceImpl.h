//
// Created by babaid on 05.10.24.
//

#ifndef OPENMM_EXTERNALPUREMDFORCEIMPL_H
#define OPENMM_EXTERNALPUREMDFORCEIMPL_H
#include "CustomCPPForceImpl.h"
#include "PuremdInterface.h"
#include "openmm/ExternalPuremdForce.h"
#include "openmm/Kernel.h"
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace OpenMM {

/**
 * This is the internal implementation of ExternalPuremdForce.
 */
    class ExternalPuremdForceImpl : public CustomCPPForceImpl {
    public:
      ExternalPuremdForceImpl(const ExternalPuremdForce &owner);
      ~ExternalPuremdForceImpl() = default;
      double computeForce(ContextImpl& context, const std::vector<Vec3> &positions, std::vector<Vec3>& forces) override;
      const ExternalPuremdForce& getOwner() const{
        return owner;
      }
    private:
      std::vector<char> qmSymbols;
      std::vector<char> mmSymbols;
      std::vector<int> qmParticles;
      std::vector<int> mmParticles;
      /**@private
       *
       * @param context
       * @param simBoxInfo
       */
      static void getBoxInfo(ContextImpl& context, std::vector<double>& simBoxInfo);
      const ExternalPuremdForce & owner;
      PuremdInterface Interface;
    };

} // namespace OpenMM
#endif //OPENMM_EXTERNALPUREMDFORCE_H

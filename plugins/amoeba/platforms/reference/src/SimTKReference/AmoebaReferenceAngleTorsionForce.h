#ifndef __AmoebaReferenceAngleTorsionForce_H__
#define __AmoebaReferenceAngleTorsionForce_H__

#include "RealVec.h"
#include <vector>

namespace OpenMM {

class AmoebaReferenceAngleTorsionForce {
public:
	AmoebaReferenceAngleTorsionForce() {};
	~AmoebaReferenceAngleTorsionForce() {};
	RealOpenMM calculateForceAndEnergy(int numAngleTorsions, std::vector<OpenMM::RealVec>& posData,
			const std::vector<int>& particle1,
			const std::vector<int>&  particle2,
			const std::vector<int>&  particle3,
			const std::vector<int>&  particle4,
			const std::vector<RealOpenMM>& angleCBAParameters,
			const std::vector<RealOpenMM>& angleDCBParameters,
			const std::vector<RealOpenMM>&  k1Parameters,
			const std::vector<RealOpenMM>&  k2Parameters,
			const std::vector<RealOpenMM>&  k3Parameters,
			const std::vector<RealOpenMM>&  k4Parameters,
			const std::vector<RealOpenMM>&  k5Parameters,
			const std::vector<RealOpenMM>&  k6Parameters,
			std::vector<OpenMM::RealVec>& forceData) const;
private:
	RealOpenMM calculateAngleTorsionIxn(const OpenMM::RealVec& positionAtomA, const OpenMM::RealVec& positionAtomB,
			const OpenMM::RealVec& positionAtomC, const OpenMM::RealVec& positionAtomD,
			RealOpenMM angleCBA,      RealOpenMM angleDCB,
			RealOpenMM k1Parameter, RealOpenMM k2Parameter, RealOpenMM k3Parameter,
			RealOpenMM k4Parameter, RealOpenMM k5Parameter, RealOpenMM k6Parameter,
			OpenMM::RealVec* forces) const;
	void getAddition3(std::vector<RealOpenMM>& xVector, std::vector<RealOpenMM>& yVector, std::vector<RealOpenMM>& zVector) const;
	void getAddition3(std::vector<RealOpenMM>& xVector, std::vector<RealOpenMM>& yVector) const;
	void getSubstraction3(std::vector<RealOpenMM>& xVector, std::vector<RealOpenMM>& yVector, std::vector<RealOpenMM>& zVector) const;
	void getSubstraction3(std::vector<RealOpenMM>& xVector, std::vector<RealOpenMM>& yVector) const;
	void getScale3(std::vector<RealOpenMM>& xVector, RealOpenMM scale, std::vector<RealOpenMM>& yVector) const;
	void getScale3(std::vector<RealOpenMM>& xVector, RealOpenMM scale) const;

};

} // namespace OpanMM

#endif // __AmoebaReferenceAngleTorsionForce_H__

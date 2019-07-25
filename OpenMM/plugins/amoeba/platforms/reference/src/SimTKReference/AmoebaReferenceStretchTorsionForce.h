#ifndef __AmoebaReferenceStretchTorsionForce_H__
#define __AmoebaReferenceStretchTorsionForce_H__

#include "RealVec.h"
#include <vector>

namespace OpenMM {

class AmoebaReferenceStretchTorsionForce {
public:
	AmoebaReferenceStretchTorsionForce() {};
	~AmoebaReferenceStretchTorsionForce() {};
	RealOpenMM calculateForceAndEnergy(int numStrectchTorsions, std::vector<OpenMM::RealVec>& posData,
			const std::vector<int>& particle1,
			const std::vector<int>&  particle2,
			const std::vector<int>&  particle3,
			const std::vector<int>&  particle4,
			const std::vector<RealOpenMM>& lengthBAParameters,
			const std::vector<RealOpenMM>& lengthCBParameters,
			const std::vector<RealOpenMM>& lengthDCParameters,
			const std::vector<RealOpenMM>&  k1Parameters,
			const std::vector<RealOpenMM>&  k2Parameters,
			const std::vector<RealOpenMM>&  k3Parameters,
			const std::vector<RealOpenMM>&  k4Parameters,
			const std::vector<RealOpenMM>&  k5Parameters,
			const std::vector<RealOpenMM>&  k6Parameters,
			const std::vector<RealOpenMM>&  k7Parameters,
			const std::vector<RealOpenMM>&  k8Parameters,
			const std::vector<RealOpenMM>&  k9Parameters,
			std::vector<OpenMM::RealVec>& forceData) const;
private:
	RealOpenMM calculateStretchTorsionIxn(const OpenMM::RealVec& positionAtomA, const OpenMM::RealVec& positionAtomB,
					const OpenMM::RealVec& positionAtomC, const OpenMM::RealVec& positionAtomD,
					RealOpenMM lengthBA,      RealOpenMM lengthCB,  RealOpenMM lengthDC,
					RealOpenMM k1Parameter, RealOpenMM k2Parameter, RealOpenMM k3Parameter,
					RealOpenMM k4Parameter, RealOpenMM k5Parameter, RealOpenMM k6Parameter,
					RealOpenMM k7Parameter, RealOpenMM k8Parameter, RealOpenMM k9Parameter,
					OpenMM::RealVec* forces) const;
	void getAddition3(std::vector<RealOpenMM>& xVector, std::vector<RealOpenMM>& yVector, std::vector<RealOpenMM>& zVector) const;
	void getAddition3(std::vector<RealOpenMM>& xVector, std::vector<RealOpenMM>& yVector) const;
	void getSubstraction3(std::vector<RealOpenMM>& xVector, std::vector<RealOpenMM>& yVector, std::vector<RealOpenMM>& zVector) const;
	void getSubstraction3(std::vector<RealOpenMM>& xVector, std::vector<RealOpenMM>& yVector) const;
	void getScale3(std::vector<RealOpenMM>& xVector, RealOpenMM scale, std::vector<RealOpenMM>& yVector) const;
	void getScale3(std::vector<RealOpenMM>& xVector, RealOpenMM scale) const;

};

}

#endif // _AmoebaReferenceStretchTorsionForce___

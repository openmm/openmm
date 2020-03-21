#ifndef __AmoebaReferenceAngleTorsionForce_H__
#define __AmoebaReferenceAngleTorsionForce_H__

#include "openmm/Vec3.h"
#include <vector>

namespace OpenMM {

class AmoebaReferenceAngleTorsionForce {
public:
	AmoebaReferenceAngleTorsionForce() : usePeriodic(false) {};
	~AmoebaReferenceAngleTorsionForce() {};

        void setPeriodic(OpenMM::Vec3* vectors);

	double calculateForceAndEnergy(int numAngleTorsions, std::vector<Vec3>& posData,
			const std::vector<int>& particle1,
			const std::vector<int>&  particle2,
			const std::vector<int>&  particle3,
			const std::vector<int>&  particle4,
			const std::vector<double>& angleCBAParameters,
			const std::vector<double>& angleDCBParameters,
			const std::vector<double>&  k1Parameters,
			const std::vector<double>&  k2Parameters,
			const std::vector<double>&  k3Parameters,
			const std::vector<double>&  k4Parameters,
			const std::vector<double>&  k5Parameters,
			const std::vector<double>&  k6Parameters,
			std::vector<OpenMM::Vec3>& forceData) const;
private:
        bool usePeriodic;
        Vec3 boxVectors[3];

	double calculateAngleTorsionIxn(const Vec3& positionAtomA, const Vec3& positionAtomB,
			const Vec3& positionAtomC, const Vec3& positionAtomD,
			double angleCBA,      double angleDCB,
			double k1Parameter, double k2Parameter, double k3Parameter,
			double k4Parameter, double k5Parameter, double k6Parameter,
			Vec3* forces) const;
	void getAddition3(std::vector<double>& xVector, std::vector<double>& yVector, std::vector<double>& zVector) const;
	void getAddition3(std::vector<double>& xVector, std::vector<double>& yVector) const;
	void getSubstraction3(std::vector<double>& xVector, std::vector<double>& yVector, std::vector<double>& zVector) const;
	void getSubstraction3(std::vector<double>& xVector, std::vector<double>& yVector) const;
	void getScale3(std::vector<double>& xVector, double scale, std::vector<double>& yVector) const;
	void getScale3(std::vector<double>& xVector, double scale) const;

};

} // namespace OpanMM

#endif // __AmoebaReferenceAngleTorsionForce_H__

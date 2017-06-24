#include "AmoebaReferenceForce.h"
#include "AmoebaReferenceAngleTorsionForce.h"
#include <vector>

using std::vector;
using namespace OpenMM;

void AmoebaReferenceAngleTorsionForce::getAddition3(std::vector<RealOpenMM>& xVector, std::vector<RealOpenMM>& yVector, std::vector<RealOpenMM>& zVector) const {
	zVector[0] = xVector[0] + yVector[0];
	zVector[1] = xVector[1] + yVector[1];
	zVector[2] = xVector[2] + yVector[2];
}

void AmoebaReferenceAngleTorsionForce::getAddition3(std::vector<RealOpenMM>& xVector, std::vector<RealOpenMM>& yVector) const {
	xVector[0] += yVector[0];
	xVector[1] += yVector[1];
	xVector[2] += yVector[2];
}

void AmoebaReferenceAngleTorsionForce::getSubstraction3(std::vector<RealOpenMM>& xVector, std::vector<RealOpenMM>& yVector, std::vector<RealOpenMM>& zVector) const {
	zVector[0] = xVector[0] - yVector[0];
	zVector[1] = xVector[1] - yVector[1];
	zVector[2] = xVector[2] - yVector[2];
}

void AmoebaReferenceAngleTorsionForce::getSubstraction3(std::vector<RealOpenMM>& xVector, std::vector<RealOpenMM>& yVector) const {
	xVector[0] -= yVector[0];
	xVector[1] -= yVector[1];
	xVector[2] -= yVector[2];
}

void AmoebaReferenceAngleTorsionForce::getScale3(std::vector<RealOpenMM>& xVector, RealOpenMM scale, std::vector<RealOpenMM>& yVector) const {
	yVector[0] = xVector[0] * scale;
	yVector[1] = xVector[1] * scale;
	yVector[2] = xVector[2] * scale;
}

void AmoebaReferenceAngleTorsionForce::getScale3(std::vector<RealOpenMM>& xVector, RealOpenMM scale) const {
	xVector[0] *= scale;
	xVector[1] *= scale;
	xVector[2] *= scale;
}

RealOpenMM AmoebaReferenceAngleTorsionForce::calculateForceAndEnergy(int numAngleTorsions, std::vector<OpenMM::RealVec>& posData,
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
		std::vector<OpenMM::RealVec>& forceData) const {
	RealOpenMM energy      = 0.0;
	for (unsigned int ii = 0; ii < static_cast<unsigned int>(numAngleTorsions); ii++) {
		int particle1Index      = particle1[ii];
		int particle2Index      = particle2[ii];
		int particle3Index      = particle3[ii];
		int particle4Index      = particle4[ii];
		RealOpenMM cbaAngle     = angleCBAParameters[ii];
		RealOpenMM dcbAngle     = angleDCBParameters[ii];
		RealOpenMM K1      = k1Parameters[ii];
		RealOpenMM K2      = k2Parameters[ii];
		RealOpenMM K3      = k3Parameters[ii];
		RealOpenMM K4      = k4Parameters[ii];
		RealOpenMM K5      = k5Parameters[ii];
		RealOpenMM K6      = k6Parameters[ii];
		RealVec forces[4];
		energy                 += calculateAngleTorsionIxn(posData[particle1Index], posData[particle2Index], posData[particle3Index], posData[particle4Index],
				cbaAngle, dcbAngle, K1, K2, K3, K4, K5, K6, forces);
		for(int jj = 0; jj < 3; jj++) {
			forceData[particle1Index][jj] -= forces[0][jj];
			forceData[particle2Index][jj] -= forces[1][jj];
			forceData[particle3Index][jj] -= forces[2][jj];
			forceData[particle4Index][jj] -= forces[3][jj];
		}
	}
	return energy;
}

RealOpenMM AmoebaReferenceAngleTorsionForce::calculateAngleTorsionIxn(const OpenMM::RealVec& positionAtomA, const OpenMM::RealVec& positionAtomB,
		const OpenMM::RealVec& positionAtomC, const OpenMM::RealVec& positionAtomD,
		RealOpenMM angleCBA,      RealOpenMM angleDCB,
		RealOpenMM k1Parameter, RealOpenMM k2Parameter, RealOpenMM k3Parameter,
		RealOpenMM k4Parameter, RealOpenMM k5Parameter, RealOpenMM k6Parameter,
		OpenMM::RealVec* forces) const {

	static const RealOpenMM zero          = 0.0;
        static const RealOpenMM one           = 1.0;
        static const RealOpenMM two           = 2.0;
        static const RealOpenMM three         = 3.0;

	enum { A = 0, B, C, D, LastAtomIndex };
	enum { BA = 0, CB, DC, CA, DB, BAxCB, CBxDC, TxU, DEDT, DEDU, TxBA, TxCB, UxCB, UxDC, DEDTxCB, CAxDEDT, DEDUxDC, DEDTxBA, DBxDEDU, DEDUxCB, LastDeltaIndex };

	std::vector<RealOpenMM> deltaR[LastDeltaIndex];
	for (unsigned int ii = 0; ii < LastDeltaIndex; ii++) {
		deltaR[ii].resize(3);
	}
	AmoebaReferenceForce::loadDeltaR(positionAtomA, positionAtomB, deltaR[BA]);
	AmoebaReferenceForce::loadDeltaR(positionAtomB, positionAtomC, deltaR[CB]);
	AmoebaReferenceForce::loadDeltaR(positionAtomC, positionAtomD, deltaR[DC]);

	RealOpenMM rBA2 = AmoebaReferenceForce::getNormSquared3(deltaR[BA]);
	RealOpenMM rCB2 = AmoebaReferenceForce::getNormSquared3(deltaR[CB]);
	RealOpenMM rDC2 = AmoebaReferenceForce::getNormSquared3(deltaR[DC]);
	AmoebaReferenceForce::getCrossProduct(deltaR[BA], deltaR[CB], deltaR[BAxCB]);
	AmoebaReferenceForce::getCrossProduct(deltaR[CB], deltaR[DC], deltaR[CBxDC]);
	AmoebaReferenceForce::getCrossProduct(deltaR[BAxCB], deltaR[CBxDC], deltaR[TxU]);
	RealOpenMM  rT2 = AmoebaReferenceForce::getNormSquared3(deltaR[BAxCB]);
	RealOpenMM  rU2 = AmoebaReferenceForce::getNormSquared3(deltaR[CBxDC]);
	RealOpenMM  rTrU = SQRT(rT2 * rU2);

	if (rTrU <= zero) {
		return zero;
	}

	AmoebaReferenceForce::loadDeltaR(positionAtomA, positionAtomC, deltaR[CA]);
	AmoebaReferenceForce::loadDeltaR(positionAtomB, positionAtomD, deltaR[DB]);
	RealOpenMM  rCB = SQRT(rCB2);
	RealOpenMM  cosine = AmoebaReferenceForce::getDotProduct3(deltaR[BAxCB], deltaR[CBxDC]) / rTrU;
	RealOpenMM  sine = AmoebaReferenceForce::getDotProduct3(deltaR[CB], deltaR[TxU]) / (rCB * rTrU);

        RealOpenMM  c1 = one;
        RealOpenMM  s1 = zero;
        RealOpenMM  c2 = - one;
        RealOpenMM  s2 = zero;
        RealOpenMM  c3 = one;
        RealOpenMM  s3 = zero;
        RealOpenMM  cosine2 = cosine * cosine - sine * sine;
        RealOpenMM  sine2 = two * sine * cosine;
        RealOpenMM  cosine3 = cosine * cosine2 - sine * sine2;
        RealOpenMM  sine3 = cosine * sine2 + sine * cosine2;
        RealOpenMM  phi1 = one + (cosine * c1 + sine * s1);
        RealOpenMM  phi2 = one + (cosine2 * c2 + sine2 * s2);
        RealOpenMM  phi3 = one + (cosine3 * c3 + sine3 * s3);
        RealOpenMM  dphi1 = cosine * s1 - sine * c1;
        RealOpenMM  dphi2 = two * (cosine2 * s2 - sine2 * c2);
        RealOpenMM  dphi3 = three * (cosine3 * s3 - sine3 * c3);

        RealOpenMM  v1 = k1Parameter;
	RealOpenMM  v2 = k2Parameter;
	RealOpenMM  v3 = k3Parameter;
	RealOpenMM  dot = AmoebaReferenceForce::getDotProduct3(deltaR[BA], deltaR[CB]);
	RealOpenMM  cosang = - dot / SQRT(rBA2 * rCB2);
	RealOpenMM  angle = RADIAN * ACOS(cosang);
	RealOpenMM  dt = angle - angleCBA;
	RealOpenMM  e1 = dt * (v1 * phi1 + v2 * phi2 + v3 * phi3);
	RealOpenMM  dedphi = dt * (v1 * dphi1 + v2 * dphi2 + v3 * dphi3);
	RealOpenMM  ddt = RADIAN * (v1 * phi1 + v2 * phi2 + v3 * phi3);

	AmoebaReferenceForce::getCrossProduct(deltaR[BAxCB], deltaR[CB], deltaR[DEDT]);
	RealOpenMM tempScale = dedphi / (rT2 * rCB);
	getScale3(deltaR[DEDT], tempScale);

	AmoebaReferenceForce::getCrossProduct(deltaR[CB], deltaR[CBxDC], deltaR[DEDU]);
	tempScale = dedphi / (rU2 * rCB);
	getScale3(deltaR[DEDU], tempScale);

	std::vector<RealOpenMM> subForce[LastAtomIndex];
	for (int ii = 0; ii < LastAtomIndex; ii++) {
		subForce[ii].resize(3);
	}

	RealOpenMM  terma = - ddt / (rBA2 * SQRT(rT2));
	RealOpenMM  termc = ddt / (rCB2 * SQRT(rT2));
	AmoebaReferenceForce::getCrossProduct(deltaR[BAxCB], deltaR[BA], deltaR[TxBA]);
	getScale3(deltaR[TxBA], terma);
	AmoebaReferenceForce::getCrossProduct(deltaR[DEDT], deltaR[CB], deltaR[DEDTxCB]);
	getAddition3(deltaR[TxBA], deltaR[DEDTxCB], subForce[A]);

	AmoebaReferenceForce::getCrossProduct(deltaR[BAxCB], deltaR[CB], deltaR[TxCB]);
	getScale3(deltaR[TxCB], termc);
	getSubstraction3(deltaR[TxCB], deltaR[TxBA], subForce[B]);
	AmoebaReferenceForce::getCrossProduct(deltaR[CA], deltaR[DEDT], deltaR[CAxDEDT]);
	getAddition3(subForce[B], deltaR[CAxDEDT]);
	AmoebaReferenceForce::getCrossProduct(deltaR[DEDU], deltaR[DC], deltaR[DEDUxDC]);
	getAddition3(subForce[B], deltaR[DEDUxDC]);

	AmoebaReferenceForce::getCrossProduct(deltaR[DEDT], deltaR[BA], deltaR[DEDTxBA]);
	getSubstraction3(deltaR[DEDTxBA], deltaR[TxCB], subForce[C]);
	AmoebaReferenceForce::getCrossProduct(deltaR[DB], deltaR[DEDU], deltaR[DBxDEDU]);
	getAddition3(subForce[C], deltaR[DBxDEDU]);

	AmoebaReferenceForce::getCrossProduct(deltaR[DEDU], deltaR[CB], subForce[D]);

	v1 = k4Parameter;
	v2 = k5Parameter;
	v3 = k6Parameter;
	dot = AmoebaReferenceForce::getDotProduct3(deltaR[CB], deltaR[DC]); 
	cosang = - dot / SQRT(rCB2 * rDC2);
	angle = RADIAN * ACOS(cosang);
	dt = angle - angleDCB;
	RealOpenMM  e2 = dt * (v1 * phi1 + v2 * phi2 + v3 * phi3);
	dedphi = dt * (v1 * dphi1 + v2 * dphi2 + v3 * dphi3);
	ddt = RADIAN * (v1 * phi1 + v2 * phi2 + v3 * phi3);

	AmoebaReferenceForce::getCrossProduct(deltaR[BAxCB], deltaR[CB], deltaR[DEDT]);
	tempScale = dedphi / (rT2 * rCB);
	getScale3(deltaR[DEDT], tempScale);
	AmoebaReferenceForce::getCrossProduct(deltaR[CB], deltaR[CBxDC], deltaR[DEDU]);
	tempScale = dedphi / (rU2 * rCB);
	getScale3(deltaR[DEDU], tempScale);

	RealOpenMM  termb = - ddt / (rCB2 * SQRT(rU2));
	RealOpenMM  termd = ddt / (rDC2 * SQRT(rU2));
	AmoebaReferenceForce::getCrossProduct(deltaR[DEDT], deltaR[CB], deltaR[DEDTxCB]);
	getAddition3(subForce[A], deltaR[DEDTxCB]);

	AmoebaReferenceForce::getCrossProduct(deltaR[CBxDC], deltaR[CB], deltaR[UxCB]);
	getScale3(deltaR[UxCB], termb);
	AmoebaReferenceForce::getCrossProduct(deltaR[CA], deltaR[DEDT], deltaR[CAxDEDT]);
	AmoebaReferenceForce::getCrossProduct(deltaR[DEDU], deltaR[DC], deltaR[DEDUxDC]);
	getAddition3(subForce[B], deltaR[UxCB]);
	getAddition3(subForce[B], deltaR[CAxDEDT]);
	getAddition3(subForce[B], deltaR[DEDUxDC]);

	getSubstraction3(subForce[C],deltaR[UxCB]);
	AmoebaReferenceForce::getCrossProduct(deltaR[CBxDC], deltaR[DC], deltaR[UxDC]);
	getScale3(deltaR[UxDC], termd);
	getAddition3(subForce[C], deltaR[UxDC]);
	AmoebaReferenceForce::getCrossProduct(deltaR[DEDT], deltaR[BA], deltaR[DEDTxBA]);
	getAddition3(subForce[C], deltaR[DEDTxBA]);
	AmoebaReferenceForce::getCrossProduct(deltaR[DB], deltaR[DEDU], deltaR[DBxDEDU]);
	getAddition3(subForce[C], deltaR[DBxDEDU]);

	getSubstraction3(subForce[D], deltaR[UxDC]);
	AmoebaReferenceForce::getCrossProduct(deltaR[DEDU], deltaR[CB], deltaR[DEDUxCB]);
	getAddition3(subForce[D], deltaR[DEDUxCB]);

	RealOpenMM  e = e1 + e2;
	for (int jj = 0; jj < LastAtomIndex; jj++) {
		forces[jj][0] = subForce[jj][0];
		forces[jj][1] = subForce[jj][1];
		forces[jj][2] = subForce[jj][2];
	}
	return e;
}

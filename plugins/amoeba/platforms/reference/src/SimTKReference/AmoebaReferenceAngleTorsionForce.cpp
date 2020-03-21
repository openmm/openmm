#include "AmoebaReferenceForce.h"
#include "AmoebaReferenceAngleTorsionForce.h"
#include "SimTKOpenMMRealType.h"
#include <vector>

using std::vector;
using namespace OpenMM;

void AmoebaReferenceAngleTorsionForce::setPeriodic(OpenMM::Vec3* vectors) {
    usePeriodic = true;
    boxVectors[0] = vectors[0];
    boxVectors[1] = vectors[1];
    boxVectors[2] = vectors[2];
}

void AmoebaReferenceAngleTorsionForce::getAddition3(std::vector<double>& xVector, std::vector<double>& yVector, std::vector<double>& zVector) const {
	zVector[0] = xVector[0] + yVector[0];
	zVector[1] = xVector[1] + yVector[1];
	zVector[2] = xVector[2] + yVector[2];
}

void AmoebaReferenceAngleTorsionForce::getAddition3(std::vector<double>& xVector, std::vector<double>& yVector) const {
	xVector[0] += yVector[0];
	xVector[1] += yVector[1];
	xVector[2] += yVector[2];
}

void AmoebaReferenceAngleTorsionForce::getSubstraction3(std::vector<double>& xVector, std::vector<double>& yVector, std::vector<double>& zVector) const {
	zVector[0] = xVector[0] - yVector[0];
	zVector[1] = xVector[1] - yVector[1];
	zVector[2] = xVector[2] - yVector[2];
}

void AmoebaReferenceAngleTorsionForce::getSubstraction3(std::vector<double>& xVector, std::vector<double>& yVector) const {
	xVector[0] -= yVector[0];
	xVector[1] -= yVector[1];
	xVector[2] -= yVector[2];
}

void AmoebaReferenceAngleTorsionForce::getScale3(std::vector<double>& xVector, double scale, std::vector<double>& yVector) const {
	yVector[0] = xVector[0] * scale;
	yVector[1] = xVector[1] * scale;
	yVector[2] = xVector[2] * scale;
}

void AmoebaReferenceAngleTorsionForce::getScale3(std::vector<double>& xVector, double scale) const {
	xVector[0] *= scale;
	xVector[1] *= scale;
	xVector[2] *= scale;
}

double AmoebaReferenceAngleTorsionForce::calculateForceAndEnergy(int numAngleTorsions, std::vector<OpenMM::Vec3>& posData,
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
		std::vector<OpenMM::Vec3>& forceData) const {
	double energy      = 0.0;
	for (unsigned int ii = 0; ii < static_cast<unsigned int>(numAngleTorsions); ii++) {
		int particle1Index      = particle1[ii];
		int particle2Index      = particle2[ii];
		int particle3Index      = particle3[ii];
		int particle4Index      = particle4[ii];
		double cbaAngle     = angleCBAParameters[ii];
		double dcbAngle     = angleDCBParameters[ii];
		double K1      = k1Parameters[ii];
		double K2      = k2Parameters[ii];
		double K3      = k3Parameters[ii];
		double K4      = k4Parameters[ii];
		double K5      = k5Parameters[ii];
		double K6      = k6Parameters[ii];
		Vec3 forces[4];
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

double AmoebaReferenceAngleTorsionForce::calculateAngleTorsionIxn(const Vec3& positionAtomA, const Vec3& positionAtomB,
		const Vec3& positionAtomC, const Vec3& positionAtomD,
		double angleCBA,      double angleDCB,
		double k1Parameter, double k2Parameter, double k3Parameter,
		double k4Parameter, double k5Parameter, double k6Parameter,
		Vec3* forces) const {

	static const double zero          = 0.0;
        static const double one           = 1.0;
        static const double two           = 2.0;
        static const double three         = 3.0;

	enum { A = 0, B, C, D, LastAtomIndex };
	enum { BA = 0, CB, DC, CA, DB, BAxCB, CBxDC, TxU, DEDT, DEDU, TxBA, TxCB, UxCB, UxDC, DEDTxCB, CAxDEDT, DEDUxDC, DEDTxBA, DBxDEDU, DEDUxCB, LastDeltaIndex };

	std::vector<double> deltaR[LastDeltaIndex];
	for (unsigned int ii = 0; ii < LastDeltaIndex; ii++) {
		deltaR[ii].resize(3);
	}
	AmoebaReferenceForce::loadDeltaR(positionAtomA, positionAtomB, deltaR[BA]);
	AmoebaReferenceForce::loadDeltaR(positionAtomB, positionAtomC, deltaR[CB]);
	AmoebaReferenceForce::loadDeltaR(positionAtomC, positionAtomD, deltaR[DC]);

	double rBA2 = AmoebaReferenceForce::getNormSquared3(deltaR[BA]);
	double rCB2 = AmoebaReferenceForce::getNormSquared3(deltaR[CB]);
	double rDC2 = AmoebaReferenceForce::getNormSquared3(deltaR[DC]);
	AmoebaReferenceForce::getCrossProduct(deltaR[BA], deltaR[CB], deltaR[BAxCB]);
	AmoebaReferenceForce::getCrossProduct(deltaR[CB], deltaR[DC], deltaR[CBxDC]);
	AmoebaReferenceForce::getCrossProduct(deltaR[BAxCB], deltaR[CBxDC], deltaR[TxU]);
	double  rT2 = AmoebaReferenceForce::getNormSquared3(deltaR[BAxCB]);
	double  rU2 = AmoebaReferenceForce::getNormSquared3(deltaR[CBxDC]);
	double  rTrU = SQRT(rT2 * rU2);

	if (rTrU <= zero) {
		return zero;
	}

	AmoebaReferenceForce::loadDeltaR(positionAtomA, positionAtomC, deltaR[CA]);
	AmoebaReferenceForce::loadDeltaR(positionAtomB, positionAtomD, deltaR[DB]);
	double  rCB = SQRT(rCB2);
	double  cosine = AmoebaReferenceForce::getDotProduct3(deltaR[BAxCB], deltaR[CBxDC]) / rTrU;
	double  sine = AmoebaReferenceForce::getDotProduct3(deltaR[CB], deltaR[TxU]) / (rCB * rTrU);

        double  c1 = one;
        double  s1 = zero;
        double  c2 = - one;
        double  s2 = zero;
        double  c3 = one;
        double  s3 = zero;
        double  cosine2 = cosine * cosine - sine * sine;
        double  sine2 = two * sine * cosine;
        double  cosine3 = cosine * cosine2 - sine * sine2;
        double  sine3 = cosine * sine2 + sine * cosine2;
        double  phi1 = one + (cosine * c1 + sine * s1);
        double  phi2 = one + (cosine2 * c2 + sine2 * s2);
        double  phi3 = one + (cosine3 * c3 + sine3 * s3);
        double  dphi1 = cosine * s1 - sine * c1;
        double  dphi2 = two * (cosine2 * s2 - sine2 * c2);
        double  dphi3 = three * (cosine3 * s3 - sine3 * c3);

        double  v1 = k1Parameter;
	double  v2 = k2Parameter;
	double  v3 = k3Parameter;
	double  dot = AmoebaReferenceForce::getDotProduct3(deltaR[BA], deltaR[CB]);
	double  cosang = - dot / SQRT(rBA2 * rCB2);
	double  angle = RADIAN * ACOS(cosang);
	double  dt = angle - angleCBA;
	double  e1 = dt * (v1 * phi1 + v2 * phi2 + v3 * phi3);
	double  dedphi = dt * (v1 * dphi1 + v2 * dphi2 + v3 * dphi3);
	double  ddt = RADIAN * (v1 * phi1 + v2 * phi2 + v3 * phi3);

	AmoebaReferenceForce::getCrossProduct(deltaR[BAxCB], deltaR[CB], deltaR[DEDT]);
	double tempScale = dedphi / (rT2 * rCB);
	getScale3(deltaR[DEDT], tempScale);

	AmoebaReferenceForce::getCrossProduct(deltaR[CB], deltaR[CBxDC], deltaR[DEDU]);
	tempScale = dedphi / (rU2 * rCB);
	getScale3(deltaR[DEDU], tempScale);

	std::vector<double> subForce[LastAtomIndex];
	for (int ii = 0; ii < LastAtomIndex; ii++) {
		subForce[ii].resize(3);
	}

	double  terma = - ddt / (rBA2 * SQRT(rT2));
	double  termc = ddt / (rCB2 * SQRT(rT2));
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
	double  e2 = dt * (v1 * phi1 + v2 * phi2 + v3 * phi3);
	dedphi = dt * (v1 * dphi1 + v2 * dphi2 + v3 * dphi3);
	ddt = RADIAN * (v1 * phi1 + v2 * phi2 + v3 * phi3);

	AmoebaReferenceForce::getCrossProduct(deltaR[BAxCB], deltaR[CB], deltaR[DEDT]);
	tempScale = dedphi / (rT2 * rCB);
	getScale3(deltaR[DEDT], tempScale);
	AmoebaReferenceForce::getCrossProduct(deltaR[CB], deltaR[CBxDC], deltaR[DEDU]);
	tempScale = dedphi / (rU2 * rCB);
	getScale3(deltaR[DEDU], tempScale);

	double  termb = - ddt / (rCB2 * SQRT(rU2));
	double  termd = ddt / (rDC2 * SQRT(rU2));
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

	double  e = e1 + e2;
	for (int jj = 0; jj < LastAtomIndex; jj++) {
		forces[jj][0] = subForce[jj][0];
		forces[jj][1] = subForce[jj][1];
		forces[jj][2] = subForce[jj][2];
	}
	return e;
}

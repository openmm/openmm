#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>

using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerAmoebaReferenceKernelFactories();

const double TOL = 1e-4;
#define RADIAN            57.29577951308

static void crossProductVector3(double* vectorX, double* vectorY, double* vectorZ) {

	vectorZ[0]  = vectorX[1]*vectorY[2] - vectorX[2]*vectorY[1];
	vectorZ[1]  = vectorX[2]*vectorY[0] - vectorX[0]*vectorY[2];
	vectorZ[2]  = vectorX[0]*vectorY[1] - vectorX[1]*vectorY[0];
	return;
}

static double dotVector3(double* vectorX, double* vectorY) {
	return vectorX[0]*vectorY[0] + vectorX[1]*vectorY[1] + vectorX[2]*vectorY[2];
}

static void computeAmoebaAngleTorsionForce(int Index,  std::vector<Vec3>& positions, AmoebaAngleTorsionForce& amoebaAngleTorsionForce,
		std::vector<Vec3>& forces, double* energy) {
	int particle1, particle2, particle3, particle4;
	double cbaAngle, dcbAngle, k1, k2, k3, k4, k5, k6;
	amoebaAngleTorsionForce.getAngleTorsionParameters(Index, particle1, particle2, particle3, particle4, cbaAngle, dcbAngle,
			k1, k2, k3, k4, k5, k6);

	enum { A = 0, B, C, D, LastAtomIndex };
	enum { BA = 0, CB, DC, CA, DB, BAxCB, CBxDC, TxU, DEDT, DEDU, TxBA, TxCB, UxCB, UxDC, DEDTxCB, CAxDEDT, DEDUxDC, DEDTxBA, DBxDEDU, DEDUxCB, LastDeltaIndex };

	double deltaR[LastDeltaIndex][3];
	double rBA2 = 0.0;
	double rCB2 = 0.0;
	double rDC2 = 0.0;

	for (int ii = 0; ii < 3; ii++) {
		deltaR[BA][ii]  = positions[particle2][ii] - positions[particle1][ii];
		rBA2           += deltaR[BA][ii]*deltaR[BA][ii];

		deltaR[CB][ii]  = positions[particle3][ii] - positions[particle2][ii];
		rCB2           += deltaR[CB][ii]*deltaR[CB][ii];

		deltaR[DC][ii] =  positions[particle4][ii] - positions[particle3][ii];
		rDC2           += deltaR[DC][ii]*deltaR[DC][ii];

		deltaR[CA][ii] =  positions[particle3][ii] - positions[particle1][ii];

		deltaR[DB][ii] =  positions[particle4][ii] - positions[particle2][ii];
	}

	crossProductVector3(deltaR[BA], deltaR[CB], deltaR[BAxCB]);
	crossProductVector3(deltaR[CB], deltaR[DC], deltaR[CBxDC]);
	crossProductVector3(deltaR[BAxCB], deltaR[CBxDC], deltaR[TxU]);
	double rT2 = dotVector3(deltaR[BAxCB], deltaR[BAxCB]);
	double rU2 = dotVector3(deltaR[CBxDC], deltaR[CBxDC]);
	double rTrU = sqrt(rT2 * rU2);

	if (rTrU <= 0.0) {
		return;
	}

	double rCB = sqrt(rCB2);
	double cosine = dotVector3(deltaR[BAxCB], deltaR[CBxDC]) / rTrU;
	double sine = dotVector3(deltaR[CB], deltaR[TxU]) / (rCB * rTrU);

	double c1 = 1.0;
	double s1 = 0.0;
	double c2 = -1.0;
	double s2 = 0.0;
	double c3 = 1.0;
	double s3 = 0.0;
	double cosine2 = cosine * cosine - sine * sine;
	double sine2 = 2.0 * sine * cosine;
	double cosine3 = cosine * cosine2 - sine * sine2;
	double sine3 = cosine * sine2 + sine * cosine2;
	double phi1 = 1.0 + (cosine * c1 + sine * s1);
	double phi2 = 1.0 + (cosine2 * c2 + sine2 * s2);
	double phi3 = 1.0 + (cosine3 * c3 + sine3 * s3);
	double dphi1 = cosine * s1 - sine * c1;
	double dphi2 = 2.0 * (cosine2 * s2 - sine2 * c2);
	double dphi3 = 3.0 * (cosine3 * s3 - sine3 * c3);

	double v1 = k1;
	double v2 = k2;
	double v3 = k3;
	double dot = dotVector3(deltaR[BA], deltaR[CB]);
	double cosang = - dot / sqrt(rBA2 * rCB2);
	double angle = RADIAN * acos(cosang);
	double dt = angle - cbaAngle;
	double e1 = dt * (v1 * phi1 + v2 * phi2 + v3 * phi3);
	double dedphi = dt * (v1 * dphi1 + v2 * dphi2 + v3 * dphi3);
	double ddt = RADIAN * (v1 * phi1 + v2 * phi2 + v3 * phi3);

	crossProductVector3(deltaR[BAxCB], deltaR[CB], deltaR[DEDT]);
	double tempScale = dedphi / (rT2 * rCB);
	for (int ii = 0; ii < 3; ii++) {
		deltaR[DEDT][ii] *= tempScale;
	}
	crossProductVector3(deltaR[CB], deltaR[CBxDC], deltaR[DEDU]);
	tempScale = dedphi / (rU2 * rCB);
	for (int ii = 0; ii < 3; ii++) {
		deltaR[DEDU][ii] *= tempScale;
	}
	double subForce[LastAtomIndex][3];
	double terma = - ddt / (rBA2 * sqrt(rT2));
	double termc = ddt / (rCB2 * sqrt(rT2));
	crossProductVector3(deltaR[BAxCB], deltaR[BA], deltaR[TxBA]);
	for (int ii = 0; ii < 3; ii++) {
		deltaR[TxBA][ii] *= terma;
	}
	crossProductVector3(deltaR[DEDT], deltaR[CB], deltaR[DEDTxCB]);
	for (int ii = 0; ii < 3; ii++) {
		subForce[A][ii] = deltaR[TxBA][ii] + deltaR[DEDTxCB][ii];
	}

	crossProductVector3(deltaR[BAxCB], deltaR[CB], deltaR[TxCB]);
	for (int ii = 0; ii < 3; ii++) {
		deltaR[TxCB][ii] *= termc;
	}
	for (int ii = 0; ii < 3; ii++) {
		subForce[B][ii] = deltaR[TxCB][ii] - deltaR[TxBA][ii];
	}
	crossProductVector3(deltaR[CA], deltaR[DEDT], deltaR[CAxDEDT]);
	for (int ii = 0; ii < 3; ii++) {
		subForce[B][ii] += deltaR[CAxDEDT][ii];
	}
	crossProductVector3(deltaR[DEDU], deltaR[DC], deltaR[DEDUxDC]);
	for (int ii = 0; ii < 3; ii++) {
		subForce[B][ii] += deltaR[DEDUxDC][ii];
	}

	crossProductVector3(deltaR[DEDT], deltaR[BA], deltaR[DEDTxBA]);
	for (int ii = 0; ii < 3; ii++) {
		subForce[C][ii] = deltaR[DEDTxBA][ii] - deltaR[TxCB][ii];
	}
	crossProductVector3(deltaR[DB], deltaR[DEDU], deltaR[DBxDEDU]);
	for (int ii = 0; ii < 3; ii++) {
		subForce[C][ii] += deltaR[DBxDEDU][ii];
	}

	crossProductVector3(deltaR[DEDU], deltaR[CB], subForce[D]);

	v1 = k4;
	v2 = k5;
	v3 = k6;
	dot = dotVector3(deltaR[CB], deltaR[DC]);
	cosang = - dot / sqrt(rCB2 * rDC2);
	angle = RADIAN * acos(cosang);
	dt = angle - dcbAngle;
	double e2 = dt * (v1 * phi1 + v2 * phi2 + v3 * phi3);
	dedphi = dt * (v1 * dphi1 + v2 * dphi2 + v3 * dphi3);
	ddt = RADIAN * (v1 * phi1 + v2 * phi2 + v3 * phi3);
	crossProductVector3(deltaR[BAxCB], deltaR[CB], deltaR[DEDT]);
	tempScale = dedphi / (rT2 * rCB);
	for (int ii = 0; ii < 3; ii++) {
		deltaR[DEDT][ii] *= tempScale;
	}
	crossProductVector3(deltaR[CB], deltaR[CBxDC], deltaR[DEDU]);
	tempScale = dedphi / (rU2 * rCB);
	for (int ii = 0; ii < 3; ii++) {
		deltaR[DEDU][ii] *= tempScale;
	}
	double termb = - ddt / (rCB2 * sqrt(rU2));
	double termd = ddt / (rDC2 * sqrt(rU2));
	crossProductVector3(deltaR[DEDT], deltaR[CB], deltaR[DEDTxCB]);
	for (int ii = 0; ii < 3; ii++) {
		subForce[A][ii] += deltaR[DEDTxCB][ii];
	}

	crossProductVector3(deltaR[CBxDC], deltaR[CB], deltaR[UxCB]);
	for (int ii = 0; ii < 3; ii++) {
		deltaR[UxCB][ii] *= termb;
	}
	crossProductVector3(deltaR[CA], deltaR[DEDT], deltaR[CAxDEDT]);
	crossProductVector3(deltaR[DEDU], deltaR[DC], deltaR[DEDUxDC]);
	for (int ii = 0; ii < 3; ii++) {
		subForce[B][ii] += deltaR[UxCB][ii];
		subForce[B][ii] += deltaR[CAxDEDT][ii];
		subForce[B][ii] += deltaR[DEDUxDC][ii];
	}

	for (int ii = 0; ii < 3; ii++) {
		subForce[C][ii] -= deltaR[UxCB][ii];
	}
	crossProductVector3(deltaR[CBxDC], deltaR[DC], deltaR[UxDC]);
	for (int ii = 0; ii < 3; ii++) {
		deltaR[UxDC][ii] *= termd;
	}
	for (int ii = 0; ii < 3; ii++) {
		subForce[C][ii] += deltaR[UxDC][ii];
	}
	crossProductVector3(deltaR[DEDT], deltaR[BA], deltaR[DEDTxBA]);
	for (int ii = 0; ii < 3; ii++) {
		subForce[C][ii] += deltaR[DEDTxBA][ii];
	}
	crossProductVector3(deltaR[DB], deltaR[DEDU], deltaR[DBxDEDU]);
	for (int ii = 0; ii < 3; ii++) {
		subForce[C][ii] += deltaR[DBxDEDU][ii];
	}

	for (int ii = 0; ii < 3; ii++) {
		subForce[D][ii] -= deltaR[UxDC][ii];
	}
	crossProductVector3(deltaR[DEDU], deltaR[CB], deltaR[DEDUxCB]);
	for (int ii = 0; ii < 3; ii++) {
		subForce[D][ii] += deltaR[DEDUxCB][ii];
	}

	double e = e1 + e2;

        forces[particle1][0]       -= subForce[A][0];
	forces[particle1][1]       -= subForce[A][1];
	forces[particle1][2]       -= subForce[A][2];

	forces[particle2][0]       -= subForce[B][0];
	forces[particle2][1]       -= subForce[B][1];
	forces[particle2][2]       -= subForce[B][2];

	forces[particle3][0]       -= subForce[C][0];
	forces[particle3][1]       -= subForce[C][1];
	forces[particle3][2]       -= subForce[C][2];

	forces[particle4][0]       -= subForce[D][0];
	forces[particle4][1]       -= subForce[D][1];
	forces[particle4][2]       -= subForce[D][2];

	*energy += e;
}

static void computeAmoebaAngleTorsionForces(Context& context, AmoebaAngleTorsionForce& amoebaAngleTorsionForce,
		std::vector<Vec3>& expectedForces, double* expectedEnergy) {
	// get positions and zero forces

	State state                 = context.getState(State::Positions);
	std::vector<Vec3> positions = state.getPositions();
	expectedForces.resize(positions.size());

	for (unsigned int ii = 0; ii < expectedForces.size(); ii++) {
		expectedForces[ii][0] = expectedForces[ii][1] = expectedForces[ii][2] = 0.0;
	}

	// calculates forces/energy

	*expectedEnergy = 0.0;
	for (int ii = 0; ii < amoebaAngleTorsionForce.getNumAngleTorsions(); ii++) {
		computeAmoebaAngleTorsionForce(ii, positions, amoebaAngleTorsionForce, expectedForces, expectedEnergy);
	}
}

void compareWithExpectedForceAndEnergy(Context& context, AmoebaAngleTorsionForce& amoebaAngleTorsionForce,
		double tolerance, const std::string& idString) {
	std::vector<Vec3> expectedForces;
	double expectedEnergy;
	computeAmoebaAngleTorsionForces(context, amoebaAngleTorsionForce, expectedForces, &expectedEnergy);
	State state                      = context.getState(State::Forces | State::Energy);

	const std::vector<Vec3> forces   = state.getForces();
	for (unsigned int ii = 0; ii < forces.size(); ii++) {
		ASSERT_EQUAL_VEC(expectedForces[ii], forces[ii], tolerance);
	}
	ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), tolerance);
}

void testOneAngleTorsion() {
	System system;
	int numberOfParticles = 4;
	for (int ii = 0; ii < numberOfParticles; ii++) {
		system.addParticle(1.0);
	}

	LangevinIntegrator integrator(0.0, 0.1, 0.01);

	AmoebaAngleTorsionForce* amoebaAngleTorsionForce = new AmoebaAngleTorsionForce();
	double cbaAngle = 108.30;
	double dcbAngle = 108.30;
	double k1 = 56.4;
	double k2 = -30.1;
	double k3 = 15.9;
	double k4 = 97.1;
	double k5 = -66.8;
	double k6 = 14.8;

	amoebaAngleTorsionForce->addAngleTorsion(0, 1, 2, 3, cbaAngle, dcbAngle, k1, k2, k3, k4, k5, k6);
	system.addForce(amoebaAngleTorsionForce);
	ASSERT(!amoebaAngleTorsionForce->usesPeriodicBoundaryConditions());
	ASSERT(!system.usesPeriodicBoundaryConditions());
	Context context(system, integrator, Platform::getPlatformByName("Reference"));

	std::vector<Vec3> positions(numberOfParticles);

	positions[0] = Vec3(-2.202925,   -1.205586,    1.055576);
	positions[1] = Vec3(-2.108630,   -0.382859,   -0.084112);
	positions[2] = Vec3(-1.142172,   -0.983465,   -0.890852);
	positions[3] = Vec3(0.171092,   -0.517551,   -0.531557);

	context.setPositions(positions);
	compareWithExpectedForceAndEnergy(context, *amoebaAngleTorsionForce, TOL, "testOneAngleTorsion");

	amoebaAngleTorsionForce->setAngleTorsionParameters(0, 0, 1, 2, 3, 1.1*cbaAngle, 1.2*dcbAngle, 1.4*k1, 1.4*k2, 1.4*k3, 1.4*k4, 1.4*k5, 1.4*k6);
	bool exceptionThrown = false;
	try {
		// This should throw an exception.
		compareWithExpectedForceAndEnergy(context, *amoebaAngleTorsionForce, TOL, "testOneAngleTorsion");
	}

	catch (std::exception ex) {
		exceptionThrown = true;
	}

	ASSERT(exceptionThrown);
	amoebaAngleTorsionForce->updateParametersInContext(context);
	compareWithExpectedForceAndEnergy(context, *amoebaAngleTorsionForce, TOL, "testOneAngleTorsion");
}

int main(int numberOfArguments, char* argv[]) {

	try {
		std::cout << "TestReferenceAmoebaAngleTorsionForce running test..." << std::endl;
		registerAmoebaReferenceKernelFactories();
		testOneAngleTorsion();
	}
	catch(const std::exception& e) {
		std::cout << "exception: " << e.what() << std::endl;
		std::cout << "FAIL - ERROR.  Test failed." << std::endl;
		return 1;
	}
		std::cout << "Done" << std::endl;
	return 0;
}

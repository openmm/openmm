#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/AmoebaAngleTorsionForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

extern "C" void registerAmoebaSerializationProxies();

void testSerialization() {
	// Create a Force.

	AmoebaAngleTorsionForce force1;
	force1.setForceGroup(3);
	force1.addAngleTorsion(0, 1, 2, 3, 64.6, 100.0, 3.3, 90.9, 489.2, 2.3, 3.2, 490.0);
	force1.addAngleTorsion(4, 3, 2, 3, 6.6, 50.0, 5.3, 9.9, 899.2, 2.4, 3.2, 40.0);
	force1.addAngleTorsion(5, 8, 7, 6, 4.67, 18.0, 3.8, 190.9, 389.2, 2.9, 3.3, 90.0);

	// Serialize and then deserialize it.

	stringstream buffer;
	XmlSerializer::serialize<AmoebaAngleTorsionForce>(&force1, "Force", buffer);
	AmoebaAngleTorsionForce* copy = XmlSerializer::deserialize<AmoebaAngleTorsionForce>(buffer);

	// Compare the two forces to see if they are identical.

	AmoebaAngleTorsionForce& force2 = *copy;
	ASSERT_EQUAL(force1.getForceGroup(), force2.getForceGroup());
	ASSERT_EQUAL(force1.getNumAngleTorsions(), force2.getNumAngleTorsions());
	for (unsigned int ii = 0; ii < static_cast<unsigned int>(force1.getNumAngleTorsions()); ii++) {
		int p11, p12, p13, p14;
		int p21, p22, p23, p24;
		double angleCBA1, angleDCB1;
		double angleCBA2, angleDCB2;
		double k11, k12, k13, k14, k15, k16;
		double k21, k22, k23, k24, k25, k26;
		force1.getAngleTorsionParameters(ii, p11, p12, p13, p14, angleCBA1, angleDCB1, k11, k12, k13, k14, k15, k16);
		force2.getAngleTorsionParameters(ii, p21, p22, p23, p24, angleCBA2, angleDCB2, k21, k22, k23, k24, k25, k26);

		ASSERT_EQUAL(p11, p21);
		ASSERT_EQUAL(p12, p22);
		ASSERT_EQUAL(p13, p23);
		ASSERT_EQUAL(p14, p24);
		ASSERT_EQUAL(angleCBA1, angleCBA2);
		ASSERT_EQUAL(angleDCB1, angleDCB2);
		ASSERT_EQUAL(k11, k21);
		ASSERT_EQUAL(k12, k22);
		ASSERT_EQUAL(k13, k23);
		ASSERT_EQUAL(k14, k24);
		ASSERT_EQUAL(k15, k25);
		ASSERT_EQUAL(k16, k26);
	}
}

int main() {
	try {
		registerAmoebaSerializationProxies();
		testSerialization();
	}

	catch(const exception& e) {
		cout << "exception: " << e.what() << endl;
		return 1;
	}

	cout << "Done" << endl;
	return 0;
}

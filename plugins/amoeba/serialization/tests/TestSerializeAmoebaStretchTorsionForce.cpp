#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/AmoebaStretchTorsionForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

extern "C" void registerAmoebaSerializationProxies();

void testSerialization() {
	// Create a Force.

	AmoebaStretchTorsionForce force1;
	force1.setForceGroup(3);
	force1.addStretchTorsion(0, 1, 2, 3, 1.0, 1.2, 1.5, 64.6, 100.0, 3.3, 90.9, 489.2, 2.3, 3.2, 490.0, 4.5);
	force1.addStretchTorsion(4, 3, 2, 3, 1.4, 1.7, 2.5, 6.6, 50.0, 5.3, 9.9, 899.2, 2.4, 3.2, 40.0, 0.5);
	force1.addStretchTorsion(5, 8, 7, 6, 1.9, 4.2, 3.5, 4.67, 18.0, 3.8, 190.9, 389.2, 2.9, 3.3, 90.0, 84.5);

	// Serialize and then deserialize it.

	stringstream buffer;
	XmlSerializer::serialize<AmoebaStretchTorsionForce>(&force1, "Force", buffer);
	AmoebaStretchTorsionForce* copy = XmlSerializer::deserialize<AmoebaStretchTorsionForce>(buffer);

	// Compare the two forces to see if they are identical.
	AmoebaStretchTorsionForce& force2 = *copy;
	ASSERT_EQUAL(force1.getForceGroup(), force2.getForceGroup());
	ASSERT_EQUAL(force1.getNumStretchTorsions(), force2.getNumStretchTorsions());
	for (unsigned int ii = 0; ii < static_cast<unsigned int>(force1.getNumStretchTorsions()); ii++) {
		int p11, p12, p13, p14;
		int p21, p22, p23, p24;
		double lengthBA1, lengthCB1, lengthDC1;
		double lengthBA2, lengthCB2, lengthDC2;
		double k11, k12, k13, k14, k15, k16, k17, k18, k19;
		double k21, k22, k23, k24, k25, k26, k27, k28, k29;

		force1.getStretchTorsionParameters(ii, p11, p12, p13, p14, lengthBA1, lengthCB1, lengthDC1, k11, k12, k13, k14, k15, k16, k17, k18, k19);
		force2.getStretchTorsionParameters(ii, p21, p22, p23, p24, lengthBA2, lengthCB2, lengthDC2, k21, k22, k23, k24, k25, k26, k27, k28, k29);

		ASSERT_EQUAL(p11, p21);
		ASSERT_EQUAL(p12, p22);
		ASSERT_EQUAL(p13, p23);
		ASSERT_EQUAL(p14, p24);
		ASSERT_EQUAL(lengthBA1, lengthBA2);
		ASSERT_EQUAL(lengthCB1, lengthCB2);
		ASSERT_EQUAL(lengthDC1, lengthDC2);
		ASSERT_EQUAL(k11, k21);
		ASSERT_EQUAL(k12, k22);
		ASSERT_EQUAL(k13, k23);
		ASSERT_EQUAL(k14, k24);
		ASSERT_EQUAL(k15, k25);
		ASSERT_EQUAL(k16, k26);
		ASSERT_EQUAL(k17, k27);
		ASSERT_EQUAL(k18, k28);
		ASSERT_EQUAL(k19, k29);
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

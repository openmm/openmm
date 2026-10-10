/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * Metal Platform tests: copyright (c) 2026 Chun-Chi Hung.                     *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 * This program is distributed WITHOUT ANY WARRANTY; without even the        *
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. *
 * See the GNU Lesser General Public License for more details.                *
 * You should have received a copy of the GNU Lesser General Public License  *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#include "openmm/Context.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/Platform.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/internal/AssertionUtilities.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

/** @brief Check dynamically loaded registration and single-GPU property validation. */
int main() {
    try {
        Platform::loadPluginLibrary(METAL_LIBRARY_PATH);
        Platform* metal = nullptr;
        for (int i = 0; i < Platform::getNumPlatforms(); i++)
            if (Platform::getPlatform(i).getName() == "Metal")
                metal = &Platform::getPlatform(i);
        if (metal == nullptr) {
            cout << "No supported Metal device is available" << endl;
            return 77;
        }
        ASSERT(!metal->supportsDoublePrecision());
        ASSERT_EQUAL(1, metal->getDevices().size());
        ASSERT(metal->getDevices({{"Precision", "double"}}).empty());
        ASSERT(metal->getDevices({{"DeviceIndex", "0,0"}}).empty());
        ASSERT(metal->getDevices({{"UseCpuPme", "true"}}).empty());
        System system;
        system.addParticle(1);
        CustomExternalForce* force = new CustomExternalForce("x*x+y*y+z*z");
        force->addParticle(0);
        system.addForce(force);
        for (auto property : {make_pair("Precision", "mixed"), make_pair("DeviceIndex", "0,0"),
                make_pair("UseCpuPme", "true")}) {
            bool rejected = false;
            try {
                VerletIntegrator integrator(0.001);
                Context context(system, integrator, *metal, {{property.first, property.second}});
            }
            catch (const OpenMMException&) {
                rejected = true;
            }
            ASSERT(rejected);
        }
        VerletIntegrator integrator(0.001);
        Context context(system, integrator, *metal);
        context.setPositions({Vec3(1, 2, 3)});
        ASSERT_EQUAL("single", metal->getPropertyValue(context, "Precision"));
        ASSERT_EQUAL("false", metal->getPropertyValue(context, "UseCpuPme"));
        ASSERT_EQUAL_TOL(14, context.getState(State::Energy).getPotentialEnergy(), 1e-6);
        cout << "Metal Platform loading and properties passed" << endl;
    }
    catch (const exception& error) {
        cerr << error.what() << endl;
        return 1;
    }
    return 0;
}

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/AndersenThermostat.h"
#include "openmm/CMAPTorsionForce.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/GBVIForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/MonteCarloBarostat.h"
#include "openmm/NonbondedForce.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/RBTorsionForce.h"
#include "openmm/System.h"
#include "openmm/serialization/SerializationProxy.h"
#include "openmm/serialization/AndersenThermostatProxy.h"
#include "openmm/serialization/CMAPTorsionForceProxy.h"
#include "openmm/serialization/CMMotionRemoverProxy.h"
#include "openmm/serialization/GBSAOBCForceProxy.h"
#include "openmm/serialization/GBVIForceProxy.h"
#include "openmm/serialization/HarmonicAngleForceProxy.h"
#include "openmm/serialization/HarmonicBondForceProxy.h"
#include "openmm/serialization/MonteCarloBarostatProxy.h"
#include "openmm/serialization/NonbondedForceProxy.h"
#include "openmm/serialization/PeriodicTorsionForceProxy.h"
#include "openmm/serialization/RBTorsionForceProxy.h"
#include "openmm/serialization/SystemProxy.h"

#if defined(WIN32)
    #include <windows.h>
    extern "C" void initOpenMMPlugin();
    BOOL WINAPI DllMain(HANDLE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved) {
        if (ul_reason_for_call == DLL_PROCESS_ATTACH)
            registerSerializationProxies();
        return TRUE;
    }
#else
    extern "C" void __attribute__((constructor)) registerSerializationProxies();
#endif

using namespace OpenMM;

extern "C" void registerSerializationProxies() {
    SerializationProxy::registerProxy(typeid(AndersenThermostat), new AndersenThermostatProxy());
    SerializationProxy::registerProxy(typeid(CMAPTorsionForce), new CMAPTorsionForceProxy());
    SerializationProxy::registerProxy(typeid(CMMotionRemover), new CMMotionRemoverProxy());
    SerializationProxy::registerProxy(typeid(GBSAOBCForce), new GBSAOBCForceProxy());
    SerializationProxy::registerProxy(typeid(GBVIForce), new GBVIForceProxy());
    SerializationProxy::registerProxy(typeid(HarmonicAngleForce), new HarmonicAngleForceProxy());
    SerializationProxy::registerProxy(typeid(HarmonicBondForce), new HarmonicBondForceProxy());
    SerializationProxy::registerProxy(typeid(MonteCarloBarostat), new MonteCarloBarostatProxy());
    SerializationProxy::registerProxy(typeid(NonbondedForce), new NonbondedForceProxy());
    SerializationProxy::registerProxy(typeid(PeriodicTorsionForce), new PeriodicTorsionForceProxy());
    SerializationProxy::registerProxy(typeid(RBTorsionForce), new RBTorsionForceProxy());
    SerializationProxy::registerProxy(typeid(System), new SystemProxy());
}
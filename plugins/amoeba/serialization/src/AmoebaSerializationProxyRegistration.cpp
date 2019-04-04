/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
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

#ifdef WIN32
#include <windows.h>
#include <sstream>
#else
#include <dlfcn.h>
#include <dirent.h>
#include <cstdlib>
#endif

#include "openmm/OpenMMException.h"

#include "openmm/AmoebaGeneralizedKirkwoodForce.h"
#include "openmm/AmoebaBondForce.h"
#include "openmm/AmoebaAngleForce.h"
#include "openmm/AmoebaInPlaneAngleForce.h"
#include "openmm/AmoebaMultipoleForce.h"
#include "openmm/AmoebaOutOfPlaneBendForce.h"
#include "openmm/AmoebaPiTorsionForce.h"
#include "openmm/AmoebaStretchBendForce.h"
#include "openmm/AmoebaTorsionTorsionForce.h"
#include "openmm/AmoebaVdwForce.h"
#include "openmm/AmoebaWcaDispersionForce.h"
#include "openmm/HippoNonbondedForce.h"

#include "openmm/serialization/SerializationProxy.h"

#include "openmm/serialization/AmoebaGeneralizedKirkwoodForceProxy.h"
#include "openmm/serialization/AmoebaBondForceProxy.h"
#include "openmm/serialization/AmoebaAngleForceProxy.h"
#include "openmm/serialization/AmoebaInPlaneAngleForceProxy.h"
#include "openmm/serialization/AmoebaMultipoleForceProxy.h"
#include "openmm/serialization/AmoebaOutOfPlaneBendForceProxy.h"
#include "openmm/serialization/AmoebaPiTorsionForceProxy.h"
#include "openmm/serialization/AmoebaStretchBendForceProxy.h"
#include "openmm/serialization/AmoebaTorsionTorsionForceProxy.h"
#include "openmm/serialization/AmoebaVdwForceProxy.h"
#include "openmm/serialization/AmoebaWcaDispersionForceProxy.h"
#include "openmm/serialization/HippoNonbondedForceProxy.h"

#if defined(WIN32)
    #include <windows.h>
    extern "C" OPENMM_EXPORT_AMOEBA void registerAmoebaSerializationProxies();
    BOOL WINAPI DllMain(HANDLE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved) {
        if (ul_reason_for_call == DLL_PROCESS_ATTACH)
            registerAmoebaSerializationProxies();
        return TRUE;
    }
#else
    extern "C" void __attribute__((constructor)) registerAmoebaSerializationProxies();
#endif

using namespace OpenMM;

extern "C" OPENMM_EXPORT_AMOEBA void registerAmoebaSerializationProxies() {
    SerializationProxy::registerProxy(typeid(AmoebaGeneralizedKirkwoodForce),         new AmoebaGeneralizedKirkwoodForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaBondForce),                new AmoebaBondForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaAngleForce),               new AmoebaAngleForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaInPlaneAngleForce),        new AmoebaInPlaneAngleForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaMultipoleForce),                   new AmoebaMultipoleForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaOutOfPlaneBendForce),              new AmoebaOutOfPlaneBendForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaPiTorsionForce),                   new AmoebaPiTorsionForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaStretchBendForce),                 new AmoebaStretchBendForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaTorsionTorsionForce),              new AmoebaTorsionTorsionForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaVdwForce),                         new AmoebaVdwForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaWcaDispersionForce),               new AmoebaWcaDispersionForceProxy());
    SerializationProxy::registerProxy(typeid(HippoNonbondedForce),                    new HippoNonbondedForceProxy());
}

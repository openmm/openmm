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
#include "openmm/AmoebaHarmonicBondForce.h"
#include "openmm/AmoebaHarmonicAngleForce.h"
#include "openmm/AmoebaHarmonicInPlaneAngleForce.h"
#include "openmm/AmoebaMultipoleForce.h"
#include "openmm/AmoebaOutOfPlaneBendForce.h"
#include "openmm/AmoebaPiTorsionForce.h"
#include "openmm/AmoebaStretchBendForce.h"
#include "openmm/AmoebaTorsionForce.h"
#include "openmm/AmoebaTorsionTorsionForce.h"
#include "openmm/AmoebaUreyBradleyForce.h"
#include "openmm/AmoebaVdwForce.h"
#include "openmm/AmoebaWcaDispersionForce.h"

#include "openmm/serialization/SerializationProxy.h"

#include "openmm/serialization/AmoebaGeneralizedKirkwoodForceProxy.h"
#include "openmm/serialization/AmoebaHarmonicBondForceProxy.h"
#include "openmm/serialization/AmoebaHarmonicAngleForceProxy.h"
#include "openmm/serialization/AmoebaHarmonicInPlaneAngleForceProxy.h"
#include "openmm/serialization/AmoebaMultipoleForceProxy.h"
#include "openmm/serialization/AmoebaOutOfPlaneBendForceProxy.h"
#include "openmm/serialization/AmoebaPiTorsionForceProxy.h"
#include "openmm/serialization/AmoebaStretchBendForceProxy.h"
#include "openmm/serialization/AmoebaTorsionForceProxy.h"
#include "openmm/serialization/AmoebaTorsionTorsionForceProxy.h"
#include "openmm/serialization/AmoebaUreyBradleyForceProxy.h"
#include "openmm/serialization/AmoebaVdwForceProxy.h"
#include "openmm/serialization/AmoebaWcaDispersionForceProxy.h"

#if defined(WIN32)
    #include <windows.h>
    extern "C" void registerAmoebaSerializationProxies();
    BOOL WINAPI DllMain(HANDLE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved) {
        if (ul_reason_for_call == DLL_PROCESS_ATTACH)
            registerAmoebaSerializationProxies();
        return TRUE;
    }
#else
    extern "C" void __attribute__((constructor)) registerAmoebaSerializationProxies();
#endif

using namespace OpenMM;

extern "C" void registerAmoebaSerializationProxies() {

    // force libOpenMMSerialization to load before libOpenMMSerialization
    // required in order guarantee initialization of
    //    map<const string, const SerializationProxy*> SerializationProxy::proxiesByType;
    //    map<const string, const SerializationProxy*> SerializationProxy::proxiesByName;

    // is there a cleaner solution?

    std::string file = "libOpenMMSerialization.so";
#ifdef WIN32
    // Tell Windows not to bother the user with ugly error boxes.
    const UINT oldErrorMode = SetErrorMode(SEM_FAILCRITICALERRORS);
    HMODULE handle = LoadLibrary(file.c_str());
    SetErrorMode(oldErrorMode); // Restore previous error mode.
    if (handle == NULL) {
        std::string message;
        std::stringstream(message) << "Error loading library " << file << ": " << GetLastError();
        throw OpenMMException(message);
    }   
#else
    void *handle = dlopen(file.c_str(), RTLD_LAZY | RTLD_GLOBAL);
    if (handle == NULL)
        throw OpenMMException("Error loading library "+file+": "+dlerror());

#endif

    SerializationProxy::registerProxy(typeid(AmoebaGeneralizedKirkwoodForce),         new AmoebaGeneralizedKirkwoodForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaHarmonicBondForce),                new AmoebaHarmonicBondForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaHarmonicAngleForce),               new AmoebaHarmonicAngleForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaHarmonicInPlaneAngleForce),        new AmoebaHarmonicInPlaneAngleForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaMultipoleForce),                   new AmoebaMultipoleForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaOutOfPlaneBendForce),              new AmoebaOutOfPlaneBendForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaPiTorsionForce),                   new AmoebaPiTorsionForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaStretchBendForce),                 new AmoebaStretchBendForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaTorsionForce),                     new AmoebaTorsionForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaTorsionTorsionForce),              new AmoebaTorsionTorsionForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaUreyBradleyForce),                 new AmoebaUreyBradleyForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaVdwForce),                         new AmoebaVdwForceProxy());
    SerializationProxy::registerProxy(typeid(AmoebaWcaDispersionForce),               new AmoebaWcaDispersionForceProxy());
}

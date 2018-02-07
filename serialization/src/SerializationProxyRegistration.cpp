/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2016 Stanford University and the Authors.      *
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
#include "openmm/BrownianIntegrator.h"
#include "openmm/CMAPTorsionForce.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/CompoundIntegrator.h"
#include "openmm/CustomAngleForce.h"
#include "openmm/CustomBondForce.h"
#include "openmm/CustomCentroidBondForce.h"
#include "openmm/CustomCompoundBondForce.h"
#include "openmm/CustomCVForce.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/CustomGBForce.h"
#include "openmm/CustomHbondForce.h"
#include "openmm/CustomIntegrator.h"
#include "openmm/CustomManyParticleForce.h"
#include "openmm/CustomNonbondedForce.h"
#include "openmm/CustomTorsionForce.h"
#include "openmm/GayBerneForce.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/MonteCarloAnisotropicBarostat.h"
#include "openmm/MonteCarloBarostat.h"
#include "openmm/MonteCarloMembraneBarostat.h"
#include "openmm/NonbondedForce.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/RBTorsionForce.h"
#include "openmm/RMSDForce.h"
#include "openmm/System.h"
#include "openmm/TabulatedFunction.h"
#include "openmm/VariableLangevinIntegrator.h"
#include "openmm/VariableVerletIntegrator.h"
#include "openmm/VerletIntegrator.h"

#include "openmm/serialization/SerializationProxy.h"
#include "openmm/serialization/BrownianIntegratorProxy.h"
#include "openmm/serialization/AndersenThermostatProxy.h"
#include "openmm/serialization/CMAPTorsionForceProxy.h"
#include "openmm/serialization/CMMotionRemoverProxy.h"
#include "openmm/serialization/CompoundIntegratorProxy.h"
#include "openmm/serialization/CustomAngleForceProxy.h"
#include "openmm/serialization/CustomBondForceProxy.h"
#include "openmm/serialization/CustomCentroidBondForceProxy.h"
#include "openmm/serialization/CustomCompoundBondForceProxy.h"
#include "openmm/serialization/CustomCVForceProxy.h"
#include "openmm/serialization/CustomExternalForceProxy.h"
#include "openmm/serialization/CustomGBForceProxy.h"
#include "openmm/serialization/CustomHbondForceProxy.h"
#include "openmm/serialization/CustomIntegratorProxy.h"
#include "openmm/serialization/CustomManyParticleForceProxy.h"
#include "openmm/serialization/CustomNonbondedForceProxy.h"
#include "openmm/serialization/CustomTorsionForceProxy.h"
#include "openmm/serialization/GayBerneForceProxy.h"
#include "openmm/serialization/GBSAOBCForceProxy.h"
#include "openmm/serialization/HarmonicAngleForceProxy.h"
#include "openmm/serialization/HarmonicBondForceProxy.h"
#include "openmm/serialization/LangevinIntegratorProxy.h"
#include "openmm/serialization/MonteCarloAnisotropicBarostatProxy.h"
#include "openmm/serialization/MonteCarloBarostatProxy.h"
#include "openmm/serialization/MonteCarloMembraneBarostatProxy.h"
#include "openmm/serialization/NonbondedForceProxy.h"
#include "openmm/serialization/PeriodicTorsionForceProxy.h"
#include "openmm/serialization/RBTorsionForceProxy.h"
#include "openmm/serialization/RMSDForceProxy.h"
#include "openmm/serialization/StateProxy.h"
#include "openmm/serialization/SystemProxy.h"
#include "openmm/serialization/TabulatedFunctionProxies.h"
#include "openmm/serialization/VariableLangevinIntegratorProxy.h"
#include "openmm/serialization/VariableVerletIntegratorProxy.h"
#include "openmm/serialization/VerletIntegratorProxy.h"

#if defined(WIN32)
    #include <windows.h>
    extern "C" void registerSerializationProxies();
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
    SerializationProxy::registerProxy(typeid(BrownianIntegrator), new BrownianIntegratorProxy());
    SerializationProxy::registerProxy(typeid(CMAPTorsionForce), new CMAPTorsionForceProxy());
    SerializationProxy::registerProxy(typeid(CMMotionRemover), new CMMotionRemoverProxy());
    SerializationProxy::registerProxy(typeid(CompoundIntegrator), new CompoundIntegratorProxy());
    SerializationProxy::registerProxy(typeid(Continuous1DFunction), new Continuous1DFunctionProxy());
    SerializationProxy::registerProxy(typeid(Continuous2DFunction), new Continuous2DFunctionProxy());
    SerializationProxy::registerProxy(typeid(Continuous3DFunction), new Continuous3DFunctionProxy());
    SerializationProxy::registerProxy(typeid(CustomAngleForce), new CustomAngleForceProxy());
    SerializationProxy::registerProxy(typeid(CustomBondForce), new CustomBondForceProxy());
    SerializationProxy::registerProxy(typeid(CustomCentroidBondForce), new CustomCentroidBondForceProxy());
    SerializationProxy::registerProxy(typeid(CustomCompoundBondForce), new CustomCompoundBondForceProxy());
    SerializationProxy::registerProxy(typeid(CustomCVForce), new CustomCVForceProxy());
    SerializationProxy::registerProxy(typeid(CustomExternalForce), new CustomExternalForceProxy());
    SerializationProxy::registerProxy(typeid(CustomGBForce), new CustomGBForceProxy());
    SerializationProxy::registerProxy(typeid(CustomHbondForce), new CustomHbondForceProxy());
    SerializationProxy::registerProxy(typeid(CustomIntegrator), new CustomIntegratorProxy());
    SerializationProxy::registerProxy(typeid(CustomManyParticleForce), new CustomManyParticleForceProxy());
    SerializationProxy::registerProxy(typeid(CustomNonbondedForce), new CustomNonbondedForceProxy());
    SerializationProxy::registerProxy(typeid(CustomTorsionForce), new CustomTorsionForceProxy());
    SerializationProxy::registerProxy(typeid(Discrete1DFunction), new Discrete1DFunctionProxy());
    SerializationProxy::registerProxy(typeid(Discrete2DFunction), new Discrete2DFunctionProxy());
    SerializationProxy::registerProxy(typeid(Discrete3DFunction), new Discrete3DFunctionProxy());
    SerializationProxy::registerProxy(typeid(GayBerneForce), new GayBerneForceProxy());
    SerializationProxy::registerProxy(typeid(GBSAOBCForce), new GBSAOBCForceProxy());
    SerializationProxy::registerProxy(typeid(HarmonicAngleForce), new HarmonicAngleForceProxy());
    SerializationProxy::registerProxy(typeid(HarmonicBondForce), new HarmonicBondForceProxy());
    SerializationProxy::registerProxy(typeid(LangevinIntegrator), new LangevinIntegratorProxy());
    SerializationProxy::registerProxy(typeid(MonteCarloAnisotropicBarostat), new MonteCarloAnisotropicBarostatProxy());
    SerializationProxy::registerProxy(typeid(MonteCarloBarostat), new MonteCarloBarostatProxy());
    SerializationProxy::registerProxy(typeid(MonteCarloMembraneBarostat), new MonteCarloMembraneBarostatProxy());
    SerializationProxy::registerProxy(typeid(NonbondedForce), new NonbondedForceProxy());
    SerializationProxy::registerProxy(typeid(PeriodicTorsionForce), new PeriodicTorsionForceProxy());
    SerializationProxy::registerProxy(typeid(RBTorsionForce), new RBTorsionForceProxy());
    SerializationProxy::registerProxy(typeid(RMSDForce), new RMSDForceProxy());
    SerializationProxy::registerProxy(typeid(System), new SystemProxy());
    SerializationProxy::registerProxy(typeid(State), new StateProxy());
    SerializationProxy::registerProxy(typeid(VariableLangevinIntegrator), new VariableLangevinIntegratorProxy());
    SerializationProxy::registerProxy(typeid(VariableVerletIntegrator), new VariableVerletIntegratorProxy());
    SerializationProxy::registerProxy(typeid(VerletIntegrator), new VerletIntegratorProxy());
}
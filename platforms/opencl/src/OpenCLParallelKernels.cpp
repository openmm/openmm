/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "OpenCLParallelKernels.h"

using namespace OpenMM;
using namespace std;

OpenCLParallelCalcForcesAndEnergyKernel::OpenCLParallelCalcForcesAndEnergyKernel(string name, const Platform& platform, OpenCLPlatform::PlatformData& data) :
        CalcForcesAndEnergyKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcForcesAndEnergyKernel(name, platform, *data.contexts[i])));
}

void OpenCLParallelCalcForcesAndEnergyKernel::initialize(const System& system) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system);
}

void OpenCLParallelCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForce, bool includeEnergy) {
    // Copy coordinates over to each device.
    
    OpenCLContext& mainContext = *data.contexts[0];
    mainContext.getPosq().download();
    for (int i = 1; i < (int) data.contexts.size(); i++)
        data.contexts[i]->getPosq().upload(mainContext.getPosq().getHostBuffer());
    
    // Execute the kernel on each device.
    
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).beginComputation(context, includeForce, includeEnergy);
}

double OpenCLParallelCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForce, bool includeEnergy) {
    double energy = 0.0;
    for (int i = 0; i < (int) kernels.size(); i++)
        energy += getKernel(i).finishComputation(context, includeForce, includeEnergy);
    if (includeForce) {
        // Sum the forces from all devices.
        
        for (int i = 0; i < (int) data.contexts.size(); i++)
            data.contexts[i]->getForce().download();
        OpenCLArray<mm_float4>& forces = data.contexts[0]->getForce();
        for (int i = 1; i < (int) data.contexts.size(); i++) {
            OpenCLArray<mm_float4>& contextForces = data.contexts[i]->getForce();
            for (int j = 0; j < forces.getSize(); j++) {
                mm_float4& f1 = forces[j];
                const mm_float4& f2 = contextForces[j];
                f1.x += f2.x;
                f1.y += f2.y;
                f1.z += f2.z;
            }
        }
        forces.upload();
    }
    return energy;
}

OpenCLParallelCalcHarmonicBondForceKernel::OpenCLParallelCalcHarmonicBondForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) :
        CalcHarmonicBondForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcHarmonicBondForceKernel(name, platform, *data.contexts[i], system)));
}

void OpenCLParallelCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double OpenCLParallelCalcHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).execute(context, includeForces, includeEnergy);
}

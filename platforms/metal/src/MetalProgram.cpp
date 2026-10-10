/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2026 Stanford University and the Authors.           *
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

#include "MetalProgram.h"
#include "MetalKernel.h"

using namespace OpenMM;
using namespace std;

MetalProgram::MetalProgram(MetalContext& context, MTL::Library* library) : context(context), library(library) {
}

ComputeKernel MetalProgram::createKernel(const string& name) {
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    MTL::Function* function = library->newFunction(NS::String::string(name.c_str(), NS::UTF8StringEncoding));
    if (function == nullptr) {
        pool->release();
        throw OpenMMException("Error creating kernel "+name+": no such function");
    }
    NS::Error* error = nullptr;
    MTL::ComputePipelineState* pipeline = context.getDevice().newComputePipelineState(function, &error);
    pool->release();
    if (pipeline == nullptr) {
        string message = "Error creating pipeline state for kernel "+name;
        if (error != nullptr)
            message += error->localizedDescription()->utf8String();
        function->release();
        throw OpenMMException(message);
    }
    return shared_ptr<ComputeKernelImpl>(new MetalKernel(context, pipeline, name));
}
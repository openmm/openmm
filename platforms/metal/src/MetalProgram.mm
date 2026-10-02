/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM CUDA Platform.                                      *
 * Source: platforms/cuda/src/CudaProgram.cpp                                 *
 *                                                                            *
 * Original CUDA Platform code:                                               *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 *                                                                            *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the               *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#include "MetalProgram.h"
#include "MetalContext.h"
#include "MetalKernel.h"
#include "MetalKernelSources.h"
#include "MetalLanguagePolicy.h"
#include "MetalSourceAdapter.h"
#import <Metal/Metal.h>
#include <mutex>

using namespace OpenMM;
using namespace std;

struct MetalProgram::Impl {
    MetalContext& context;
    id<MTLLibrary> libraries[2];
    string source;
    map<string, string> defines;
    bool commonSource, strictMath;
    mutex libraryMutex;

    Impl(MetalContext& context, const string& source, const map<string, string>& defines, bool commonSource, bool strictMath) :
            context(context), source(source), defines(defines), commonSource(commonSource), strictMath(strictMath) {
        libraries[0] = nil;
        libraries[1] = nil;
        this->defines.erase("OPENMM_METAL_FLOAT_ACCUMULATORS");
    }

    /** Compile the other ABI only when a kernel first needs it. */
    id<MTLLibrary> getLibrary(bool floating) {
        lock_guard<mutex> lock(libraryMutex);
        int index = commonSource && floating ? 1 : 0;
        if (libraries[index] == nil) {
            if (!commonSource || source.empty())
                throw OpenMMException("Metal program has no Common source for an alternate accumulator ABI");
            const bool fastMath = OPENMM_METAL_FAST_MATH && !floating && !strictMath;
            map<string, string> variantDefines = defines;
            variantDefines["OPENMM_METAL_USE_FAST_MATH"] = fastMath ? "1" : "0";
            string code;
            for (const auto& define : variantDefines)
                code += "#define "+define.first+" "+define.second+"\n";
            if (floating)
                code += "#define OPENMM_METAL_FLOAT_ACCUMULATORS 1\n";
            code += MetalKernelSources::mathPolicy;
            code += MetalKernelSources::common+MetalKernelSources::gbsaTransport+MetalSourceAdapter::translate(source, floating);
            MTLCompileOptions* options = [[MTLCompileOptions alloc] init];
            options.languageVersion = MetalLanguagePolicy::languageVersion(context.getMetalLanguageVersion());
            if (@available(macOS 15.0, *)) {
                options.mathFloatingPointFunctions = fastMath ? MTLMathFloatingPointFunctionsFast : MTLMathFloatingPointFunctionsPrecise;
            }
            else {
                // macOS 13/14 runtimes do not expose the replacement properties.
                options.fastMathEnabled = fastMath;
            }
            NSError* error = nil;
            id<MTLDevice> device = (__bridge id<MTLDevice>) context.getDevice();
            libraries[index] = [device newLibraryWithSource:[NSString stringWithUTF8String:code.c_str()]
                    options:options error:&error];
            if (libraries[index] == nil)
                throw OpenMMException("Error compiling Metal accumulator variant: "+
                        (error == nil ? string("unknown error") : string(error.localizedDescription.UTF8String)));
        }
        return libraries[index];
    }
};

MetalProgram::MetalProgram(MetalContext& context, void* library, bool commonSource, const string& source,
        const map<string, string>& defines, bool strictMath) :
        impl(new Impl(context, source, defines, commonSource, strictMath)), context(context), commonSource(commonSource) {
    int index = commonSource && context.getUseFloatingPointAccumulators() ? 1 : 0;
    impl->libraries[index] = (__bridge id<MTLLibrary>) library;
}

MetalProgram::~MetalProgram() {
}

ComputeKernel MetalProgram::createKernel(const string& name) {
    @autoreleasepool {
        // This shared capture retains both strong library slots even if the
        // caller releases the ComputeProgram before its kernels.
        shared_ptr<Impl> libraries = impl;
        auto lookup = [libraries](bool floating) -> void* {
            return (__bridge void*) libraries->getLibrary(floating);
        };
        int pipelineMaximum = 0;
#if OPENMM_METAL_TUNE_FORCE_PIPELINE_MAX_THREADS
        // Names alone are not enough: synthetic kernels may reuse an entry
        // point with unrelated launch geometry.  The tiled production programs
        // declare their shared force geometry with this numeric source define.
        auto workgroup = impl->defines.find("FORCE_WORK_GROUP_SIZE");
        if (commonSource && workgroup != impl->defines.end() &&
                (name == "computeNonbonded" || name == "computeBornSum" || name == "computeGBSAForce1")) {
            const string& value = workgroup->second;
            if (value != "64" && value != "128" && value != "256")
                throw OpenMMException("Unsupported Metal tiled force workgroup size for pipeline tuning");
            pipelineMaximum = OPENMM_METAL_FORCE_PIPELINE_MAX_THREADS;
            if (pipelineMaximum < stoi(value))
                throw OpenMMException("Metal force pipeline maximum is smaller than the shader's force workgroup size");
        }
#endif
        return ComputeKernel(new MetalKernel(context, name, commonSource, lookup, pipelineMaximum));
    }
}

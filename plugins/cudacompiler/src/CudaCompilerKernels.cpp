/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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

#include "CudaCompilerKernels.h"
#include "openmm/OpenMMException.h"
#include <sstream>
#include <nvrtc.h>

using namespace OpenMM;
using namespace std;

#define CHECK_RESULT(result, prefix) \
    if (result != NVRTC_SUCCESS) { \
        stringstream m; \
        m<<prefix<<": "<<getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }

static string getErrorString(nvrtcResult result) {
    return nvrtcGetErrorString(result);
}

string CudaRuntimeCompilerKernel::createModule(const string& source, const string& flags, CudaContext& cu) {
    // Split the command line flags into an array of options.
    
    stringstream flagsStream(flags);
    string flag;
    vector<string> splitFlags;
    while (flagsStream >> flag)
        splitFlags.push_back(flag);
    int numOptions = splitFlags.size();
    vector<const char*> options(numOptions);
    for (int i = 0; i < numOptions; i++)
        options[i] = &splitFlags[i][0];
    
    // Compile the program to PTX.
    
    nvrtcProgram program;
    CHECK_RESULT(nvrtcCreateProgram(&program, source.c_str(), "", 0, NULL, NULL), "Error creating program");
    try {
        nvrtcResult result = nvrtcCompileProgram(program, options.size(), &options[0]);
        if (result != NVRTC_SUCCESS) {
            size_t logSize;
            nvrtcGetProgramLogSize(program, &logSize);
            vector<char> log(logSize);
            nvrtcGetProgramLog(program, &log[0]);
            throw OpenMMException("Error compiling program: "+string(&log[0]));
        }
        size_t ptxSize;
        nvrtcGetPTXSize(program, &ptxSize);
        vector<char> ptx(ptxSize);
        nvrtcGetPTX(program, &ptx[0]);
        nvrtcDestroyProgram(&program);
        return string(&ptx[0]);
    }
    catch (...) {
        nvrtcDestroyProgram(&program);
        throw;
    }
}

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2015 Stanford University and the Authors.      *
 * Portions copyright (c) 2021 Advanced Micro Devices, Inc.                   *
 * Authors:                                                                   *
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

#include "HipFFT3D.h"
#include "HipContext.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>

using namespace OpenMM;
using namespace std;

HipFFT3D::HipFFT3D(HipContext& context, int xsize, int ysize, int zsize, bool realToComplex, hipStream_t stream, HipArray& in, HipArray& out) :
        context(context), stream(stream) {

    deviceIndex = context.getDeviceIndex();
    inputBuffer = in.getDevicePointer();
    outputBuffer = out.getDevicePointer();
    size_t valueSize = context.getUseDoublePrecision() ? sizeof(double) : sizeof(float);
    inputBufferSize = zsize * ysize * xsize * valueSize;
    if (realToComplex) {
        outputBufferSize = (zsize/2 + 1) * ysize * xsize * valueSize * 2;
    }
    else {
        outputBufferSize = zsize * ysize * xsize * valueSize;
    }

    VkFFTConfiguration configuration = {};
    configuration.performR2C = realToComplex;
    configuration.device = &deviceIndex;
    configuration.num_streams = 1;
    configuration.stream = &this->stream;
    configuration.doublePrecision = context.getUseDoublePrecision();

    configuration.FFTdim = 3;
    configuration.size[0] = zsize;
    configuration.size[1] = ysize;
    configuration.size[2] = xsize;

    configuration.inverseReturnToInputBuffer = true;
    configuration.isInputFormatted = true;
    configuration.inputBufferSize = &inputBufferSize;
    configuration.inputBuffer = &inputBuffer;
    configuration.inputBufferStride[0] = zsize;
    configuration.inputBufferStride[1] = configuration.inputBufferStride[0] * ysize;
    configuration.inputBufferStride[2] = configuration.inputBufferStride[1] * xsize;

    configuration.bufferSize = &outputBufferSize;
    configuration.buffer = &outputBuffer;
    configuration.bufferStride[0] = realToComplex ? (zsize/2 + 1) : zsize;
    configuration.bufferStride[1] = configuration.bufferStride[0] * ysize;
    configuration.bufferStride[2] = configuration.bufferStride[1] * xsize;

    // Combine all parameters into a unique key
    stringstream info;
    int runtimeVersion;
    (void)hipRuntimeGetVersion(&runtimeVersion);
    info << runtimeVersion;
    info << " " << VkFFTGetVersion();
    info << " " << xsize << " " << ysize << " " << zsize;
    info << " " << realToComplex << " " << context.getUseDoublePrecision();

    string cacheFile = context.getCacheFileName(info.str());

    bool hasCache = false;
    vector<char> cacheContent;

    ifstream cache(cacheFile.c_str(), ios::in | ios::binary);
    if (cache.is_open()) {
        cacheContent.insert(cacheContent.begin(), istreambuf_iterator<char>(cache), istreambuf_iterator<char>());
        cache.close();
        hasCache = true;
        // There is an existing cache, load VkFFT kernels from it
        configuration.loadApplicationFromString = 1;
        configuration.loadApplicationString = cacheContent.data();
    }
    else {
        // There is no existing cache, request saving
        configuration.saveApplicationToString = 1;
    }

    app = new VkFFTApplication();
    VkFFTResult fftResult = initializeVkFFT(app, configuration);
    if (fftResult != VKFFT_SUCCESS) {
        throw OpenMMException("Error executing initializeVkFFT: "+context.intToString(fftResult));
    }

    if (!hasCache) {
        // There is no existing cache, create it
        string outputFile = context.getTempFileName() + ".vkfftcache";
        try {
            ofstream out(outputFile.c_str(), ios::out | ios::binary);
            out.write(reinterpret_cast<char*>(app->saveApplicationString), size_t(app->applicationStringSize));
            out.close();
            if (!out.fail()) {
                if (rename(outputFile.c_str(), cacheFile.c_str()) != 0)
                    remove(outputFile.c_str());
            }
        }
        catch (...) {
            // An error occurred.  Possibly we don't have permission to write to the temp directory.
        }
    }
}

HipFFT3D::~HipFFT3D() {
    deleteVkFFT(app);
    delete app;
}

void HipFFT3D::execFFT(bool forward) {
    VkFFTResult fftResult = VkFFTAppend(app, forward ? -1 : 1, NULL);
    if (fftResult != VKFFT_SUCCESS) {
        throw OpenMMException("Error executing VkFFTAppend: "+context.intToString(fftResult));
    }
}

int HipFFT3D::findLegalDimension(int minimum) {
    if (minimum < 1)
        return 1;
    while (true) {
        // Attempt to factor the current value.

        int unfactored = minimum;
        // VkFFT supports prime factors up to 13
        for (int factor = 2; factor <= 13; factor++) {
            while (unfactored > 1 && unfactored%factor == 0)
                unfactored /= factor;
        }
        if (unfactored == 1)
            return minimum;
        minimum++;
    }
}

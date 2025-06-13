/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2025 Stanford University and the Authors.      *
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

#include "CudaFFT3D.h"
#include "CudaContext.h"

using namespace OpenMM;

CudaFFT3D::CudaFFT3D(CudaContext& context, int xsize, int ysize, int zsize, bool realToComplex) :
        context(context), realToComplex(realToComplex), hasInitialized(false) {
    cufftType type1, type2;
    if (realToComplex) {
        if (context.getUseDoublePrecision()) {
            type1 = CUFFT_D2Z;
            type2 = CUFFT_Z2D;
        }
        else {
            type1 = CUFFT_R2C;
            type2 = CUFFT_C2R;
        }
    }
    else {
        if (context.getUseDoublePrecision())
            type1 = type2 = CUFFT_Z2Z;
        else
            type1 = type2 = CUFFT_C2C;
    }
    cufftResult result = cufftPlan3d(&fftForward, xsize, ysize, zsize, type1);
    if (result != CUFFT_SUCCESS)
        throw OpenMMException("Error initializing FFT: "+context.intToString(result));
    result = cufftPlan3d(&fftBackward, xsize, ysize, zsize, type2);
    if (result != CUFFT_SUCCESS)
        throw OpenMMException("Error initializing FFT: "+context.intToString(result));
    hasInitialized = true;
}

CudaFFT3D::~CudaFFT3D() {
    if (hasInitialized) {
        cufftDestroy(fftForward);
        cufftDestroy(fftBackward);
    }
}

void CudaFFT3D::execFFT(ArrayInterface& in, ArrayInterface& out, bool forward) {
    CUdeviceptr in2 = context.unwrap(in).getDevicePointer();
    CUdeviceptr out2 = context.unwrap(out).getDevicePointer();
    cufftResult result;
    if (forward) {
        cufftSetStream(fftForward, context.getCurrentStream());
        if (realToComplex) {
            if (context.getUseDoublePrecision())
                result = cufftExecD2Z(fftForward, (double*) in2, (double2*) out2);
            else
                result = cufftExecR2C(fftForward, (float*) in2, (float2*) out2);
        }
        else {
            if (context.getUseDoublePrecision())
                result = cufftExecZ2Z(fftForward, (double2*) in2, (double2*) out2, CUFFT_FORWARD);
            else
                result = cufftExecC2C(fftForward, (float2*) in2, (float2*) out2, CUFFT_FORWARD);
        }
    }
    else {
        cufftSetStream(fftBackward, context.getCurrentStream());
        if (realToComplex) {
            if (context.getUseDoublePrecision())
                result = cufftExecZ2D(fftBackward, (double2*) in2, (double*) out2);
            else
                result = cufftExecC2R(fftBackward, (float2*) in2, (float*) out2);
        }
        else {
            if (context.getUseDoublePrecision())
                result = cufftExecZ2Z(fftBackward, (double2*) in2, (double2*) out2, CUFFT_INVERSE);
            else
                result = cufftExecC2C(fftBackward, (float2*) in2, (float2*) out2, CUFFT_INVERSE);
        }
    }
    if (result != CUFFT_SUCCESS)
        throw OpenMMException("Error executing FFT: "+context.intToString(result));
}

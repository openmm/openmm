/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2020 Stanford University and the Authors.      *
 * Portions copyright (c) 2021 Advanced Micro Devices, Inc.                   *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "AmoebaHipKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/AmoebaGeneralizedKirkwoodForceImpl.h"
#include "openmm/internal/AmoebaMultipoleForceImpl.h"
#include "openmm/internal/AmoebaWcaDispersionForceImpl.h"
#include "openmm/internal/AmoebaTorsionTorsionForceImpl.h"
#include "openmm/internal/AmoebaVdwForceImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "HipBondedUtilities.h"
#include "HipFFT3D.h"
#include "HipForceInfo.h"
#include "HipKernelSources.h"
#include "SimTKOpenMMRealType.h"
#include "jama_lu.h"

#include <algorithm>
#include <cmath>
#ifdef _MSC_VER
#include <windows.h>
#endif

using namespace OpenMM;
using namespace std;

static void setPeriodicBoxArgs(ComputeContext& cc, ComputeKernel kernel, int index) {
    Vec3 a, b, c;
    cc.getPeriodicBoxVectors(a, b, c);
    if (cc.getUseDoublePrecision()) {
        kernel->setArg(index++, mm_double4(a[0], b[1], c[2], 0.0));
        kernel->setArg(index++, mm_double4(1.0/a[0], 1.0/b[1], 1.0/c[2], 0.0));
        kernel->setArg(index++, mm_double4(a[0], a[1], a[2], 0.0));
        kernel->setArg(index++, mm_double4(b[0], b[1], b[2], 0.0));
        kernel->setArg(index, mm_double4(c[0], c[1], c[2], 0.0));
    }
    else {
        kernel->setArg(index++, mm_float4((float) a[0], (float) b[1], (float) c[2], 0.0f));
        kernel->setArg(index++, mm_float4(1.0f/(float) a[0], 1.0f/(float) b[1], 1.0f/(float) c[2], 0.0f));
        kernel->setArg(index++, mm_float4((float) a[0], (float) a[1], (float) a[2], 0.0f));
        kernel->setArg(index++, mm_float4((float) b[0], (float) b[1], (float) b[2], 0.0f));
        kernel->setArg(index, mm_float4((float) c[0], (float) c[1], (float) c[2], 0.0f));
    }
}

/* -------------------------------------------------------------------------- *
 *                             AmoebaMultipole                                *
 * -------------------------------------------------------------------------- */

HipCalcAmoebaMultipoleForceKernel::~HipCalcAmoebaMultipoleForceKernel() {
    cc.setAsCurrent();
    if (hasInitializedFFT)
        hipfftDestroy(fft);
}

void HipCalcAmoebaMultipoleForceKernel::initialize(const System& system, const AmoebaMultipoleForce& force) {
    CommonCalcAmoebaMultipoleForceKernel::initialize(system, force);
    if (usePME) {
        hipfftResult result = hipfftPlan3d(&fft, gridSizeX, gridSizeY, gridSizeZ, cc.getUseDoublePrecision() ? HIPFFT_Z2Z : HIPFFT_C2C);
        if (result != HIPFFT_SUCCESS)
            throw OpenMMException("Error initializing FFT: "+cc.intToString(result));
        hasInitializedFFT = true;
    }
}

void HipCalcAmoebaMultipoleForceKernel::computeFFT(bool forward) {
    HipArray& grid1 = dynamic_cast<HipContext&>(cc).unwrap(pmeGrid1);
    HipArray& grid2 = dynamic_cast<HipContext&>(cc).unwrap(pmeGrid2);
    if (forward) {
        if (cc.getUseDoublePrecision()) {
            hipfftResult result = hipfftExecZ2Z(fft, (double2*) grid1.getDevicePointer(), (double2*) grid2.getDevicePointer(), HIPFFT_FORWARD);
            if (result != HIPFFT_SUCCESS)
                throw OpenMMException("Error executing FFT: "+cc.intToString(result));
        } else {
            hipfftResult result = hipfftExecC2C(fft, (float2*) grid1.getDevicePointer(), (float2*) grid2.getDevicePointer(), HIPFFT_FORWARD);
            if (result != HIPFFT_SUCCESS)
                throw OpenMMException("Error executing FFT: "+cc.intToString(result));
        }
    }
    else {
        if (cc.getUseDoublePrecision()) {
            hipfftResult result = hipfftExecZ2Z(fft, (double2*) grid2.getDevicePointer(), (double2*) grid1.getDevicePointer(), HIPFFT_BACKWARD);
            if (result != HIPFFT_SUCCESS)
                throw OpenMMException("Error executing FFT: "+cc.intToString(result));
        } else {
            hipfftResult result = hipfftExecC2C(fft, (float2*) grid2.getDevicePointer(), (float2*) grid1.getDevicePointer(), HIPFFT_BACKWARD);
            if (result != HIPFFT_SUCCESS)
                throw OpenMMException("Error executing FFT: "+cc.intToString(result));
        }
    }
}

/* -------------------------------------------------------------------------- *
 *                           HippoNonbondedForce                              *
 * -------------------------------------------------------------------------- */

HipCalcHippoNonbondedForceKernel::~HipCalcHippoNonbondedForceKernel() {
    cc.setAsCurrent();
    if (sort != NULL)
        delete sort;
    if (hasInitializedFFT) {
        hipfftDestroy(fftForward);
        hipfftDestroy(fftBackward);
        hipfftDestroy(dfftForward);
        hipfftDestroy(dfftBackward);
    }
}

void HipCalcHippoNonbondedForceKernel::initialize(const System& system, const HippoNonbondedForce& force) {
    CommonCalcHippoNonbondedForceKernel::initialize(system, force);
    if (usePME) {
        sort = new HipSort(cu, new SortTrait(), cc.getNumAtoms());
        hipfftResult result = hipfftPlan3d(&fftForward, gridSizeX, gridSizeY, gridSizeZ, cc.getUseDoublePrecision() ? HIPFFT_D2Z : HIPFFT_R2C);
        if (result != HIPFFT_SUCCESS)
            throw OpenMMException("Error initializing FFT: "+cc.intToString(result));
        result = hipfftPlan3d(&fftBackward, gridSizeX, gridSizeY, gridSizeZ, cc.getUseDoublePrecision() ? HIPFFT_Z2D : HIPFFT_C2R);
        if (result != HIPFFT_SUCCESS)
            throw OpenMMException("Error initializing FFT: "+cc.intToString(result));
        result = hipfftPlan3d(&dfftForward, dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ, cc.getUseDoublePrecision() ? HIPFFT_D2Z : HIPFFT_R2C);
        if (result != HIPFFT_SUCCESS)
            throw OpenMMException("Error initializing FFT: "+cc.intToString(result));
        result = hipfftPlan3d(&dfftBackward, dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ, cc.getUseDoublePrecision() ? HIPFFT_Z2D : HIPFFT_C2R);
        if (result != HIPFFT_SUCCESS)
            throw OpenMMException("Error initializing FFT: "+cc.intToString(result));
        hasInitializedFFT = true;
    }
}

void HipCalcHippoNonbondedForceKernel::computeFFT(bool forward, bool dispersion) {
    HipArray& grid1 = dynamic_cast<HipContext&>(cc).unwrap(pmeGrid1);
    HipArray& grid2 = dynamic_cast<HipContext&>(cc).unwrap(pmeGrid2);
    if (forward) {
        hipfftHandle fft = dispersion ? dfftForward : fftForward;
        if (cc.getUseDoublePrecision()) {
            hipfftResult result = hipfftExecD2Z(fft, (double*) grid1.getDevicePointer(), (double2*) grid2.getDevicePointer());
            if (result != HIPFFT_SUCCESS)
                throw OpenMMException("Error executing FFT: "+cc.intToString(result));
        } else {
            hipfftResult result = hipfftExecR2C(fft, (float*) grid1.getDevicePointer(), (float2*) grid2.getDevicePointer());
            if (result != HIPFFT_SUCCESS)
                throw OpenMMException("Error executing FFT: "+cc.intToString(result));
        }
    }
    else {
        hipfftHandle fft = dispersion ? dfftBackward : fftBackward;
        if (cc.getUseDoublePrecision()) {
            hipfftResult result = hipfftExecZ2D(fft, (double2*) grid2.getDevicePointer(), (double*) grid1.getDevicePointer());
            if (result != HIPFFT_SUCCESS)
                throw OpenMMException("Error executing FFT: "+cc.intToString(result));
        } else {
            hipfftResult result = hipfftExecC2R(fft, (float2*) grid2.getDevicePointer(), (float*) grid1.getDevicePointer());
            if (result != HIPFFT_SUCCESS)
                throw OpenMMException("Error executing FFT: "+cc.intToString(result));
        }
    }
}

void HipCalcHippoNonbondedForceKernel::sortGridIndex() {
    sort->sort(dynamic_cast<HipContext&>(cc).unwrap(pmeAtomGridIndex));
}

#ifndef __OPENMM_CUDAFFT3D_H__
#define __OPENMM_CUDAFFT3D_H__

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

#include "openmm/common/FFT3D.h"
#include "openmm/common/ArrayInterface.h"
#include <cuda.h>
#include <cufft.h>

namespace OpenMM {

class CudaContext;

/**
 * This class performs three dimensional Fast Fourier Transforms.  It is implemented
 * using cuFFT.
 *
 * FFTs tend to be most efficient when the size of each dimension is a product of
 * small prime factors.  You can call findLegalFFTDimension() on the ComputeContext
 * to determine the smallest size that satisfies this requirement and is greater
 * than or equal to a specified minimum size.
 *
 * Note that this class performs an unnormalized transform.  That means that if you perform
 * a forward transform followed immediately by an inverse transform, the effect is to
 * multiply every value of the original data set by the total number of data points.
 */

class OPENMM_EXPORT_COMMON CudaFFT3D : public FFT3DImpl {
public:
    /**
     * Create a CudaFFT3D object for performing transforms of a particular size.
     *
     * @param context the context in which to perform calculations
     * @param xsize   the first dimension of the data sets on which FFTs will be performed
     * @param ysize   the second dimension of the data sets on which FFTs will be performed
     * @param zsize   the third dimension of the data sets on which FFTs will be performed
     * @param realToComplex  if true, a real-to-complex transform will be done.  Otherwise, it is complex-to-complex.
     */
    CudaFFT3D(CudaContext& context, int xsize, int ysize, int zsize, bool realToComplex=false);
    ~CudaFFT3D();
    /**
     * Perform a Fourier transform.  The transform cannot be done in-place: the input and output
     * arrays must be different.  Also, the input array is used as workspace, so its contents
     * are destroyed.  This also means that both arrays must be large enough to hold complex values,
     * even when performing a real-to-complex transform.
     *
     * When performing a real-to-complex transform, the output data is of size xsize*ysize*(zsize/2+1)
     * and contains only the non-redundant elements.
     *
     * @param in       the data to transform, ordered such that in[x*ysize*zsize + y*zsize + z] contains element (x, y, z)
     * @param out      on exit, this contains the transformed data
     * @param forward  true to perform a forward transform, false to perform an inverse transform
     */
    void execFFT(ArrayInterface& in, ArrayInterface& out, bool forward = true);
private:
    CudaContext& context;
    cufftHandle fftForward;
    cufftHandle fftBackward;
    bool realToComplex, hasInitialized;
};

} // namespace OpenMM

#endif // __OPENMM_CUDAFFT3D_H__

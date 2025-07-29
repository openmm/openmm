#ifndef __OPENMM_FFT3D_H__
#define __OPENMM_FFT3D_H__

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

#include "openmm/common/ArrayInterface.h"
#include <memory>

namespace OpenMM {

/**
 * This class defines a uniform API for three dimensional Fast Fourier Transforms.
 * Each platform provides its own implementation.  Instances can be created by
 * calling createFFT() on a ComputeContext.
 *
 * FFTs tend to be most efficient when the size of each dimension is a product of
 * small prime factors.  You can call findLegalFFTDimension() on the ComputeContext
 * to determine the smallest size that satisfies this requirement and is greater
 * than or equal to a specified minimum size.
 *
 * Note that this class performs an unnormalized transform.  That means that if you perform
 * a forward transform followed immediately by an inverse transform, the effect is to
 * multiply every value of the original data set by the total number of data points.
 * 
 * Instead of referring to this class directly, it is best to use a FFT3D, which is
 * a typedef for a shared_ptr to a FFT3DImpl.  This allows you to treat it as having
 * value semantics, and frees you from having to manage memory.  
 */

class OPENMM_EXPORT_COMMON FFT3DImpl {
public:
    virtual ~FFT3DImpl() {
    }
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
    virtual void execFFT(ArrayInterface& in, ArrayInterface& out, bool forward=true) = 0;
};

typedef std::shared_ptr<FFT3DImpl> FFT3D;

} // namespace OpenMM

#endif // __OPENMM_FFT3D_H__

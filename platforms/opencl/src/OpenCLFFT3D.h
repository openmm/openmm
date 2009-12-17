#ifndef __OPENMM_OPENCLFFT3D_H__
#define __OPENMM_OPENCLFFT3D_H__

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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

#include "OpenCLArray.h"

namespace OpenMM {

class OpenCLFFT3D {
public:
    OpenCLFFT3D(OpenCLContext& context, int xsize, int ysize, int zsize);
    ~OpenCLFFT3D();
    void execFFT(OpenCLArray<mm_float2>& data, bool forward = true);
private:
    cl::Kernel createKernel(int size);
    int xsize, ysize, zsize;
    OpenCLContext& context;
    cl::Kernel xkernel, ykernel, zkernel;
};

} // namespace OpenMM

#endif // __OPENMM_OPENCLFFT3D_H__
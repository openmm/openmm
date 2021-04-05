/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2021 Stanford University and the Authors.      *
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

#include "AmoebaOpenCLKernels.h"

using namespace OpenMM;
using namespace std;

/* -------------------------------------------------------------------------- *
 *                             AmoebaMultipole                                *
 * -------------------------------------------------------------------------- */

OpenCLCalcAmoebaMultipoleForceKernel::~OpenCLCalcAmoebaMultipoleForceKernel() {
    if (fft != NULL)
        delete fft;
}

void OpenCLCalcAmoebaMultipoleForceKernel::initialize(const System& system, const AmoebaMultipoleForce& force) {
    CommonCalcAmoebaMultipoleForceKernel::initialize(system, force);
    if (usePME)
        fft = new OpenCLFFT3D(dynamic_cast<OpenCLContext&>(cc), gridSizeX, gridSizeY, gridSizeZ, false);
}

void OpenCLCalcAmoebaMultipoleForceKernel::computeForwardFFT() {
    OpenCLArray& grid = dynamic_cast<OpenCLContext&>(cc).unwrap(pmeGrid);
    fft->execFFT(grid, grid, true);
}

void OpenCLCalcAmoebaMultipoleForceKernel::computeInverseFFT() {
    OpenCLArray& grid = dynamic_cast<OpenCLContext&>(cc).unwrap(pmeGrid);
    fft->execFFT(grid, grid, false);
}

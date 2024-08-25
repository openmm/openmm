/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2021 Stanford University and the Authors.      *
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
#include "openmm/common/ContextSelector.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/AmoebaGeneralizedKirkwoodForceImpl.h"
#include "openmm/internal/AmoebaMultipoleForceImpl.h"
#include "openmm/internal/AmoebaWcaDispersionForceImpl.h"
#include "openmm/internal/AmoebaTorsionTorsionForceImpl.h"
#include "openmm/internal/AmoebaVdwForceImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "HipBondedUtilities.h"
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

/* -------------------------------------------------------------------------- *
 *                             AmoebaMultipole                                *
 * -------------------------------------------------------------------------- */

HipCalcAmoebaMultipoleForceKernel::~HipCalcAmoebaMultipoleForceKernel() {
    ContextSelector selector(cc);
    if (fft != NULL)
        delete fft;
}

void HipCalcAmoebaMultipoleForceKernel::initialize(const System& system, const AmoebaMultipoleForce& force) {
    CommonCalcAmoebaMultipoleForceKernel::initialize(system, force);
    if (usePME) {
        ContextSelector selector(cc);
        HipArray& grid1 = cu.unwrap(pmeGrid1);
        HipArray& grid2 = cu.unwrap(pmeGrid2);
        fft = cu.createFFT(gridSizeX, gridSizeY, gridSizeZ, false, cu.getCurrentStream(), grid1, grid2);
    }
}

void HipCalcAmoebaMultipoleForceKernel::computeFFT(bool forward) {
    fft->execFFT(forward);
}

/* -------------------------------------------------------------------------- *
 *                           HippoNonbondedForce                              *
 * -------------------------------------------------------------------------- */

HipCalcHippoNonbondedForceKernel::~HipCalcHippoNonbondedForceKernel() {
    ContextSelector selector(cc);
    if (sort != NULL)
        delete sort;
    if (fft != NULL)
        delete fft;
    if (dfft != NULL)
        delete dfft;
}

void HipCalcHippoNonbondedForceKernel::initialize(const System& system, const HippoNonbondedForce& force) {
    CommonCalcHippoNonbondedForceKernel::initialize(system, force);
    if (usePME) {
        ContextSelector selector(cc);
        sort = new HipSort(cu, new SortTrait(), cc.getNumAtoms());
        HipArray& grid1 = cu.unwrap(pmeGrid1);
        HipArray& grid2 = cu.unwrap(pmeGrid2);
        fft = cu.createFFT(gridSizeX, gridSizeY, gridSizeZ, true, cu.getCurrentStream(), grid1, grid2);
        dfft = cu.createFFT(dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ, true, cu.getCurrentStream(), grid1, grid2);
    }
}

void HipCalcHippoNonbondedForceKernel::computeFFT(bool forward, bool dispersion) {
    if (dispersion) {
        dfft->execFFT(forward);
    }
    else {
        fft->execFFT(forward);
    }
}

void HipCalcHippoNonbondedForceKernel::sortGridIndex() {
    sort->sort(dynamic_cast<HipContext&>(cc).unwrap(pmeAtomGridIndex));
}

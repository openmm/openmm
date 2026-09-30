/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "openmm/kernels.h"

using namespace OpenMM;

CalcForcesAndEnergyKernel::~CalcForcesAndEnergyKernel() = default;

UpdateStateDataKernel::~UpdateStateDataKernel() = default;

ApplyConstraintsKernel::~ApplyConstraintsKernel() = default;

VirtualSitesKernel::~VirtualSitesKernel() = default;

MinimizeKernel::~MinimizeKernel() = default;

CalcHarmonicBondForceKernel::~CalcHarmonicBondForceKernel() = default;

CalcCustomBondForceKernel::~CalcCustomBondForceKernel() = default;

CalcHarmonicAngleForceKernel::~CalcHarmonicAngleForceKernel() = default;

CalcCustomAngleForceKernel::~CalcCustomAngleForceKernel() = default;

CalcPeriodicTorsionForceKernel::~CalcPeriodicTorsionForceKernel() = default;

CalcRBTorsionForceKernel::~CalcRBTorsionForceKernel() = default;

CalcCMAPTorsionForceKernel::~CalcCMAPTorsionForceKernel() = default;

CalcCustomTorsionForceKernel::~CalcCustomTorsionForceKernel() = default;

CalcNonbondedForceKernel::~CalcNonbondedForceKernel() = default;

CalcConstantPotentialForceKernel::~CalcConstantPotentialForceKernel() = default;

CalcCustomNonbondedForceKernel::~CalcCustomNonbondedForceKernel() = default;

CalcGBSAOBCForceKernel::~CalcGBSAOBCForceKernel() = default;

CalcCustomGBForceKernel::~CalcCustomGBForceKernel() = default;

CalcCustomExternalForceKernel::~CalcCustomExternalForceKernel() = default;

CalcCustomHbondForceKernel::~CalcCustomHbondForceKernel() = default;

CalcCustomCentroidBondForceKernel::~CalcCustomCentroidBondForceKernel() = default;

CalcCustomCompoundBondForceKernel::~CalcCustomCompoundBondForceKernel() = default;

CalcCustomManyParticleForceKernel::~CalcCustomManyParticleForceKernel() = default;

CalcGayBerneForceKernel::~CalcGayBerneForceKernel() = default;

CalcLCPOForceKernel::~CalcLCPOForceKernel() = default;

CalcCustomCVForceKernel::~CalcCustomCVForceKernel() = default;

CalcRMSDForceKernel::~CalcRMSDForceKernel() = default;

CalcRGForceKernel::~CalcRGForceKernel() = default;

CalcOrientationRestraintForceKernel::~CalcOrientationRestraintForceKernel() = default;

IntegrateVerletStepKernel::~IntegrateVerletStepKernel() = default;

IntegrateNoseHooverStepKernel::~IntegrateNoseHooverStepKernel() = default;

IntegrateLangevinMiddleStepKernel::~IntegrateLangevinMiddleStepKernel() = default;

IntegrateBrownianStepKernel::~IntegrateBrownianStepKernel() = default;

IntegrateVariableLangevinStepKernel::~IntegrateVariableLangevinStepKernel() = default;

IntegrateVariableVerletStepKernel::~IntegrateVariableVerletStepKernel() = default;

IntegrateCustomStepKernel::~IntegrateCustomStepKernel() = default;

IntegrateDPDStepKernel::~IntegrateDPDStepKernel() = default;

IntegrateQTBStepKernel::~IntegrateQTBStepKernel() = default;

ApplyAndersenThermostatKernel::~ApplyAndersenThermostatKernel() = default;

ApplyMonteCarloBarostatKernel::~ApplyMonteCarloBarostatKernel() = default;

RemoveCMMotionKernel::~RemoveCMMotionKernel() = default;

CalcPmeReciprocalForceKernel::~CalcPmeReciprocalForceKernel() = default;

CalcDispersionPmeReciprocalForceKernel::~CalcDispersionPmeReciprocalForceKernel() = default;

CalcATMForceKernel::~CalcATMForceKernel() = default;

CalcCustomCPPForceKernel::~CalcCustomCPPForceKernel() = default;

CalcPythonForceKernel::~CalcPythonForceKernel() = default;

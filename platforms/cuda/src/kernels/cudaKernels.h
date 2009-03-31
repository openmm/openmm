/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Scott Le Grand, Peter Eastman                                     *
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

#include "gputypes.h"

// Initialization
extern void kClearForces(gpuContext gpu);
extern void kCalculateObcGbsaBornSum(gpuContext gpu);
extern void kReduceObcGbsaBornSum(gpuContext gpu);
extern void kGenerateRandoms(gpuContext gpu);

// Main loop
extern void kCalculateCDLJObcGbsaForces1(gpuContext gpu);
extern void kCalculateCDLJForces(gpuContext gpu);
extern void kReduceObcGbsaBornForces(gpuContext gpu);
extern void kCalculateObcGbsaForces2(gpuContext gpu);
extern void kCalculateLocalForces(gpuContext gpu);
extern void kCalculateAndersenThermostat(gpuContext gpu);
extern void kReduceBornSumAndForces(gpuContext gpu);
extern void kUpdatePart1(gpuContext gpu);
extern void kApplyFirstShake(gpuContext gpu);
extern void kApplyFirstSettle(gpuContext gpu);
extern void kApplyFirstLincs(gpuContext gpu);
extern void kUpdatePart2(gpuContext gpu);
extern void kApplySecondShake(gpuContext gpu);
extern void kApplySecondSettle(gpuContext gpu);
extern void kApplySecondLincs(gpuContext gpu);
extern void kVerletUpdatePart1(gpuContext gpu);
extern void kVerletUpdatePart2(gpuContext gpu);
extern void kBrownianUpdatePart1(gpuContext gpu);
extern void kBrownianUpdatePart2(gpuContext gpu);

// Extras
extern void kReduceForces(gpuContext gpu);
extern void kClearBornForces(gpuContext gpu);

// Initializers
extern void SetCalculateCDLJObcGbsaForces1Sim(gpuContext gpu);
extern void GetCalculateCDLJObcGbsaForces1Sim(gpuContext gpu);
extern void SetCalculateCDLJForcesSim(gpuContext gpu);
extern void GetCalculateCDLJForcesSim(gpuContext gpu);
extern void SetCalculateLocalForcesSim(gpuContext gpu);
extern void GetCalculateLocalForcesSim(gpuContext gpu);
extern void SetCalculateObcGbsaBornSumSim(gpuContext gpu);
extern void GetCalculateObcGbsaBornSumSim(gpuContext gpu);
extern void SetCalculateObcGbsaForces2Sim(gpuContext gpu);
extern void GetCalculateObcGbsaForces2Sim(gpuContext gpu);
extern void SetCalculateAndersenThermostatSim(gpuContext gpu);
extern void GetCalculateAndersenThermostatSim(gpuContext gpu);
extern void SetForcesSim(gpuContext gpu);
extern void GetForcesSim(gpuContext gpu);
extern void SetUpdateShakeHSim(gpuContext gpu);
extern void GetUpdateShakeHSim(gpuContext gpu);
extern void SetSettleSim(gpuContext gpu);
extern void GetSettleSim(gpuContext gpu);
extern void SetLincsSim(gpuContext gpu);
extern void GetLincsSim(gpuContext gpu);
extern void SetVerletUpdateSim(gpuContext gpu);
extern void GetVerletUpdateSim(gpuContext gpu);
extern void SetBrownianUpdateSim(gpuContext gpu);
extern void GetBrownianUpdateSim(gpuContext gpu);
extern void SetRandomSim(gpuContext gpu);
extern void GetRandomSim(gpuContext gpu);

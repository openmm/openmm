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
extern void kApplyFirstCShake(gpuContext gpu);
extern void kApplyFirstSettle(gpuContext gpu);
extern void kApplyFirstLincs(gpuContext gpu);
extern void kUpdatePart2(gpuContext gpu);
extern void kApplySecondShake(gpuContext gpu);
extern void kApplySecondCShake(gpuContext gpu);
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
extern void SetCShakeSim(gpuContext gpu);
extern void GetCShakeSim(gpuContext gpu);
extern void SetLincsSim(gpuContext gpu);
extern void GetLincsSim(gpuContext gpu);
extern void SetVerletUpdateSim(gpuContext gpu);
extern void GetVerletUpdateSim(gpuContext gpu);
extern void SetBrownianUpdateSim(gpuContext gpu);
extern void GetBrownianUpdateSim(gpuContext gpu);
extern void SetRandomSim(gpuContext gpu);
extern void GetRandomSim(gpuContext gpu);

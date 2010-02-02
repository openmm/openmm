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
extern void kClearEnergy(gpuContext gpu);
extern void kClearBornForces(gpuContext gpu);
extern void kClearObcGbsaBornSum(gpuContext gpu);
extern void kCalculateObcGbsaBornSum(gpuContext gpu);
extern void kReduceObcGbsaBornSum(gpuContext gpu);
extern void kCalculateGBVIBornSum(gpuContext gpu);
extern void kReduceGBVIBornSum(gpuContext gpu);
extern void kClearGBVIBornSum( gpuContext gpu );
extern void kGenerateRandoms(gpuContext gpu);

// Main loop
extern void kCalculateCDLJObcGbsaForces1(gpuContext gpu);
extern void kCalculateCDLJGBVIForces1(gpuContext gpu);
extern void kCalculateCDLJForces(gpuContext gpu);
extern void kCalculateCustomBondForces(gpuContext gpu);
extern void kCalculateCustomAngleForces(gpuContext gpu);
extern void kCalculateCustomTorsionForces(gpuContext gpu);
extern void kCalculateCustomExternalForces(gpuContext gpu);
extern void kCalculateCustomNonbondedForces(gpuContext gpu, bool neighborListValid);
extern void kReduceObcGbsaBornForces(gpuContext gpu);
extern void kCalculateObcGbsaForces2(gpuContext gpu);
extern void kCalculateGBVIForces2(gpuContext gpu);
extern void kCalculateLocalForces(gpuContext gpu);
extern void kCalculateAndersenThermostat(gpuContext gpu);
extern void kReduceBornSumAndForces(gpuContext gpu);
extern void kApplyFirstShake(gpuContext gpu);
extern void kApplyFirstCCMA(gpuContext gpu);
extern void kApplyFirstSettle(gpuContext gpu);
extern void kApplySecondShake(gpuContext gpu);
extern void kApplySecondCCMA(gpuContext gpu);
extern void kApplySecondSettle(gpuContext gpu);
extern void kLangevinUpdatePart1(gpuContext gpu);
extern void kLangevinUpdatePart2(gpuContext gpu);
extern void kSelectLangevinStepSize(gpuContext gpu, float maxTimeStep);
extern void kVerletUpdatePart1(gpuContext gpu);
extern void kVerletUpdatePart2(gpuContext gpu);
extern void kSelectVerletStepSize(gpuContext gpu, float maxTimeStep);
extern void kBrownianUpdatePart1(gpuContext gpu);
extern void kBrownianUpdatePart2(gpuContext gpu);

// Extras
extern void kReduceForces(gpuContext gpu);
extern double kReduceEnergy(gpuContext gpu);

// Initializers
extern void SetCalculateCDLJObcGbsaForces1Sim(gpuContext gpu);
extern void GetCalculateCDLJObcGbsaForces1Sim(gpuContext gpu);
extern void SetCalculateCDLJForcesSim(gpuContext gpu);
extern void GetCalculateCDLJForcesSim(gpuContext gpu);
extern void SetCalculateCustomBondForcesSim(gpuContext gpu);
extern void GetCalculateCustomBondForcesSim(gpuContext gpu);
extern void SetCalculateCustomAngleForcesSim(gpuContext gpu);
extern void GetCalculateCustomAngleForcesSim(gpuContext gpu);
extern void SetCalculateCustomTorsionForcesSim(gpuContext gpu);
extern void GetCalculateCustomTorsionForcesSim(gpuContext gpu);
extern void SetCalculateCustomExternalForcesSim(gpuContext gpu);
extern void GetCalculateCustomExternalForcesSim(gpuContext gpu);
extern void SetCalculateCustomNonbondedForcesSim(gpuContext gpu);
extern void GetCalculateCustomNonbondedForcesSim(gpuContext gpu);
extern void SetCalculateLocalForcesSim(gpuContext gpu);
extern void GetCalculateLocalForcesSim(gpuContext gpu);
extern void SetCalculateObcGbsaBornSumSim(gpuContext gpu);
extern void GetCalculateObcGbsaBornSumSim(gpuContext gpu);
extern void SetCalculateGBVIBornSumSim(gpuContext gpu);
extern void GetCalculateGBVIBornSumSim(gpuContext gpu);
extern void SetCalculateObcGbsaForces2Sim(gpuContext gpu);
extern void GetCalculateObcGbsaForces2Sim(gpuContext gpu);
extern void SetCalculateGBVIForces2Sim(gpuContext gpu);
extern void GetCalculateGBVIForces2Sim(gpuContext gpu);
extern void SetCalculateAndersenThermostatSim(gpuContext gpu);
extern void GetCalculateAndersenThermostatSim(gpuContext gpu);
extern void SetCalculatePMESim(gpuContext gpu);
extern void GetCalculatePMESim(gpuContext gpu);
extern void SetForcesSim(gpuContext gpu);
extern void GetForcesSim(gpuContext gpu);
extern void SetShakeHSim(gpuContext gpu);
extern void GetShakeHSim(gpuContext gpu);
extern void SetLangevinUpdateSim(gpuContext gpu);
extern void GetLangevinUpdateSim(gpuContext gpu);
extern void SetSettleSim(gpuContext gpu);
extern void GetSettleSim(gpuContext gpu);
extern void SetCCMASim(gpuContext gpu);
extern void GetCCMASim(gpuContext gpu);
extern void SetVerletUpdateSim(gpuContext gpu);
extern void GetVerletUpdateSim(gpuContext gpu);
extern void SetBrownianUpdateSim(gpuContext gpu);
extern void GetBrownianUpdateSim(gpuContext gpu);
extern void SetRandomSim(gpuContext gpu);
extern void GetRandomSim(gpuContext gpu);
extern void SetCustomBondForceExpression(const Expression<256>& expression);
extern void SetCustomBondEnergyExpression(const Expression<256>& expression);
extern void SetCustomBondGlobalParams(const std::vector<float>&  paramValues);
extern void SetCustomAngleForceExpression(const Expression<256>& expression);
extern void SetCustomAngleEnergyExpression(const Expression<256>& expression);
extern void SetCustomAngleGlobalParams(const std::vector<float>&  paramValues);
extern void SetCustomTorsionForceExpression(const Expression<256>& expression);
extern void SetCustomTorsionEnergyExpression(const Expression<256>& expression);
extern void SetCustomTorsionGlobalParams(const std::vector<float>&  paramValues);
extern void SetCustomExternalForceExpressions(const Expression<256>& expressionX, const Expression<256>& expressionY, const Expression<256>& expressionZ);
extern void SetCustomExternalEnergyExpression(const Expression<256>& expression);
extern void SetCustomExternalGlobalParams(const std::vector<float>& paramValues);
extern void SetCustomNonbondedForceExpression(const Expression<256>& expression);
extern void SetCustomNonbondedEnergyExpression(const Expression<256>& expression);
extern void SetCustomNonbondedGlobalParams(const std::vector<float>& paramValues);

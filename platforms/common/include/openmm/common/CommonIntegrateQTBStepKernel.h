#ifndef OPENMM_COMMONINTEGRATEQTBSTEPKERNEL_H_
#define OPENMM_COMMONINTEGRATEQTBSTEPKERNEL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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

#include "openmm/common/ComputeArray.h"
#include "openmm/common/ComputeContext.h"
#include "openmm/Platform.h"
#include "openmm/kernels.h"

namespace OpenMM {

/**
 * This kernel is invoked by QTBIntegrator to take one time step.
 */
class CommonIntegrateQTBStepKernel : public IntegrateQTBStepKernel {
public:
    CommonIntegrateQTBStepKernel(std::string name, const Platform& platform, ComputeContext& cc) : IntegrateQTBStepKernel(name, platform), cc(cc),
        hasInitializedKernels(false) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the QTBIntegrator this kernel will be used for
     */
    void initialize(const System& system, const QTBIntegrator& integrator);
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the QTBIntegrator this kernel is being used for
     */
    void execute(ContextImpl& context, const QTBIntegrator& integrator);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the QTBIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const QTBIntegrator& integrator);
    /**
     * Get the adapted friction coefficients for a particle.
     * 
     * @param context    the context in which to execute this kernel
     * @param particle   the index of the particle for which to get the friction
     * @param friction   the adapted friction coefficients used in generating the
     *                   random force
     */
    void getAdaptedFriction(ContextImpl& context, int particle, std::vector<double>& friction) const;
    /**
     * Set the adapted friction coefficients for a particle.  This affects the
     * specified particle, and all others that have the same type.
     * 
     * @param context    the context in which to execute this kernel
     * @param particle   the index of the particle for which to get the friction
     * @param friction   the adapted friction coefficients used in generating the
     *                   random force
     */
    void setAdaptedFriction(ContextImpl& context, int particle, const std::vector<double>& friction);
    /**
     * Write the adapted friction to a checkpoint.
     * 
     * @param context    the context in which to execute this kernel
     * @param stream     the stream to write the checkpoint to
     */
    void createCheckpoint(ContextImpl& context, std::ostream& stream) const;
    /**
     * Load the adapted friction from a checkpoint.
     * 
     * @param context    the context in which to execute this kernel
     * @param stream     the stream to read the checkpoint from
     */
    void loadCheckpoint(ContextImpl& context, std::istream& stream);
private:
    std::string createFFT(int size, int inputIndex, int& outputIndex, bool forward);
    ComputeContext& cc;
    int segmentLength, stepIndex, numFreq, numTypes;
    double prevTemp, dt, friction;
    std::vector<int> particleTypeVec;
    ComputeArray oldDelta, noise, randomForce, segmentVelocity, thetad, cutoffFunction, workspace;
    ComputeArray particleType, typeParticleCount, typeAdaptationRate, adaptedFriction, dfdt;
    ComputeKernel kernel1, kernel2, kernel3, noiseKernel, forceKernel, adapt1Kernel, adapt2Kernel;
    bool hasInitializedKernels;
};

} // namespace OpenMM

#endif /*OPENMM_COMMONINTEGRATEQTBSTEPKERNEL_H_*/

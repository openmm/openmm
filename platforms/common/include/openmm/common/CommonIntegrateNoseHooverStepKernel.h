#ifndef OPENMM_COMMONINTEGRATENOSEHOOVERSTEPKERNEL_H_
#define OPENMM_COMMONINTEGRATENOSEHOOVERSTEPKERNEL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
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
#include <map>
#include <string>
#include <vector>

namespace OpenMM {

/*
 * This kernel is invoked by NoseHooverIntegrator to take one time step.
 */
class CommonIntegrateNoseHooverStepKernel : public IntegrateNoseHooverStepKernel {
public:
    CommonIntegrateNoseHooverStepKernel(std::string name, const Platform& platform, ComputeContext& cc) :
                                  IntegrateNoseHooverStepKernel(name, platform), cc(cc), hasInitializedKernels(false),
                                  hasInitializedKineticEnergyKernel(false), hasInitializedHeatBathEnergyKernel(false),
                                  hasInitializedScaleVelocitiesKernel(false), hasInitializedPropagateKernel(false) {}
    ~CommonIntegrateNoseHooverStepKernel() {}
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the NoseHooverIntegrator this kernel will be used for
     */
    void initialize(const System& system, const NoseHooverIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the VerletIntegrator this kernel is being used for
     */
    void execute(ContextImpl& context, const NoseHooverIntegrator& integrator);
    /**
     * Compute the kinetic energy.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the NoseHooverIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const NoseHooverIntegrator& integrator);
    /**
     * Execute the kernel that propagates the Nose Hoover chain and determines the velocity scale factor.
     * 
     * @param context  the context in which to execute this kernel
     * @param noseHooverChain the object describing the chain to be propagated.
     * @param kineticEnergy the {center of mass, relative} kineticEnergies of the particles being thermostated by this chain.
     * @param timeStep the time step used by the integrator.
     * @return the velocity scale factor to apply to the particles associated with this heat bath.
     */
    std::pair<double, double> propagateChain(ContextImpl& context, const NoseHooverChain &noseHooverChain, std::pair<double, double> kineticEnergy, double timeStep);
    /**
     * Execute the kernal that computes the total (kinetic + potential) heat bath energy.
     *
     * @param context the context in which to execute this kernel
     * @param noseHooverChain the chain whose energy is to be determined.
     * @return the total heat bath energy.
     */
    double computeHeatBathEnergy(ContextImpl& context, const NoseHooverChain &noseHooverChain);
    /**
     * Execute the kernel that computes the kinetic energy for a subset of atoms,
     * or the relative kinetic energy of Drude particles with respect to their parent atoms
     *
     * @param context the context in which to execute this kernel
     * @param noseHooverChain the chain whose energy is to be determined.
     * @param downloadValue whether the computed value should be downloaded and returned.
     */
    std::pair<double, double> computeMaskedKineticEnergy(ContextImpl& context, const NoseHooverChain &noseHooverChain, bool downloadValue);
    /**
     * Execute the kernel that scales the velocities of particles associated with a nose hoover chain
     *
     * @param context the context in which to execute this kernel
     * @param noseHooverChain the chain whose energy is to be determined.
     * @param scaleFactor the multiplicative factor by which {absolute, relative} velocities are scaled.
     */
    void scaleVelocities(ContextImpl& context, const NoseHooverChain &noseHooverChain, std::pair<double, double> scaleFactor);
    /**
     * Write the chain states to a checkpoint.
     */
    void createCheckpoint(ContextImpl& context, std::ostream& stream) const;
    /**
     * Load the chain states from a checkpoint.
     */
    void loadCheckpoint(ContextImpl& context, std::istream& stream);
    /**
     * Get the internal states of all chains.
     * 
     * @param context       the context for which to get the states
     * @param positions     element [i][j] contains the position of bead j for chain i
     * @param velocities    element [i][j] contains the velocity of bead j for chain i
     */
    void getChainStates(ContextImpl& context, std::vector<std::vector<double> >& positions, std::vector<std::vector<double> >& velocities) const;
    /**
     * Set the internal states of all chains.
     * 
     * @param context       the context for which to get the states
     * @param positions     element [i][j] contains the position of bead j for chain i
     * @param velocities    element [i][j] contains the velocity of bead j for chain i
     */
    void setChainStates(ContextImpl& context, const std::vector<std::vector<double> >& positions, const std::vector<std::vector<double> >& velocities);
private:
    ComputeContext& cc;
    float prevMaxPairDistance;
    ComputeArray maxPairDistanceBuffer, pairListBuffer, atomListBuffer, pairTemperatureBuffer, oldDelta;
    std::map<int, ComputeArray> chainState;
    ComputeKernel kernel1, kernel2, kernel3, kernel4, kernelHardWall;
    bool hasInitializedKernels;
    ComputeKernel reduceEnergyKernel;
    ComputeKernel computeHeatBathEnergyKernel;
    ComputeKernel computeAtomsKineticEnergyKernel;
    ComputeKernel computePairsKineticEnergyKernel;
    ComputeKernel scaleAtomsVelocitiesKernel;
    ComputeKernel scalePairsVelocitiesKernel;
    ComputeArray energyBuffer, scaleFactorBuffer, kineticEnergyBuffer, chainMasses, chainForces, heatBathEnergy;
    std::map<int, ComputeArray> atomlists, pairlists;
    std::map<int, ComputeKernel> propagateKernels;
    bool hasInitializedPropagateKernel;
    bool hasInitializedKineticEnergyKernel;
    bool hasInitializedHeatBathEnergyKernel;
    bool hasInitializedScaleVelocitiesKernel;
};

} // namespace OpenMM

#endif /*OPENMM_COMMONINTEGRATENOSEHOOVERSTEPKERNEL_H_*/

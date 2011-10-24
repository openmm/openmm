#ifndef OPENMM_FREE_ENERGY_CUDA_KERNELS_H_
#define OPENMM_FREE_ENERGY_CUDA_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "CudaPlatform.h"
#include "openmm/kernels.h"
#include "kernels/gputypes.h"
#include "openmm/System.h"
#include "OpenMMFreeEnergy.h"
#include "openmm/freeEnergyKernels.h"
#include "FreeEnergyCudaData.h"

namespace OpenMM {

/**
 * This kernel is invoked by NonbondedSoftcoreForce to calculate the forces acting on the system.
 */
class CudaFreeEnergyCalcNonbondedSoftcoreForceKernel : public CalcNonbondedSoftcoreForceKernel {

public:
    CudaFreeEnergyCalcNonbondedSoftcoreForceKernel(std::string name, const Platform& platform, FreeEnergyCudaData& data, System& system) :
             CalcNonbondedSoftcoreForceKernel(name, platform), data(data), system(system) {

        numExceptions        = 0;
        numParticles         = 0;
        bIncludeGBSA         = false;
        bIncludeGBVI         = false;
        includeSoftcore      = false;
        log                  = NULL;
        data.incrementKernelCount();
    }

    ~CudaFreeEnergyCalcNonbondedSoftcoreForceKernel();

    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the NonbondedForce this kernel will be used for
     */
    void initialize(const System& system, const NonbondedSoftcoreForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Get flag signalling whether GBSA/OBC force is included
     *
     * @return flag
     */
    bool getIncludeGBSA( void ) const;
    /**
     * Set flag signalling whether GBSA/OBC force is included
     *
     * @param inputIncludeGBSA input flag value
     */
    void setIncludeGBSA( bool inputIncludeGBSA );
    /**
     * Get flag signalling whether GB/VI force is included
     *
     * @return flag
     */
    bool getIncludeGBVI( void ) const;
    /**
     * Set flag signalling whether GB/VI force is included
     *
     * @param inputIncludeGBVI input flag value
     */
    void setIncludeGBVI( bool inputIncludeGBVI );
    /**
     * Get flag signalling whether softcore force is included
     *
     * @return flag
     */
    int getIncludeSoftcore( void ) const;
    /**
     * Set flag signalling whether GB/VI force is included
     *
     * @param inputIncludeGBVI input flag value
     */
    void setIncludeSoftcore( int inputSoftcore);
    /**
     * Get number of exceptions
     *
     * @return number of exceptions
     */
    int getNumExceptions( void ) const;
private:
    FreeEnergyCudaData& data;
    class ForceInfo;
    int numParticles;
    System& system;
    bool bIncludeGBSA;
    bool bIncludeGBVI;
    int includeSoftcore;
    int numExceptions;
    FILE* log;
};

/**
 * This kernel is invoked by GBSAOBCForce to calculate the forces acting on the system.
 */
class CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel : public CalcGBSAOBCSoftcoreForceKernel {
public:
    CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel(std::string name, const Platform& platform, FreeEnergyCudaData& data) :
       CalcGBSAOBCSoftcoreForceKernel(name, platform), data(data) {
        log                  = NULL;
        data.incrementKernelCount();
    }
    ~CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the GBSAOBCForce this kernel will be used for
     */
    void initialize(const System& system, const GBSAOBCSoftcoreForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
private:
    FreeEnergyCudaData& data;
    class ForceInfo;
    FILE* log;
};

/**
 * This kernel is invoked by GBVIForce to calculate the forces acting on the system.
 */
class CudaFreeEnergyCalcGBVISoftcoreForceKernel : public CalcGBVISoftcoreForceKernel {
public:
    CudaFreeEnergyCalcGBVISoftcoreForceKernel(std::string name, const Platform& platform, FreeEnergyCudaData& data) :
         CalcGBVISoftcoreForceKernel(name, platform), data(data) {

        log                  = NULL;
        quinticScaling       = 0;
        data.incrementKernelCount();
    }

    ~CudaFreeEnergyCalcGBVISoftcoreForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system      the System this kernel will be applied to
     * @param force       the GBVIForce this kernel will be used for
     * @param scaledRadii the scaled radii (Eq. 5 of Labute paper)
     */
    void initialize(const System& system, const GBVISoftcoreForce& force, const std::vector<double> & scaledRadii);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);

private:
    FreeEnergyCudaData& data;
    class ForceInfo;
    FILE* log;
    int quinticScaling;
};

} // namespace OpenMM

#endif /*OPENMM_FREE_ENERGY_CUDA_KERNELS_H_*/

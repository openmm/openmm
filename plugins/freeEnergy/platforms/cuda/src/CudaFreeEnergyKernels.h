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
#include "kernels/GpuNonbondedSoftcore.h"
#include "kernels/GpuLJ14Softcore.h"
#include "kernels/GpuObcGbsaSoftcore.h"
#include "kernels/GpuGBVISoftcore.h"

//#define FreeEnergyDebug

namespace OpenMM {

/**
 * This kernel is invoked by NonbondedSoftcoreForce to calculate the forces acting on the system.
 */
class CudaFreeEnergyCalcNonbondedSoftcoreForceKernel : public CalcNonbondedSoftcoreForceKernel {
public:
    CudaFreeEnergyCalcNonbondedSoftcoreForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, System& system) :
             CalcNonbondedSoftcoreForceKernel(name, platform), data(data), system(system) {
        gpuNonbondedSoftcore = NULL;
        gpuLJ14Softcore      = NULL;
#ifdef FreeEnergyDebug
        log                  = stderr;
#else
        log                  = NULL;
#endif
        setSim               = 0;
        numExceptions        = 0;
        numParticles         = 0;
        bIncludeGBSA         = false;
        bIncludeGBVI         = false;
        includeSoftcore      = false;
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
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    void executeForces(ContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the NonbondedForce
     */
    double executeEnergy(ContextImpl& context);
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
    /**
     * Get GpuLJ14Softcore
     *
     * @return GpuLJ14Softcore object
     */
    GpuLJ14Softcore* getGpuLJ14Softcore( void ) const;
    /**
     * Set GpuLJ14Softcore
     *
     * @param GpuLJ14Softcore object
     */
    void setGpuLJ14Softcore( GpuLJ14Softcore* gpuLJ14Softcore );
private:
    CudaPlatform::PlatformData& data;
    int numParticles;
    System& system;
    GpuNonbondedSoftcore* gpuNonbondedSoftcore;
    GpuLJ14Softcore* gpuLJ14Softcore;
    bool bIncludeGBSA;
    bool bIncludeGBVI;
    int includeSoftcore;
    int numExceptions;
    FILE* log;
    int setSim;
};

/**
 * This kernel is invoked by GBSAOBCForce to calculate the forces acting on the system.
 */
class CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel : public CalcGBSAOBCSoftcoreForceKernel {
public:
    CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data) :
       CalcGBSAOBCSoftcoreForceKernel(name, platform), data(data) {
#ifdef FreeEnergyDebug
        log                  = stderr;
#else
        log                  = NULL;
#endif
        setSim               = 0;
        gpuObcGbsaSoftcore   = NULL;
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
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    void executeForces(ContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the GBSAOBCForce
     */
    double executeEnergy(ContextImpl& context);
private:
    CudaPlatform::PlatformData& data;
    FILE* log;
    int setSim;
    GpuObcGbsaSoftcore* gpuObcGbsaSoftcore;
};

/**
 * This kernel is invoked by GBVIForce to calculate the forces acting on the system.
 */
class CudaFreeEnergyCalcGBVISoftcoreForceKernel : public CalcGBVISoftcoreForceKernel {
public:
    CudaFreeEnergyCalcGBVISoftcoreForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data) :
         CalcGBVISoftcoreForceKernel(name, platform), data(data) {
#ifdef FreeEnergyDebug
        log                  = stderr;
#else
        log                  = NULL;
#endif
        setSim               = 0;
        quinticScaling       = 0;
        gpuGBVISoftcore      = NULL;
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
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    void executeForces(ContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the GBVIForce
     */
    double executeEnergy(ContextImpl& context);

    /**
     * Apply quintic scaling for Born radii
     * 
     * @return nonzero value if scaling is to be applied
     */
    int getQuinticScaling(void) const;

    /**
     * Set flag for quintic scaling for Born radii
     * 
     * @param nonzero value if scaling is to be applied
     */
    void setQuinticScaling(int quinticScaling );

private:
    CudaPlatform::PlatformData& data;
    GpuGBVISoftcore* gpuGBVISoftcore;
    FILE* log;
    int setSim;
    int quinticScaling;
};

} // namespace OpenMM

#endif /*OPENMM_FREE_ENERGY_CUDA_KERNELS_H_*/

#ifndef AMOEBA_OPENMM_KERNELS_H_
#define AMOEBA_OPENMM_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                             OpenMMAmoeba                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2018 Stanford University and the Authors.      *
 * Authors:                                                                   *
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

#include "OpenMMAmoeba.h"
#include "openmm/KernelImpl.h"
#include "openmm/System.h"
#include "openmm/Platform.h"

#include <set>
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This kernel is invoked by AmoebaBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcAmoebaBondForceKernel : public KernelImpl {

public:

    static std::string Name() {
        return "CalcAmoebaBondForce";
    }

    CalcAmoebaBondForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaBondForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaBondForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaBondForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const AmoebaBondForce& force) = 0;
};

/**
 * This kernel is invoked by AmoebaAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcAmoebaAngleForceKernel : public KernelImpl {

public:

    static std::string Name() {
        return "CalcAmoebaAngleForce";
    }

    CalcAmoebaAngleForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaAngleForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaAngleForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaAngleForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const AmoebaAngleForce& force) = 0;
};

/**
 * This kernel is invoked by AmoebaInPlaneAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcAmoebaInPlaneAngleForceKernel : public KernelImpl {

public:

    static std::string Name() {
        return "CalcAmoebaInPlaneAngleForce";
    }

    CalcAmoebaInPlaneAngleForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaInPlaneAngleForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaInPlaneAngleForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaInPlaneAngleForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const AmoebaInPlaneAngleForce& force) = 0;
};

/**
 * This kernel is invoked by AmoebaTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcAmoebaPiTorsionForceKernel : public KernelImpl {

public:

    static std::string Name() {
        return "CalcAmoebaPiTorsionForce";
    }

    CalcAmoebaPiTorsionForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the PiTorsionForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaPiTorsionForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaPiTorsionForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const AmoebaPiTorsionForce& force) = 0;
};

/**
 * This kernel is invoked by AmoebaTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcAmoebaStretchBendForceKernel : public KernelImpl {

public:

    static std::string Name() {
        return "CalcAmoebaStretchBendForce";
    }

    CalcAmoebaStretchBendForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the StretchBendForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaStretchBendForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaStretchBendForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const AmoebaStretchBendForce& force) = 0;
};

/**
 * This kernel is invoked by AmoebaTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcAmoebaOutOfPlaneBendForceKernel : public KernelImpl {

public:

    static std::string Name() {
        return "CalcAmoebaOutOfPlaneBendForce";
    }

    CalcAmoebaOutOfPlaneBendForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the OutOfPlaneBendForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaOutOfPlaneBendForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaOutOfPlaneBendForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const AmoebaOutOfPlaneBendForce& force) = 0;
};

/**
 * This kernel is invoked by AmoebaTorsionTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcAmoebaTorsionTorsionForceKernel : public KernelImpl {

public:

    static std::string Name() {
        return "CalcAmoebaTorsionTorsionForce";
    }

    CalcAmoebaTorsionTorsionForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the TorsionTorsionForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaTorsionTorsionForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
};

/**
 * This kernel is invoked by AmoebaMultipoleForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcAmoebaMultipoleForceKernel : public KernelImpl {

public:

    static std::string Name() {
        return "CalcAmoebaMultipoleForce";
    }

    CalcAmoebaMultipoleForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the MultipoleForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaMultipoleForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;

    virtual void getLabFramePermanentDipoles(ContextImpl& context, std::vector<Vec3>& dipoles) = 0;
    virtual void getInducedDipoles(ContextImpl& context, std::vector<Vec3>& dipoles) = 0;
    virtual void getTotalDipoles(ContextImpl& context, std::vector<Vec3>& dipoles) = 0;

    virtual void getElectrostaticPotential(ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                           std::vector< double >& outputElectrostaticPotential) = 0;

    virtual void getSystemMultipoleMoments(ContextImpl& context, std::vector< double >& outputMultipoleMoments) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaMultipoleForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const AmoebaMultipoleForce& force) = 0;

    /**
     * Get the parameters being used for PME.
     *
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    virtual void getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const = 0;
};

/**
 * This kernel is invoked by AmoebaGeneralizedKirkwoodForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcAmoebaGeneralizedKirkwoodForceKernel : public KernelImpl {

public:

    static std::string Name() {
        return "CalcAmoebaGeneralizedKirkwoodForce";
    }

    CalcAmoebaGeneralizedKirkwoodForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the GBSAOBCForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaGeneralizedKirkwoodForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaGeneralizedKirkwoodForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const AmoebaGeneralizedKirkwoodForce& force) = 0;
};


/**
 * This kernel is invoked by AmoebaVdwForce to calculate the vdw forces acting on the system and the vdw energy of the system.
 */
class CalcAmoebaVdwForceKernel : public KernelImpl {
public:

    static std::string Name() {
        return "CalcAmoebaVdwForce";
    }

    CalcAmoebaVdwForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the GBSAOBCForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaVdwForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaVdwForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const AmoebaVdwForce& force) = 0;
};

/**
 * This kernel is invoked by AmoebaWcaDispersionForce to calculate the WCA dispersion forces acting on the system and the WCA dispersion energy of the system.
 */
class CalcAmoebaWcaDispersionForceKernel : public KernelImpl {

public:

    static std::string Name() {
        return "CalcAmoebaWcaDispersionForce";
    }

    CalcAmoebaWcaDispersionForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the GBSAOBCForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaWcaDispersionForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaWcaDispersionForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const AmoebaWcaDispersionForce& force) = 0;
};

/**
 * This kernel is invoked by HippoNonbondedForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcHippoNonbondedForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcHippoNonbondedForce";
    }
    CalcHippoNonbondedForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the MultipoleForce this kernel will be used for
     */
    virtual void initialize(const System& system, const HippoNonbondedForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    virtual void getLabFramePermanentDipoles(ContextImpl& context, std::vector<Vec3>& dipoles) = 0;
    virtual void getInducedDipoles(ContextImpl& context, std::vector<Vec3>& dipoles) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaMultipoleForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const HippoNonbondedForce& force) = 0;
    /**
     * Get the parameters being used for PME.
     *
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    virtual void getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const = 0;
    /**
     * Get the parameters being used for dispersion PME.
     *
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    virtual void getDPMEParameters(double& alpha, int& nx, int& ny, int& nz) const = 0;
};

} // namespace OpenMM

#endif /*AMOEBA_OPENMM_KERNELS_H*/

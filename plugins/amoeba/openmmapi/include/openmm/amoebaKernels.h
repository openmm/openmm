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
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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
 * This kernel is invoked by AmoebaHarmonicBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcAmoebaHarmonicBondForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcAmoebaHarmonicBondForce";
    }
    CalcAmoebaHarmonicBondForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaHarmonicBondForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaHarmonicBondForce& force) = 0;
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

class CalcAmoebaUreyBradleyForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcAmoebaUreyBradleyForce";
    }

    CalcAmoebaUreyBradleyForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaUreyBradleyForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaUreyBradleyForce& force) = 0;
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
 * This kernel is invoked by AmoebaHarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcAmoebaHarmonicAngleForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcAmoebaHarmonicAngleForce";
    }
    CalcAmoebaHarmonicAngleForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the HarmonicAngleForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaHarmonicAngleForce& force) = 0;
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
 * This kernel is invoked by AmoebaHarmonicInPlaneAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcAmoebaHarmonicInPlaneAngleForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcAmoebaHarmonicInPlaneAngleForce";
    }
    CalcAmoebaHarmonicInPlaneAngleForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the HarmonicInPlaneAngleForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaHarmonicInPlaneAngleForce& force) = 0;
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
 * This kernel is invoked by AmoebaTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcAmoebaTorsionForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcAmoebaTorsionForce";
    }
    CalcAmoebaTorsionForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the TorsionForce this kernel will be used for
     */
    virtual void initialize(const System& system, const AmoebaTorsionForce& force) = 0;
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
};

} // namespace OpenMM

#endif /*AMOEBA_OPENMM_KERNELS_H*/

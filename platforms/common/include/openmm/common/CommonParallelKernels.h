#ifndef OPENMM_COMMONPARALLELKERNELS_H_
#define OPENMM_COMMONPARALLELKERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2024 Stanford University and the Authors.      *
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

#include "openmm/common/ComputeContext.h"
#include "openmm/common/CommonKernels.h"

namespace OpenMM {

    /**
     * This kernel is invoked by HarmonicBondForce to calculate the forces acting on the system and the energy of the system.
     */
    class CommonParallelCalcHarmonicBondForceKernel : public CalcHarmonicBondForceKernel {
    public:
        CommonParallelCalcHarmonicBondForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system);

        CommonCalcHarmonicBondForceKernel& getKernel(int index) {
            return dynamic_cast<CommonCalcHarmonicBondForceKernel&> (kernels[index].getImpl());
        }
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the HarmonicBondForce this kernel will be used for
         */
        void initialize(const System& system, const HarmonicBondForce& force);
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
         * Copy changed parameters over to a context.
         *
         * @param firstBond  the index of the first bond whose parameters might have changed
         * @param lastBond   the index of the last bond whose parameters might have changed
         */
        void copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force, int firstBond, int lastBond);
    private:
        class Task;
        ComputeContext& cc;
        std::vector<Kernel> kernels;
    };

    /**
     * This kernel is invoked by CustomBondForce to calculate the forces acting on the system and the energy of the system.
     */
    class CommonParallelCalcCustomBondForceKernel : public CalcCustomBondForceKernel {
    public:
        CommonParallelCalcCustomBondForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system);

        CommonCalcCustomBondForceKernel& getKernel(int index) {
            return dynamic_cast<CommonCalcCustomBondForceKernel&> (kernels[index].getImpl());
        }
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the CustomBondForce this kernel will be used for
         */
        void initialize(const System& system, const CustomBondForce& force);
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
         * Copy changed parameters over to a context.
         *
         * @param firstBond  the index of the first bond whose parameters might have changed
         * @param lastBond   the index of the last bond whose parameters might have changed
         */
        void copyParametersToContext(ContextImpl& context, const CustomBondForce& force, int firstBond, int lastBond);
    private:
        class Task;
        ComputeContext& cc;
        std::vector<Kernel> kernels;
    };

    /**
     * This kernel is invoked by HarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
     */
    class CommonParallelCalcHarmonicAngleForceKernel : public CalcHarmonicAngleForceKernel {
    public:
        CommonParallelCalcHarmonicAngleForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system);

        CommonCalcHarmonicAngleForceKernel& getKernel(int index) {
            return dynamic_cast<CommonCalcHarmonicAngleForceKernel&> (kernels[index].getImpl());
        }
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the HarmonicAngleForce this kernel will be used for
         */
        void initialize(const System& system, const HarmonicAngleForce& force);
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
         * Copy changed parameters over to a context.
         *
         * @param firstAngle the index of the first bond whose parameters might have changed
         * @param lastAngle  the index of the last bond whose parameters might have changed
         */
        void copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force, int firstAngle, int lastAngle);
    private:
        class Task;
        ComputeContext& cc;
        std::vector<Kernel> kernels;
    };

    /**
     * This kernel is invoked by CustomAngleForce to calculate the forces acting on the system and the energy of the system.
     */
    class CommonParallelCalcCustomAngleForceKernel : public CalcCustomAngleForceKernel {
    public:
        CommonParallelCalcCustomAngleForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system);

        CommonCalcCustomAngleForceKernel& getKernel(int index) {
            return dynamic_cast<CommonCalcCustomAngleForceKernel&> (kernels[index].getImpl());
        }
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the CustomAngleForce this kernel will be used for
         */
        void initialize(const System& system, const CustomAngleForce& force);
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
         * Copy changed parameters over to a context.
         *
         * @param firstAngle the index of the first bond whose parameters might have changed
         * @param lastAngle  the index of the last bond whose parameters might have changed
         */
        void copyParametersToContext(ContextImpl& context, const CustomAngleForce& force, int firstAngle, int lastAngle);
    private:
        class Task;
        ComputeContext& cc;
        std::vector<Kernel> kernels;
    };

    /**
     * This kernel is invoked by PeriodicTorsionForce to calculate the forces acting on the system and the energy of the system.
     */
    class CommonParallelCalcPeriodicTorsionForceKernel : public CalcPeriodicTorsionForceKernel {
    public:
        CommonParallelCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system);

        CommonCalcPeriodicTorsionForceKernel& getKernel(int index) {
            return dynamic_cast<CommonCalcPeriodicTorsionForceKernel&> (kernels[index].getImpl());
        }
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the PeriodicTorsionForce this kernel will be used for
         */
        void initialize(const System& system, const PeriodicTorsionForce& force);
        /**
         * Execute the kernel to calculate the forces and/or energy.
         *
         * @param context        the context in which to execute this kernel
         * @param includeForces  true if forces should be calculated
         * @param includeEnergy  true if the energy should be calculated
         * @return the potential energy due to the force
         */
        double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
        class Task;
        /**
         * Copy changed parameters over to a context.
         *
         * @param context      the context to copy parameters to
         * @param force        the PeriodicTorsionForce to copy the parameters from
         * @param firstTorsion the index of the first torsion whose parameters might have changed
         * @param lastTorsion  the index of the last torsion whose parameters might have changed
         */
        void copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force, int firstTorsion, int lastTorsion);
    private:
        ComputeContext& cc;
        std::vector<Kernel> kernels;
    };

    /**
     * This kernel is invoked by RBTorsionForce to calculate the forces acting on the system and the energy of the system.
     */
    class CommonParallelCalcRBTorsionForceKernel : public CalcRBTorsionForceKernel {
    public:
        CommonParallelCalcRBTorsionForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system);

        CommonCalcRBTorsionForceKernel& getKernel(int index) {
            return dynamic_cast<CommonCalcRBTorsionForceKernel&> (kernels[index].getImpl());
        }
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the RBTorsionForce this kernel will be used for
         */
        void initialize(const System& system, const RBTorsionForce& force);
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
         * Copy changed parameters over to a context.
         *
         * @param context    the context to copy parameters to
         * @param force      the RBTorsionForce to copy the parameters from
         */
        void copyParametersToContext(ContextImpl& context, const RBTorsionForce& force);
    private:
        class Task;
        ComputeContext& cc;
        std::vector<Kernel> kernels;
    };

    /**
     * This kernel is invoked by CMAPTorsionForce to calculate the forces acting on the system and the energy of the system.
     */
    class CommonParallelCalcCMAPTorsionForceKernel : public CalcCMAPTorsionForceKernel {
    public:
        CommonParallelCalcCMAPTorsionForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system);

        CommonCalcCMAPTorsionForceKernel& getKernel(int index) {
            return dynamic_cast<CommonCalcCMAPTorsionForceKernel&> (kernels[index].getImpl());
        }
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the CMAPTorsionForce this kernel will be used for
         */
        void initialize(const System& system, const CMAPTorsionForce& force);
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
         * Copy changed parameters over to a context.
         *
         * @param context    the context to copy parameters to
         * @param force      the CMAPTorsionForce to copy the parameters from
         */
        void copyParametersToContext(ContextImpl& context, const CMAPTorsionForce& force);
    private:
        class Task;
        ComputeContext& cc;
        std::vector<Kernel> kernels;
    };

    /**
     * This kernel is invoked by CustomTorsionForce to calculate the forces acting on the system and the energy of the system.
     */
    class CommonParallelCalcCustomTorsionForceKernel : public CalcCustomTorsionForceKernel {
    public:
        CommonParallelCalcCustomTorsionForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system);

        CommonCalcCustomTorsionForceKernel& getKernel(int index) {
            return dynamic_cast<CommonCalcCustomTorsionForceKernel&> (kernels[index].getImpl());
        }
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the CustomTorsionForce this kernel will be used for
         */
        void initialize(const System& system, const CustomTorsionForce& force);
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
         * Copy changed parameters over to a context.
         *
         * @param context      the context to copy parameters to
         * @param force        the CustomTorsionForce to copy the parameters from
         * @param firstTorsion the index of the first torsion whose parameters might have changed
         * @param lastTorsion  the index of the last torsion whose parameters might have changed
         */
        void copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force, int firstTorsion, int lastTorsion);
    private:
        class Task;
        ComputeContext& cc;
        std::vector<Kernel> kernels;
    };

    /**
     * This kernel is invoked by CustomNonbondedForce to calculate the forces acting on the system.
     */
    class CommonParallelCalcCustomNonbondedForceKernel : public CalcCustomNonbondedForceKernel {
    public:
        CommonParallelCalcCustomNonbondedForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system);

        CommonCalcCustomNonbondedForceKernel& getKernel(int index) {
            return dynamic_cast<CommonCalcCustomNonbondedForceKernel&> (kernels[index].getImpl());
        }
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the CustomNonbondedForce this kernel will be used for
         */
        void initialize(const System& system, const CustomNonbondedForce& force);
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
         * Copy changed parameters over to a context.
         *
         * @param context    the context to copy parameters to
         * @param force      the CustomNonbondedForce to copy the parameters from
         * @param firstParticle  the index of the first particle whose parameters might have changed
         * @param lastParticle   the index of the last particle whose parameters might have changed
         */
        void copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force, int firstParticle, int lastParticle);
    private:
        class Task;
        ComputeContext& cc;
        std::vector<Kernel> kernels;
    };

    /**
     * This kernel is invoked by CustomExternalForce to calculate the forces acting on the system and the energy of the system.
     */
    class CommonParallelCalcCustomExternalForceKernel : public CalcCustomExternalForceKernel {
    public:
        CommonParallelCalcCustomExternalForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system);

        CommonCalcCustomExternalForceKernel& getKernel(int index) {
            return dynamic_cast<CommonCalcCustomExternalForceKernel&> (kernels[index].getImpl());
        }
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the CustomExternalForce this kernel will be used for
         */
        void initialize(const System& system, const CustomExternalForce& force);
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
         * Copy changed parameters over to a context.
         *
         * @param context        the context to copy parameters to
         * @param force          the CustomExternalForce to copy the parameters from
         * @param firstParticle  the index of the first particle whose parameters might have changed
         * @param lastParticle   the index of the last particle whose parameters might have changed
         */
        void copyParametersToContext(ContextImpl& context, const CustomExternalForce& force, int firstParticle, int lastParticle);
    private:
        class Task;
        ComputeContext& cc;
        std::vector<Kernel> kernels;
    };

    /**
     * This kernel is invoked by CustomHbondForce to calculate the forces acting on the system.
     */
    class CommonParallelCalcCustomHbondForceKernel : public CalcCustomHbondForceKernel {
    public:
        CommonParallelCalcCustomHbondForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system);

        CommonCalcCustomHbondForceKernel& getKernel(int index) {
            return dynamic_cast<CommonCalcCustomHbondForceKernel&> (kernels[index].getImpl());
        }
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the CustomHbondForce this kernel will be used for
         */
        void initialize(const System& system, const CustomHbondForce& force);
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
         * Copy changed parameters over to a context.
         *
         * @param context    the context to copy parameters to
         * @param force      the CustomHbondForce to copy the parameters from
         */
        void copyParametersToContext(ContextImpl& context, const CustomHbondForce& force);
    private:
        class Task;
        ComputeContext& cc;
        std::vector<Kernel> kernels;
    };

    /**
     * This kernel is invoked by CustomCompoundBondForce to calculate the forces acting on the system.
     */
    class CommonParallelCalcCustomCompoundBondForceKernel : public CalcCustomCompoundBondForceKernel {
    public:
        CommonParallelCalcCustomCompoundBondForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system);

        CommonCalcCustomCompoundBondForceKernel& getKernel(int index) {
            return dynamic_cast<CommonCalcCustomCompoundBondForceKernel&> (kernels[index].getImpl());
        }
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the CustomCompoundBondForce this kernel will be used for
         */
        void initialize(const System& system, const CustomCompoundBondForce& force);
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
         * Copy changed parameters over to a context.
         *
         * @param context    the context to copy parameters to
         * @param force      the CustomCompoundBondForce to copy the parameters from
         */
        void copyParametersToContext(ContextImpl& context, const CustomCompoundBondForce& force);
    private:
        class Task;
        ComputeContext& cc;
        std::vector<Kernel> kernels;
    };

} // namespace OpenMM

#endif /*OPENMM_COMMONPARALLELKERNELS_H_*/

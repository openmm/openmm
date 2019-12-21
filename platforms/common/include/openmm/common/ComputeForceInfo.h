#ifndef OPENMM_COMPUTEFORCEINFO_H_
#define OPENMM_COMPUTEFORCEINFO_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2019 Stanford University and the Authors.       *
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

#include "openmm/common/windowsExportCommon.h"
#include <vector>

namespace OpenMM {

/**
 * ComputeForceInfo objects describe information about the behavior and requirements of
 * a force.  They exist primarily to help a ComputeContext determine how particles can be
 * reordered without affecting forces.  Force kernels create them during initialization
 * and add them to the ComputeContext by calling addForce().
 */

class OPENMM_EXPORT_COMMON ComputeForceInfo {
public:
    ComputeForceInfo() {
    }
    /**
     * Get whether or not two particles have identical force field parameters.
     */
    virtual bool areParticlesIdentical(int particle1, int particle2);
    /**
     * Get the number of particle groups defined by this force.
     */
    virtual int getNumParticleGroups();
    /**
     * Get the list of particles in a particular group.
     */
    virtual void getParticlesInGroup(int index, std::vector<int>& particles);
    /**
     * Get whether two particle groups are identical.
     */
    virtual bool areGroupsIdentical(int group1, int group2);
};

} // namespace OpenMM

#endif /*OPENMM_COMPUTEFORCEINFO_H_*/

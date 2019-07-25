#ifndef OPENMM_OPENCLFORCEINFO_H_
#define OPENMM_OPENCLFORCEINFO_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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

#include "windowsExportOpenCL.h"
#include <vector>

namespace OpenMM {

/**
 * This class is used by the OpenCL implementation of a Force class to convey information
 * about the behavior and requirements of that force.
 */

class OPENMM_EXPORT_OPENCL OpenCLForceInfo {
public:
    OpenCLForceInfo(int requiredForceBuffers) : requiredForceBuffers(requiredForceBuffers) {
    }
    /**
     * Get the number of force buffers this force requires.
     */
    int getRequiredForceBuffers() {
        return requiredForceBuffers;
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
private:
    int requiredForceBuffers;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLFORCEINFO_H_*/

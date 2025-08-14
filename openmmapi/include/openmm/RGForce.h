#ifndef OPENMM_RGFORCE_H_
#define OPENMM_RGFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "Force.h"
#include "Vec3.h"
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This is a force whose energy equals the radius of gyration (Rg) of a set of
 * particles.  It is intended for use with CustomCVForce.  You will not normally
 * want a force that exactly equals the radius of gyration, but there are many
 * situations where it is useful to have a restraining or biasing force that
 * depends on Rg in some way.
 */

class OPENMM_EXPORT RGForce : public Force {
public:
    /**
     * Create an RGForce.
     *
     * @param particles    the indices of the particles to use when computing
     *                     Rg.  If this is empty (the default), all particles in
     *                     the system will be used.
     */
    explicit RGForce(const std::vector<int>& particles=std::vector<int>());
    /**
     * Get the indices of the particles to use when computing Rg.  If this
     * is empty, all particles in the system will be used.
     */
    const std::vector<int>& getParticles() const {
        return particles;
    }
    /**
     * Set the indices of the particles to use when computing Rg.  If this
     * is empty, all particles in the system will be used.
     */
    void setParticles(const std::vector<int>& particles);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }
protected:
    ForceImpl* createImpl() const;
private:
    std::vector<int> particles;
};

} // namespace OpenMM

#endif /*OPENMM_RGFORCE_H_*/


/* Portions copyright (c) 2009-2013 Stanford University and Simbios.
 * Contributors: Peter Eastman
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef __ReferenceConstraintAlgorithm_H__
#define __ReferenceConstraintAlgorithm_H__

#include "openmm/Vec3.h"
#include "openmm/internal/windowsExport.h"
#include "SimTKOpenMMRealType.h"
#include <vector>

namespace OpenMM {

/**
 * This abstract class defines the interface which constraint algorithms must implement.
 */
class OPENMM_EXPORT ReferenceConstraintAlgorithm {
public:

    virtual ~ReferenceConstraintAlgorithm() {};

    /**
     * Apply the constraint algorithm.
     * 
     * @param atomCoordinates  the original atom coordinates
     * @param atomCoordinatesP the new atom coordinates
     * @param inverseMasses    1/mass
     * @param tolerance        the constraint tolerance
     */
    virtual void apply(std::vector<OpenMM::Vec3>& atomCoordinates,
                     std::vector<OpenMM::Vec3>& atomCoordinatesP, std::vector<double>& inverseMasses, double tolerance) = 0;

    /**
     * Apply the constraint algorithm to velocities.
     * 
     * @param atomCoordinates  the atom coordinates
     * @param atomCoordinatesP the velocities to modify
     * @param inverseMasses    1/mass
     * @param tolerance        the constraint tolerance
     */
    virtual void applyToVelocities(std::vector<OpenMM::Vec3>& atomCoordinates,
                     std::vector<OpenMM::Vec3>& velocities, std::vector<double>& inverseMasses, double tolerance) = 0;
};

} // namespace OpenMM

#endif // __ReferenceConstraintAlgorithm_H__

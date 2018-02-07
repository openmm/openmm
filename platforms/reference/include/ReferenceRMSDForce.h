/* Portions copyright (c) 2018 Stanford University and Simbios.
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

#ifndef __ReferenceRMSDForce_H__
#define __ReferenceRMSDForce_H__

#include "openmm/RMSDForce.h"
#include <vector>

namespace OpenMM {

class ReferenceRMSDForce {
private:
    std::vector<OpenMM::Vec3> referencePos;
    std::vector<int> particles;

public:
    /**
     * Constructor
     */
    ReferenceRMSDForce(std::vector<OpenMM::Vec3>& referencePos, std::vector<int>& particles);

    /**
     * Destructor
     */
    ~ReferenceRMSDForce();

    /**
     * Calculate the interaction.
     * 
     * @param atomCoordinates    atom coordinates
     * @param forces             the forces are added to this
     * @return the energy of the interaction
     */
   double calculateIxn(std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& forces) const;
};

} // namespace OpenMM

#endif // __ReferenceRMSDForce_H__

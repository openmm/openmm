#ifndef OPENMM_REFERENCECONSTANTPOTENTIAL14_H_
#define OPENMM_REFERENCECONSTANTPOTENTIAL14_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2006-2025 Stanford University and the Authors.      *
 * Authors: Pande Group, Evan Pretti                                          *
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

#include "ReferenceBondIxn.h"
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

class OPENMM_EXPORT ReferenceConstantPotential14 : public ReferenceBondIxn {

public:
     ReferenceConstantPotential14();
     ~ReferenceConstantPotential14();

    /**
     * Sets the force to use periodic boundary conditions with the specified
     * vectors.
     */
    void setPeriodic(OpenMM::Vec3* vectors);

    /**
     * Calculates 1-4 nonbonded interactions (i.e., exceptions with a non-zero
     * charge product that should behave effectively as bonded interactions.
     * parameters should contain a single item (the charge product).
     */
    void calculateBondIxn(std::vector<int>& atomIndices, std::vector<OpenMM::Vec3>& atomCoordinates,
                          std::vector<double>& parameters, std::vector<OpenMM::Vec3>& forces,
                          double* totalEnergy, double* energyParamDerivs);

private:
    bool periodic;
    OpenMM::Vec3 periodicBoxVectors[3];
};

} // namespace OpenMM

#endif // OPENMM_REFERENCECONSTANTPOTENTIAL14_H_

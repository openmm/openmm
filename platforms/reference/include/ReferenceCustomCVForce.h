
/* Portions copyright (c) 2017 Stanford University and Simbios.
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

#ifndef __ReferenceCustomCVForce_H__
#define __ReferenceCustomCVForce_H__

#include "openmm/CustomCVForce.h"
#include "openmm/internal/ContextImpl.h"
#include "lepton/ExpressionProgram.h"
#include <map>
#include <string>
#include <vector>

namespace OpenMM {

class ReferenceCustomCVForce {
private:
    Lepton::ExpressionProgram energyExpression;
    std::vector<std::string> variableNames, paramDerivNames;
    std::vector<Lepton::ExpressionProgram> variableDerivExpressions;
    std::vector<Lepton::ExpressionProgram> paramDerivExpressions;

public:
    /**
     * Constructor
     */
    ReferenceCustomCVForce(const OpenMM::CustomCVForce& force);

    /**
     * Destructor
     */
    ~ReferenceCustomCVForce();

    /**
     * Update any tabulated functions used by the force.  This is called when the user calls
     * updateParametersInContext().
     */
    void updateTabulatedFunctions(const OpenMM::CustomCVForce& force);

    /**
     * Calculate the interaction.
     * 
     * @param innerContext       the context created by the force for evaluating collective variables
     * @param atomCoordinates    atom coordinates
     * @param globalParameters   the values of global parameters
     * @param forces             the forces are added to this
     * @param totalEnergy        the energy is added to this
     * @param energyParamDerivs  parameter derivatives are added to this
     */
   void calculateIxn(ContextImpl& innerContext, std::vector<OpenMM::Vec3>& atomCoordinates,
                     const std::map<std::string, double>& globalParameters,
                     std::vector<OpenMM::Vec3>& forces, double* totalEnergy, std::map<std::string, double>& energyParamDerivs) const;
};

} // namespace OpenMM

#endif // __ReferenceCustomCVForce_H__

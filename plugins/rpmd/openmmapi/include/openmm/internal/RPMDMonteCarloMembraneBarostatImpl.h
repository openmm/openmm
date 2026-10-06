#ifndef OPENMM_RPMDMONTECARLOMEMBRANEBAROSTATIMPL_H_
#define OPENMM_RPMDMONTECARLOMEMBRANEBAROSTATIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2010-2026 Stanford University and the Authors.      *
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

#include "openmm/RPMDMonteCarloMembraneBarostat.h"
#include "openmm/RPMDUpdater.h"
#include "openmm/Kernel.h"
#include "openmm/Vec3.h"
#include "sfmt/SFMT.h"
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This is the internal implementation of RPMDMonteCarloMembraneBarostat.
 */

class RPMDMonteCarloMembraneBarostatImpl : public RPMDUpdater {
public:
    RPMDMonteCarloMembraneBarostatImpl(const RPMDMonteCarloMembraneBarostat& owner);
    void initialize(ContextImpl& context);
    const RPMDMonteCarloMembraneBarostat& getOwner() const {
        return owner;
    }
    void updateRPMDState(ContextImpl& context);
    void updateContextState(ContextImpl& context, bool& forcesInvalid) {
        // This is unused, since the updating is done in updateRPMDState().
    }
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
        // This force doesn't apply forces to particles.
        return 0.0;
    }
    std::map<std::string, double> getDefaultParameters();
    std::vector<std::string> getKernelNames();
private:
    const RPMDMonteCarloMembraneBarostat& owner;
    int step, numAttempted[3], numAccepted[3];
    double volumeScale[3];
    std::vector<std::vector<Vec3> > savedPositions;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_RPMDMONTECARLOMEMBRANEBAROSTATIMPL_H_*/

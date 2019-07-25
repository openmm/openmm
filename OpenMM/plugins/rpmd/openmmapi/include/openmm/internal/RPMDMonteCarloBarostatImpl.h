#ifndef OPENMM_RPMDMONTECARLOBAROSTATIMPL_H_
#define OPENMM_RPMDMONTECARLOBAROSTATIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2015 Stanford University and the Authors.      *
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

#include "openmm/RPMDMonteCarloBarostat.h"
#include "openmm/RPMDUpdater.h"
#include "openmm/Kernel.h"
#include "openmm/Vec3.h"
#include "sfmt/SFMT.h"
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This is the internal implementation of RPMDMonteCarloBarostat.
 */

class RPMDMonteCarloBarostatImpl : public RPMDUpdater {
public:
    RPMDMonteCarloBarostatImpl(const RPMDMonteCarloBarostat& owner);
    void initialize(ContextImpl& context);
    const RPMDMonteCarloBarostat& getOwner() const {
        return owner;
    }
    void updateRPMDState(ContextImpl& context);
    void updateContextState(ContextImpl& context) {
        // This is unused, since the updating is done in updateRPMDState().
    }
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
        // This force doesn't apply forces to particles.
        return 0.0;
    }
    std::map<std::string, double> getDefaultParameters();
    std::vector<std::string> getKernelNames();
private:
    const RPMDMonteCarloBarostat& owner;
    int step, numAttempted, numAccepted;
    double volumeScale;
    OpenMM_SFMT::SFMT random;
    std::vector<std::vector<Vec3> > savedPositions;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_RPMDMONTECARLOBAROSTATIMPL_H_*/

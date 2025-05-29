#ifndef OPENMM_LANGEVININTEGRATOR_H_
#define OPENMM_LANGEVININTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
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

#include "LangevinMiddleIntegrator.h"
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This is an Integrator which simulates a System using Langevin dynamics.  It
 * is identical to LangevinMiddleIntegrator.  (It simply subclasses it without
 * adding or changing anything.)  The two different classes exist only for
 * historical reasons.
 */

class OPENMM_EXPORT LangevinIntegrator : public LangevinMiddleIntegrator {
public:
    /**
     * Create a LangevinIntegrator.
     * 
     * @param temperature    the temperature of the heat bath (in Kelvin)
     * @param frictionCoeff  the friction coefficient which couples the system to the heat bath (in inverse picoseconds)
     * @param stepSize       the step size with which to integrate the system (in picoseconds)
     */
    LangevinIntegrator(double temperature, double frictionCoeff, double stepSize);
};

} // namespace OpenMM

#endif /*OPENMM_LANGEVININTEGRATOR_H_*/

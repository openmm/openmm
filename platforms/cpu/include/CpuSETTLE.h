#ifndef OPENMM_CPUSETTLE_H_
#define OPENMM_CPUSETTLE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2017 Stanford University and the Authors.      *
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

#include "ReferenceSETTLEAlgorithm.h"
#include "windowsExportCpu.h"
#include "openmm/System.h"
#include "openmm/internal/ThreadPool.h"
#include <vector>

namespace OpenMM {

/**
 * This class uses multiple ReferenceSETTLEAlgorithm objects to execute the algorithm in parallel.
 */
class OPENMM_EXPORT_CPU CpuSETTLE : public ReferenceConstraintAlgorithm {
public:
    CpuSETTLE(const System& system, const ReferenceSETTLEAlgorithm& settle, ThreadPool& threads);
    ~CpuSETTLE();

    /**
     * Apply the constraint algorithm.
     * 
     * @param atomCoordinates  the original atom coordinates
     * @param atomCoordinatesP the new atom coordinates
     * @param inverseMasses    1/mass
     * @param tolerance        the constraint tolerance
     */
    void apply(std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& atomCoordinatesP, std::vector<double>& inverseMasses, double tolerance);

    /**
     * Apply the constraint algorithm to velocities.
     * 
     * @param atomCoordinates  the atom coordinates
     * @param atomCoordinatesP the velocities to modify
     * @param inverseMasses    1/mass
     * @param tolerance        the constraint tolerance
     */
    void applyToVelocities(std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& velocities, std::vector<double>& inverseMasses, double tolerance);
private:
    std::vector<ReferenceSETTLEAlgorithm*> threadSettle;
    ThreadPool& threads;
};

} // namespace OpenMM

#endif /*OPENMM_CPUSETTLE_H_*/

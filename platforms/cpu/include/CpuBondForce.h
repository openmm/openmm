#ifndef OPENMM_CPUBONDFORCE_H_
#define OPENMM_CPUBONDFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014-2018 Stanford University and the Authors.      *
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

#include "ReferenceBondIxn.h"
#include "windowsExportCpu.h"
#include "openmm/internal/ThreadPool.h"
#include <list>
#include <set>
#include <vector>

namespace OpenMM {

/**
 * This class parallelizes the calculation of bonded forces.
 */
class OPENMM_EXPORT_CPU CpuBondForce {
public:
    CpuBondForce();
    /**
     * Analyze the set of bonds and decide which to compute with each thread.
     */
    void initialize(int numAtoms, int numBonds, int numAtomsPerBond, std::vector<std::vector<int> >& bondAtoms, ThreadPool& threads);
    /**
     * Compute the forces from all bonds.
     */
    void calculateForce(std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<std::vector<double> >& parameters, std::vector<OpenMM::Vec3>& forces, 
            double* totalEnergy, ReferenceBondIxn& referenceBondIxn);
    /**
     * This routine contains the code executed by each thread.
     */
    void threadComputeForce(ThreadPool& threads, int threadIndex, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<std::vector<double> >& parameters,
            std::vector<OpenMM::Vec3>& forces, double* totalEnergy, ReferenceBondIxn& referenceBondIxn);
private:
    bool canAssignBond(int bond, int thread, std::vector<int>& atomThread);
    void assignBond(int bond, int thread, std::vector<int>& atomThread, std::vector<int>& bondThread, std::vector<std::set<int> >& atomBonds, std::list<int>& candidateBonds);
    int numBonds, numAtomsPerBond;
    std::vector<int>* bondAtoms;
    ThreadPool* threads;
    std::vector<std::vector<int> > threadBonds;
    std::vector<int> extraBonds;
};

} // namespace OpenMM

#endif /*OPENMM_CPUBONDFORCE_H_*/

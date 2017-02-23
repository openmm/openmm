
/* Portions copyright (c) 2013-2017 Stanford University and Simbios.
 * Authors: Peter Eastman
 * Contributors: 
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

#ifndef __CPU_LANGEVIN_DYNAMICS_H__
#define __CPU_LANGEVIN_DYNAMICS_H__

#include "ReferenceStochasticDynamics.h"
#include "CpuRandom.h"
#include "openmm/internal/ThreadPool.h"
#include "sfmt/SFMT.h"

namespace OpenMM {

class CpuLangevinDynamics : public ReferenceStochasticDynamics {
public:
    /**
     * Constructor.
     *
     * @param numberOfAtoms  number of atoms
     * @param deltaT         delta t for dynamics
     * @param friction       friction coefficient
     * @param temperature    temperature
     * @param threads        thread pool for parallelizing computation
     * @param random         random number generator
     */
    CpuLangevinDynamics(int numberOfAtoms, double deltaT, double friction, double temperature, OpenMM::ThreadPool& threads, OpenMM::CpuRandom& random);

    /**
     * Destructor.
     */
    ~CpuLangevinDynamics();

    /**
     * First update step.
     * 
     * @param numberOfAtoms       number of atoms
     * @param atomCoordinates     atom coordinates
     * @param velocities          velocities
     * @param forces              forces
     * @param inverseMasses       inverse atom masses
     * @param xPrime              xPrime
     */
    void updatePart1(int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& velocities,
                     std::vector<OpenMM::Vec3>& forces, std::vector<double>& inverseMasses, std::vector<OpenMM::Vec3>& xPrime);
      
    /**
     * Second update step.
     * 
     * @param numberOfAtoms       number of atoms
     * @param atomCoordinates     atom coordinates
     * @param velocities          velocities
     * @param forces              forces
     * @param inverseMasses       inverse atom masses
     * @param xPrime              xPrime
     */
    void updatePart2(int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& velocities,
                     std::vector<OpenMM::Vec3>& forces, std::vector<double>& inverseMasses, std::vector<OpenMM::Vec3>& xPrime);

    /**  
     * Third update
     * 
     * @param numberOfAtoms       number of atoms
     * @param atomCoordinates     atom coordinates
     * @param velocities          velocities
     * @param inverseMasses       inverse atom masses
     */
    void updatePart3(int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& velocities,
                     std::vector<double>& inverseMasses, std::vector<OpenMM::Vec3>& xPrime);

private:
    void threadUpdate1(int threadIndex);
    void threadUpdate2(int threadIndex);
    void threadUpdate3(int threadIndex);
    OpenMM::ThreadPool& threads;
    OpenMM::CpuRandom& random;
    std::vector<OpenMM_SFMT::SFMT> threadRandom;
    // The following variables are used to make information accessible to the individual threads.
    int numberOfAtoms;
    OpenMM::Vec3* atomCoordinates;
    OpenMM::Vec3* velocities;
    OpenMM::Vec3* forces;
    double* inverseMasses;
    OpenMM::Vec3* xPrime;
};

} // namespace OpenMM

#endif // __CPU_LANGEVIN_DYNAMICS_H__

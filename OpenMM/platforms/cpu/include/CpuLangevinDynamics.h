
/* Portions copyright (c) 2013-2016 Stanford University and Simbios.
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
    class Update1Task;
    class Update2Task;
    class Update3Task;
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
    CpuLangevinDynamics(int numberOfAtoms, RealOpenMM deltaT, RealOpenMM friction, RealOpenMM temperature, OpenMM::ThreadPool& threads, OpenMM::CpuRandom& random);

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
    void updatePart1(int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates, std::vector<OpenMM::RealVec>& velocities,
                     std::vector<OpenMM::RealVec>& forces, std::vector<RealOpenMM>& inverseMasses, std::vector<OpenMM::RealVec>& xPrime);
      
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
    void updatePart2(int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates, std::vector<OpenMM::RealVec>& velocities,
                     std::vector<OpenMM::RealVec>& forces, std::vector<RealOpenMM>& inverseMasses, std::vector<OpenMM::RealVec>& xPrime);

    /**  
     * Third update
     * 
     * @param numberOfAtoms       number of atoms
     * @param atomCoordinates     atom coordinates
     * @param velocities          velocities
     * @param inverseMasses       inverse atom masses
     */
    void updatePart3(int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates, std::vector<OpenMM::RealVec>& velocities,
                     std::vector<RealOpenMM>& inverseMasses, std::vector<OpenMM::RealVec>& xPrime);

private:
    void threadUpdate1(int threadIndex);
    void threadUpdate2(int threadIndex);
    void threadUpdate3(int threadIndex);
    OpenMM::ThreadPool& threads;
    OpenMM::CpuRandom& random;
    std::vector<OpenMM_SFMT::SFMT> threadRandom;
    // The following variables are used to make information accessible to the individual threads.
    int numberOfAtoms;
    OpenMM::RealVec* atomCoordinates;
    OpenMM::RealVec* velocities;
    OpenMM::RealVec* forces;
    RealOpenMM* inverseMasses;
    OpenMM::RealVec* xPrime;
};

} // namespace OpenMM

#endif // __CPU_LANGEVIN_DYNAMICS_H__

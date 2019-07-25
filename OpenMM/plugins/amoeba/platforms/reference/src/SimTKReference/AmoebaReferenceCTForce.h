
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
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

#ifndef __AmoebaReferenceCTForce_H__
#define __AmoebaReferenceCTForce_H__

#include "RealVec.h"
#include "openmm/Vec3.h"
#include "ReferenceNeighborList.h"
#include <string>
#include <vector>

namespace OpenMM {

class AmoebaReferenceCTForce {

public:
    /**
     * This is an enumeration of the different methods that may be used for handling long range CT forces.
     */
    enum NonbondedMethod {

        /**
         * No cutoff is applied to the interactions.  The full set of N^2 interactions is computed exactly.
         * This necessarily means that periodic boundary conditions cannot be used.  This is the default.
         */

        NoCutoff = 0,

        /**
         * Interactions beyond the cutoff distance are ignored.
         */
        CutoffNonPeriodic = 1,
        /**
         * Periodic boundary conditions are used, so that each particle interacts only with the nearest periodic copy of
         * each other particle.  Interactions beyond the cutoff distance are ignored.
         */
        CutoffPeriodic = 2,
    };


    /**---------------------------------------------------------------------------------------

       Constructor

       --------------------------------------------------------------------------------------- */

    AmoebaReferenceCTForce();

    /**---------------------------------------------------------------------------------------

       Destructor

       --------------------------------------------------------------------------------------- */

    ~AmoebaReferenceCTForce() {};

    /**---------------------------------------------------------------------------------------

       Get nonbonded method

       @return nonbonded method

       --------------------------------------------------------------------------------------- */

    NonbondedMethod getNonbondedMethod() const;

    /**---------------------------------------------------------------------------------------

       Set nonbonded method

       @param nonbonded method

       --------------------------------------------------------------------------------------- */

    void setNonbondedMethod(NonbondedMethod nonbondedMethod);

    /**---------------------------------------------------------------------------------------

       Get cutoff

       @return cutoff

       --------------------------------------------------------------------------------------- */

    double getCutoff() const;

    /**---------------------------------------------------------------------------------------

       Set cutoff

       @param cutoff

       --------------------------------------------------------------------------------------- */

    void setCutoff(double cutoff);

    /**---------------------------------------------------------------------------------------

       Set box dimensions

       @param vectors    the vectors defining the periodic box

       --------------------------------------------------------------------------------------- */

    void setPeriodicBox(OpenMM::RealVec* vectors);

    RealOpenMM calculateForceAndEnergy(int numParticles, int numCTprTypes,
        const std::vector<int>& CTprTypes,
        const std::vector<OpenMM::RealVec>& particlePositions,
        const std::vector<RealOpenMM>& combinedApres, const std::vector<RealOpenMM>& combinedBexps,
        const std::vector<std::set<int> >& CTExclusions,
        std::vector<OpenMM::RealVec>& forces) const;

    RealOpenMM calculateForceAndEnergy(int numParticles, int numCTprTypes,
        const std::vector<int>& CTprTypes,
        const std::vector<OpenMM::RealVec>& particlePositions,
        const std::vector<RealOpenMM>& combinedApres, const std::vector<RealOpenMM>& combinedBexps,
        const NeighborList& neighborList,
        std::vector<OpenMM::RealVec>& forces) const;

private:
    // taper coefficient indices

    static const int C3 = 0;
    static const int C4 = 1;
    static const int C5 = 2;

    NonbondedMethod _nonbondedMethod;
    double _cutoff;
    double _taperCutoffFactor;
    double _taperCutoff;
    RealOpenMM _taperCoefficients[3];
    RealVec _periodicBoxVectors[3];

    /**---------------------------------------------------------------------------------------

       Set taper coefficients

       @param  cutoff cutoff

       --------------------------------------------------------------------------------------- */

    void setTaperCoefficients(double cutoff);

    /**---------------------------------------------------------------------------------------

       Calculate pair ixn

       @param  combindedApre       combined sigmas
       @param  combindedBexp       combined epsilons
       @param  particleIPosition   particle I position
       @param  particleJPosition   particle J position
       @param  force               output force

       @return energy for ixn

       --------------------------------------------------------------------------------------- */

    RealOpenMM calculatePairIxn(RealOpenMM combindedApre, RealOpenMM combindedBexp,
        const Vec3& particleIPosition, const Vec3& particleJPosition,
        Vec3& force) const;
};
}
// ---------------------------------------------------------------------------------------

#endif // _AmoebaReferenceCTForce___

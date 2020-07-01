
/* Portions copyright (c) 2006-2020 Stanford University and Simbios.
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

#ifndef __AmoebaReferenceVdwForce_H__
#define __AmoebaReferenceVdwForce_H__

#include "openmm/Vec3.h"
#include "openmm/AmoebaVdwForce.h"
#include "ReferenceNeighborList.h"
#include <set>
#include <string>
#include <vector>

namespace OpenMM {

class AmoebaReferenceVdwForce;
typedef double (AmoebaReferenceVdwForce::*CombiningFunction)(double x, double y) const;
typedef double (AmoebaReferenceVdwForce::*CombiningFunctionEpsilon)(double x, double y, double z, double w) const;

// ---------------------------------------------------------------------------------------

class AmoebaReferenceVdwForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    AmoebaReferenceVdwForce();
    
    void initialize(const AmoebaVdwForce& force);
 
    /**---------------------------------------------------------------------------------------
    
       Set cutoff
    
       @param cutoff
    
       --------------------------------------------------------------------------------------- */
    
    void setCutoff(double cutoff);

    /**---------------------------------------------------------------------------------------
    
       Set box dimensions
    
       @param vectors    the vectors defining the periodic box
    
       --------------------------------------------------------------------------------------- */
    
    void setPeriodicBox(OpenMM::Vec3* vectors);
 
    /**---------------------------------------------------------------------------------------
    
       Get the set of exclusions for each particle.
    
       --------------------------------------------------------------------------------------- */
    
    std::vector<std::set<int> >& getExclusions();

    /**---------------------------------------------------------------------------------------
    
       Calculate Amoeba Hal vdw ixns
    
       @param numParticles            number of particles
       @param lambda                  lambda value
       @param particlePositions       Cartesian coordinates of particles
       @param forces                  add forces to this vector
    
       @return energy
    
       --------------------------------------------------------------------------------------- */
    
    double calculateForceAndEnergy(int numParticles, double lambda, const std::vector<OpenMM::Vec3>& particlePositions,
                                   std::vector<OpenMM::Vec3>& forces) const;
         
    /**---------------------------------------------------------------------------------------
    
       Calculate Vdw ixn using neighbor list
    
       @param numParticles            number of particles
       @param lambda                  lambda value
       @param particlePositions       Cartesian coordinates of particles
       @param neighborList            neighbor list
       @param forces                  add forces to this vector
    
       @return energy
    
       --------------------------------------------------------------------------------------- */
    
    double calculateForceAndEnergy(int numParticles, double lambda, const std::vector<OpenMM::Vec3>& particlePositions, 
                                   const NeighborList& neighborList, std::vector<OpenMM::Vec3>& forces) const;
         
private:
    // taper coefficient indices
    static const int C3=0;
    static const int C4=1;
    static const int C5=2;

    AmoebaVdwForce::NonbondedMethod _nonbondedMethod;
    AmoebaVdwForce::AlchemicalMethod _alchemicalMethod;
    AmoebaVdwForce::PotentialFunction potentialFunction;
    double _n;
    double _alpha;
    double _cutoff;
    double _taperCutoffFactor;
    double _taperCutoff;
    double _taperCoefficients[3];
    std::vector<int> particleType;
    std::vector<std::vector<double> > sigmaMatrix;
    std::vector<std::vector<double> > epsilonMatrix;
    std::vector<int> indexIVs;
    std::vector<double> reductions;
    std::vector<bool> isAlchemical;
    std::vector<std::set<int> > allExclusions;
    Vec3 _periodicBoxVectors[3];

    /**---------------------------------------------------------------------------------------
    
       Set reduced positions: position used to calculate vdw interaction is moved towards 
       covalent partner
       
    
       @param  numParticles         number of particles
       @param  particlePositions    current particle positions
       @param  indexIVs             particle index of covalent partner
       @param  reductions           fraction of bond length to move particle interacting site;
                                    reductions[i] = zero, 
                                    if interacting position == particle position
       @param  reducedPositions     output: modfied or original position depending on whether
                                    reduction factor is nonzero
    
       --------------------------------------------------------------------------------------- */
    
    void setReducedPositions(int numParticles, const std::vector<Vec3>& particlePositions,
                             const std::vector<int>& indexIVs, const std::vector<double>& reductions,
                             std::vector<Vec3>& reducedPositions) const;

    /**---------------------------------------------------------------------------------------
    
       Add reduced forces to force vector
    
       @param  particleI            index of particleI
       @param  particleIV           index of particleIV
       @param  reduction            reduction factor
       @param  sign                 +1 or -1 for add/sutracting forces
       @param  force                force vector to add
       @param  forces               force vector for particles
    
       --------------------------------------------------------------------------------------- */
    
    void addReducedForce(unsigned int particleI, unsigned int particleIV,
                         double reduction, double sign,
                         Vec3& force, std::vector<OpenMM::Vec3>& forces) const;
    
    /**---------------------------------------------------------------------------------------
    
       Set taper coefficients
    
       @param  cutoff cutoff

       --------------------------------------------------------------------------------------- */
    
    void setTaperCoefficients(double cutoff);

    /**---------------------------------------------------------------------------------------
    
       Calculate pair ixn
    
       @param  combindedSigma       combined sigmas
       @param  combindedEpsilon     combined epsilons
       @param  softcore             softcore offset parameter
       @param  particleIPosition    particle I position 
       @param  particleJPosition    particle J position 
       @param  force                output force
    
       @return energy for ixn

       --------------------------------------------------------------------------------------- */
    
    double calculatePairIxn(double combindedSigma, double combindedEpsilon, double softcore,
                            const Vec3& particleIPosition, const Vec3& particleJPosition,
                            Vec3& force) const;

};

}
// ---------------------------------------------------------------------------------------

#endif // _AmoebaReferenceVdwForce___

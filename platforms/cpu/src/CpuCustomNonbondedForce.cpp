
/* Portions copyright (c) 2009-2014 Stanford University and Simbios.
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

#include <string.h>
#include <sstream>

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"
#include "CpuCustomNonbondedForce.h"
#include "gmx_atomic.h"

using namespace OpenMM;
using namespace std;

/**---------------------------------------------------------------------------------------

   CpuCustomNonbondedForce constructor

   --------------------------------------------------------------------------------------- */

CpuCustomNonbondedForce::CpuCustomNonbondedForce(const Lepton::CompiledExpression& energyExpression,
        const Lepton::CompiledExpression& forceExpression, const vector<string>& parameterNames) :
            cutoff(false), useSwitch(false), periodic(false), energyExpression(energyExpression), forceExpression(forceExpression), paramNames(parameterNames) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuCustomNonbondedForce::CpuCustomNonbondedForce";

   // ---------------------------------------------------------------------------------------

    energyR = ReferenceForce::getVariablePointer(this->energyExpression, "r");
    forceR = ReferenceForce::getVariablePointer(this->forceExpression, "r");
    for (int i = 0; i < (int) paramNames.size(); i++) {
        for (int j = 1; j < 3; j++) {
            stringstream name;
            name << paramNames[i] << j;
            energyParticleParams.push_back(ReferenceForce::getVariablePointer(this->energyExpression, name.str()));
            forceParticleParams.push_back(ReferenceForce::getVariablePointer(this->forceExpression, name.str()));
        }
    }
}

/**---------------------------------------------------------------------------------------

   CpuCustomNonbondedForce destructor

   --------------------------------------------------------------------------------------- */

CpuCustomNonbondedForce::~CpuCustomNonbondedForce( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuCustomNonbondedForce::~CpuCustomNonbondedForce";

   // ---------------------------------------------------------------------------------------

}

  /**---------------------------------------------------------------------------------------

     Set the force to use a cutoff.

     @param distance            the cutoff distance
     @param neighbors           the neighbor list to use

     --------------------------------------------------------------------------------------- */

  void CpuCustomNonbondedForce::setUseCutoff( RealOpenMM distance, const CpuNeighborList& neighbors ) {

    cutoff = true;
    cutoffDistance = distance;
    neighborList = &neighbors;
  }

/**---------------------------------------------------------------------------------------

   Restrict the force to a list of interaction groups.

   @param distance            the cutoff distance
   @param neighbors           the neighbor list to use

   --------------------------------------------------------------------------------------- */

void CpuCustomNonbondedForce::setInteractionGroups(const vector<pair<set<int>, set<int> > >& groups) {
    interactionGroups = groups;
}

/**---------------------------------------------------------------------------------------

   Set the force to use a switching function.

   @param distance            the switching distance

   --------------------------------------------------------------------------------------- */

void CpuCustomNonbondedForce::setUseSwitchingFunction( RealOpenMM distance ) {
    useSwitch = true;
    switchingDistance = distance;
}

  /**---------------------------------------------------------------------------------------

     Set the force to use periodic boundary conditions.  This requires that a cutoff has
     also been set, and the smallest side of the periodic box is at least twice the cutoff
     distance.

     @param boxSize             the X, Y, and Z widths of the periodic box

     --------------------------------------------------------------------------------------- */

  void CpuCustomNonbondedForce::setPeriodic( RealVec& boxSize ) {

    assert(cutoff);
    assert(boxSize[0] >= 2.0*cutoffDistance);
    assert(boxSize[1] >= 2.0*cutoffDistance);
    assert(boxSize[2] >= 2.0*cutoffDistance);
    periodic = true;
    periodicBoxSize[0] = boxSize[0];
    periodicBoxSize[1] = boxSize[1];
    periodicBoxSize[2] = boxSize[2];

  }


void CpuCustomNonbondedForce::calculatePairIxn(int numberOfAtoms, float* posq, vector<RealVec>& atomCoordinates,
                                             RealOpenMM** atomParameters, vector<set<int> >& exclusions,
                                             RealOpenMM* fixedParameters, const map<string, double>& globalParameters,
                                             vector<AlignedArray<float> >& threadForce, double* totalEnergy, ThreadPool& threads) {
    // Record the parameters for the threads.
    
    this->numberOfAtoms = numberOfAtoms;
    this->posq = posq;
    this->atomCoordinates = &atomCoordinates[0];
    this->atomParameters = atomParameters;
    this->exclusions = &exclusions[0];
    this->threadForce = &threadForce;
    includeEnergy = (totalEnergy != NULL);
    threadEnergy.resize(threads.getNumThreads());
    gmx_atomic_t counter;
    gmx_atomic_set(&counter, 0);
    this->atomicCounter = &counter;

    
    float* forces = &threadForce[0][0];
    
    
    for (map<string, double>::const_iterator iter = globalParameters.begin(); iter != globalParameters.end(); ++iter) {
        ReferenceForce::setVariable(ReferenceForce::getVariablePointer(energyExpression, iter->first), iter->second);
        ReferenceForce::setVariable(ReferenceForce::getVariablePointer(forceExpression, iter->first), iter->second);
    }
    fvec4 boxSize(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2], 0);
    fvec4 invBoxSize((1/periodicBoxSize[0]), (1/periodicBoxSize[1]), (1/periodicBoxSize[2]), 0);
    if (interactionGroups.size() > 0) {
        // The user has specified interaction groups, so compute only the requested interactions.
        
        for (int group = 0; group < (int) interactionGroups.size(); group++) {
            const set<int>& set1 = interactionGroups[group].first;
            const set<int>& set2 = interactionGroups[group].second;
            for (set<int>::const_iterator atom1 = set1.begin(); atom1 != set1.end(); ++atom1) {
                for (set<int>::const_iterator atom2 = set2.begin(); atom2 != set2.end(); ++atom2) {
                    if (*atom1 == *atom2 || exclusions[*atom1].find(*atom2) != exclusions[*atom1].end())
                        continue; // This is an excluded interaction.
                    if (*atom1 > *atom2 && set1.find(*atom2) != set1.end() && set2.find(*atom1) != set2.end())
                        continue; // Both atoms are in both sets, so skip duplicate interactions.
                    for (int j = 0; j < (int) paramNames.size(); j++) {
                        ReferenceForce::setVariable(energyParticleParams[j*2], atomParameters[*atom1][j]);
                        ReferenceForce::setVariable(energyParticleParams[j*2+1], atomParameters[*atom2][j]);
                        ReferenceForce::setVariable(forceParticleParams[j*2], atomParameters[*atom1][j]);
                        ReferenceForce::setVariable(forceParticleParams[j*2+1], atomParameters[*atom2][j]);
                    }
                    calculateOneIxn(*atom1, *atom2, forces, totalEnergy, boxSize, invBoxSize);
                }
            }
        }
    }
    else if (cutoff) {
        // We are using a cutoff, so get the interactions from the neighbor list.
        
        for (int blockIndex = 0; blockIndex < neighborList->getNumBlocks(); blockIndex++) {
            int blockAtom[4];
            for (int i = 0; i < 4; i++)
                blockAtom[i] = neighborList->getSortedAtoms()[4*blockIndex+i];
            const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
            const vector<char>& exclusions = neighborList->getBlockExclusions(blockIndex);
            for (int i = 0; i < (int) neighbors.size(); i++) {
                int first = neighbors[i];
                for (int j = 0; j < (int) paramNames.size(); j++) {
                    ReferenceForce::setVariable(energyParticleParams[j*2], atomParameters[first][j]);
                    ReferenceForce::setVariable(forceParticleParams[j*2], atomParameters[first][j]);
                }
                for (int k = 0; k < 4; k++) {
                    if ((exclusions[i] & (1<<k)) == 0) {
                        int second = blockAtom[k];
                        for (int j = 0; j < (int) paramNames.size(); j++) {
                            ReferenceForce::setVariable(energyParticleParams[j*2+1], atomParameters[second][j]);
                            ReferenceForce::setVariable(forceParticleParams[j*2+1], atomParameters[second][j]);
                        }
                        calculateOneIxn(first, second, forces, totalEnergy, boxSize, invBoxSize);
                    }
                }
                        
            }
        }
//        for (int i = 0; i < (int) neighborList->size(); i++) {
//            OpenMM::AtomPair pair = (*neighborList)[i];
//            for (int j = 0; j < (int) paramNames.size(); j++) {
//                ReferenceForce::setVariable(energyParticleParams[j*2], atomParameters[pair.first][j]);
//                ReferenceForce::setVariable(energyParticleParams[j*2+1], atomParameters[pair.second][j]);
//                ReferenceForce::setVariable(forceParticleParams[j*2], atomParameters[pair.first][j]);
//                ReferenceForce::setVariable(forceParticleParams[j*2+1], atomParameters[pair.second][j]);
//            }
//            calculateOneIxn(pair.first, pair.second, atomCoordinates, forces, energyByAtom, totalEnergy);
//        }
    }
    else {
        // Every particle interacts with every other one.
        
        for (int ii = 0; ii < numberOfAtoms; ii++) {
            for (int jj = ii+1; jj < numberOfAtoms; jj++) {
                if (exclusions[jj].find(ii) == exclusions[jj].end()) {
                    for (int j = 0; j < (int) paramNames.size(); j++) {
                        ReferenceForce::setVariable(energyParticleParams[j*2], atomParameters[ii][j]);
                        ReferenceForce::setVariable(energyParticleParams[j*2+1], atomParameters[jj][j]);
                        ReferenceForce::setVariable(forceParticleParams[j*2], atomParameters[ii][j]);
                        ReferenceForce::setVariable(forceParticleParams[j*2+1], atomParameters[jj][j]);
                    }
                    calculateOneIxn(ii, jj, forces, totalEnergy, boxSize, invBoxSize);
                }
            }
        }
    }
}

  /**---------------------------------------------------------------------------------------

     Calculate one pair ixn between two atoms

     @param ii               the index of the first atom
     @param jj               the index of the second atom
     @param atomCoordinates  atom coordinates
     @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
     @param forces           force array (forces added)
     @param totalEnergy      total energy

     --------------------------------------------------------------------------------------- */

void CpuCustomNonbondedForce::calculateOneIxn(int ii, int jj, float* forces, RealOpenMM* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {

    // ---------------------------------------------------------------------------------------

    static const std::string methodName = "\nCpuCustomNonbondedForce::calculateOneIxn";

    // ---------------------------------------------------------------------------------------

    // constants -- reduce Visual Studio warnings regarding conversions between float & double

    static const RealOpenMM zero        =  0.0;
    static const RealOpenMM one         =  1.0;
    static const RealOpenMM two         =  2.0;
    static const RealOpenMM three       =  3.0;
    static const RealOpenMM six         =  6.0;
    static const RealOpenMM twelve      = 12.0;
    static const RealOpenMM oneM        = -1.0;

    // get deltaR, R2, and R between 2 atoms

    fvec4 deltaR;
    fvec4 posI(posq+4*ii);
    fvec4 posJ(posq+4*jj);
    float r2;
    getDeltaR(posI, posJ, deltaR, r2, boxSize, invBoxSize);
    if (cutoff && r2 >= cutoffDistance*cutoffDistance)
        return;
    float r = sqrtf(r2);

    // accumulate forces

    ReferenceForce::setVariable(energyR, r);
    ReferenceForce::setVariable(forceR, r);
    double dEdR = forceExpression.evaluate()/r;
    double energy = energyExpression.evaluate();
    if (useSwitch) {
        if (r > switchingDistance) {
            RealOpenMM t = (r-switchingDistance)/(cutoffDistance-switchingDistance);
            RealOpenMM switchValue = 1+t*t*t*(-10+t*(15-t*6));
            RealOpenMM switchDeriv = t*t*(-30+t*(60-t*30))/(cutoffDistance-switchingDistance);
            dEdR = switchValue*dEdR + energy*switchDeriv/r;
            energy *= switchValue;
        }
    }
    fvec4 result = deltaR*dEdR;
    (fvec4(forces+4*ii)+result).store(forces+4*ii);
    (fvec4(forces+4*jj)-result).store(forces+4*jj);

    // accumulate energies

    if (totalEnergy)
        *totalEnergy += energy;
}

void CpuCustomNonbondedForce::getDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, const fvec4& boxSize, const fvec4& invBoxSize) const {
    deltaR = posJ-posI;
    if (periodic) {
        fvec4 base = round(deltaR*invBoxSize)*boxSize;
        deltaR = deltaR-base;
    }
    r2 = dot3(deltaR, deltaR);
}


/* Portions copyright (c) 2006-2013 Stanford University and Simbios.
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

#include <string.h>
#include <sstream>
#include <complex>
#include <algorithm>
#include <iostream>

#include "SimTKOpenMMUtilities.h"
#include "ReferenceDPMECoulombIxn.h"
#include "ReferenceForce.h"
#include "ReferencePME.h"
#include "openmm/OpenMMException.h"
#include "openmm/Units.h"

// In case we're using some primitive version of Visual Studio this will
// make sure that erf() and erfc() are defined.
#include "openmm/internal/MSVC_erfc.h"

using std::set;
using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceDPMECoulombIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceDPMECoulombIxn::ReferenceDPMECoulombIxn() : cutoff(false), periodic(false), pme_setup(false), dpme_setup(false) {
}

/**---------------------------------------------------------------------------------------

   ReferenceDPMECoulombIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceDPMECoulombIxn::~ReferenceDPMECoulombIxn() {
}

/**---------------------------------------------------------------------------------------

     Set the force to use a cutoff.

     @param distance            the cutoff distance
     @param neighbors           the neighbor list to use
     @param solventDielectric   the dielectric constant of the bulk solvent

     --------------------------------------------------------------------------------------- */

void ReferenceDPMECoulombIxn::setUseCutoff(RealOpenMM distance, const OpenMM::NeighborList& neighbors, RealOpenMM solventDielectric) {

    cutoff = true;
    cutoffDistance = distance;
    neighborList = &neighbors;
    krf = pow(cutoffDistance, -3.0)*(solventDielectric-1.0)/(2.0*solventDielectric+1.0);
    crf = (1.0/cutoffDistance)*(3.0*solventDielectric)/(2.0*solventDielectric+1.0);
}


/**---------------------------------------------------------------------------------------

     Set the force to use periodic boundary conditions.  This requires that a cutoff has
     also been set, and the smallest side of the periodic box is at least twice the cutoff
     distance.

     @param vectors    the vectors defining the periodic box

     --------------------------------------------------------------------------------------- */

void ReferenceDPMECoulombIxn::setPeriodic(OpenMM::RealVec* vectors) {

    assert(cutoff);
    assert(vectors[0][0] >= 2.0*cutoffDistance);
    assert(vectors[1][1] >= 2.0*cutoffDistance);
    assert(vectors[2][2] >= 2.0*cutoffDistance);
    periodic = true;
    periodicBoxVectors[0] = vectors[0];
    periodicBoxVectors[1] = vectors[1];
    periodicBoxVectors[2] = vectors[2];
}


/**---------------------------------------------------------------------------------------

     Set the force to use Particle-Mesh Ewald (PME) summation.

     @param alpha  the Ewald separation parameter
     @param gridSize the dimensions of the mesh

     --------------------------------------------------------------------------------------- */

void ReferenceDPMECoulombIxn::setupPME(RealOpenMM alpha, int meshSize[3]) {
    alphaEwald = alpha;
    meshDim[0] = meshSize[0];
    meshDim[1] = meshSize[1];
    meshDim[2] = meshSize[2];
    pme_setup = true;
}

/**---------------------------------------------------------------------------------------

     Set the force to use Particle-Mesh Ewald (PME) summation for dispersion.

     @param alpha  the Ewald dispersion separation parameter
     @param gridSize the dimensions of the dispersion mesh

     --------------------------------------------------------------------------------------- */

void ReferenceDPMECoulombIxn::setupDispersionPME(RealOpenMM dalpha, int dispersionMeshSize[3]) {
    alphaDispersionEwald = dalpha;
    dispersionMeshDim[0] = dispersionMeshSize[0];
    dispersionMeshDim[1] = dispersionMeshSize[1];
    dispersionMeshDim[2] = dispersionMeshSize[2];
    dpme_setup = true;
}

/**---------------------------------------------------------------------------------------

   Calculate PME ixn

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param atomParameters   atom parameters                             atomParameters[atomIndex][paramterIndex]
   @param exclusions       atom exclusion indices
                           exclusions[atomIndex] contains the list of exclusions for that atom
   @param fixedParameters  non atom parameters (not currently used)
   @param forces           force array (forces added)
   @param energyByAtom     atom energy
   @param totalEnergy      total energy
   @param includeDirect      true if direct space interactions should be included
   @param includeReciprocal  true if reciprocal space interactions should be included

   --------------------------------------------------------------------------------------- */

void ReferenceDPMECoulombIxn::calculatePMEIxn(int numberOfAtoms, vector<RealVec>& atomCoordinates,
                                              RealOpenMM** atomParameters, vector<set<int> >& exclusions,
                                              RealOpenMM* fixedParameters, vector<RealVec>& forces,
                                              RealOpenMM* energyByAtom, RealOpenMM* totalEnergy, bool includeDirect, bool includeReciprocal) const {
    typedef std::complex<RealOpenMM> d_complex;

    static const RealOpenMM one         =  1.0;
    static const RealOpenMM six         =  6.0;
    static const RealOpenMM twelve      = 12.0;

    RealOpenMM SQRT_PI                  = sqrt(PI_M);

    RealOpenMM totalSelfEwaldEnergy     = 0.0;
    RealOpenMM realSpaceEwaldEnergy     = 0.0;
    RealOpenMM recipEnergy              = 0.0;
    RealOpenMM vdwEnergy                = 0.0;

    // **************************************************************************************
    // SELF ENERGY
    // **************************************************************************************

    if(!pme_setup)
        throw OpenMMException("setupPME must be called before computing energies or forces.");
    if(!dpme_setup)
        throw OpenMMException("setupDispersionPME must be called before computing energies or forces.");
    if (includeReciprocal) {

        for (int atomID = 0; atomID < numberOfAtoms; atomID++) {
            // Electrostatic term
            RealOpenMM selfEwaldEnergy       = (RealOpenMM) (ONE_4PI_EPS0*atomParameters[atomID][QIndex]*atomParameters[atomID][QIndex] * alphaEwald/SQRT_PI);
            // Dispersion term
            selfEwaldEnergy -= pow(alphaDispersionEwald, 6.0) * 64.0*pow(atomParameters[atomID][SigIndex], 6.0) * pow(atomParameters[atomID][EpsIndex], 2.0) / 12.0;

            totalSelfEwaldEnergy            -= selfEwaldEnergy;
            if (energyByAtom) {
                energyByAtom[atomID]        -= selfEwaldEnergy;
            }
        }
    }

    if (totalEnergy) {
        *totalEnergy += totalSelfEwaldEnergy;
    }

    // **************************************************************************************
    // RECIPROCAL SPACE EWALD ENERGY AND FORCES
    // **************************************************************************************
    // PME

    if (includeReciprocal) {
        pme_t          pmedata; /* abstract handle for PME data */

        pme_init(&pmedata,alphaEwald,numberOfAtoms,meshDim,5,1);
        vector<RealOpenMM> charges(numberOfAtoms);
        for (int i = 0; i < numberOfAtoms; i++)
            charges[i] = atomParameters[i][QIndex];
        pme_exec(pmedata,atomCoordinates,forces,charges,periodicBoxVectors,&recipEnergy);
        pme_destroy(pmedata);
        if (totalEnergy)
            *totalEnergy += recipEnergy;

        if (energyByAtom)
            for (int n = 0; n < numberOfAtoms; n++)
                energyByAtom[n] += recipEnergy;

        // Dispersion reciprocal space terms
        pme_init(&pmedata,alphaDispersionEwald,numberOfAtoms,dispersionMeshDim,5,1);
        std::vector<RealVec> dpmeforces;

        for (int i = 0; i < numberOfAtoms; i++){
            charges[i] = 8.0*pow(atomParameters[i][SigIndex], 3.0) * atomParameters[i][EpsIndex];
            dpmeforces.push_back(RealVec());
        }
        pme_exec_dpme(pmedata,atomCoordinates,dpmeforces,charges,periodicBoxVectors,&recipEnergy);
        for (int i = 0; i < numberOfAtoms; i++){
            forces[i][0] -= 2.0*dpmeforces[i][0];
            forces[i][1] -= 2.0*dpmeforces[i][1];
            forces[i][2] -= 2.0*dpmeforces[i][2];
        }

        pme_destroy(pmedata);

        if (totalEnergy)
            *totalEnergy += recipEnergy;

        if (energyByAtom)
            for (int n = 0; n < numberOfAtoms; n++)
                energyByAtom[n] += recipEnergy;


    }



    // **************************************************************************************
    // SHORT-RANGE ENERGY AND FORCES
    // **************************************************************************************

    if (!includeDirect)
        return;
    RealOpenMM totalVdwEnergy            = 0.0f;
    RealOpenMM totalRealSpaceEwaldEnergy = 0.0f;

    for (int i = 0; i < (int) neighborList->size(); i++) {
        OpenMM::AtomPair pair = (*neighborList)[i];
        int ii = pair.first;
        int jj = pair.second;

        RealOpenMM deltaR[2][ReferenceForce::LastDeltaRIndex];
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[jj], atomCoordinates[ii], periodicBoxVectors, deltaR[0]);
        RealOpenMM r         = deltaR[0][ReferenceForce::RIndex];
        RealOpenMM inverseR  = one/(deltaR[0][ReferenceForce::RIndex]);
        RealOpenMM alphaR    = alphaEwald * r;
        RealOpenMM dalphaR   = alphaDispersionEwald * r;
        RealOpenMM dar2 = dalphaR*dalphaR;
        RealOpenMM dar4 = dar2*dar2;
        RealOpenMM dar6 = dar4*dar2;
        RealOpenMM inverseR2 = inverseR*inverseR;



        RealOpenMM dEdR      = (RealOpenMM) (ONE_4PI_EPS0 * atomParameters[ii][QIndex] * atomParameters[jj][QIndex] * inverseR * inverseR * inverseR);
        dEdR      = (RealOpenMM) (dEdR * (erfc(alphaR) + 2 * alphaR * exp (- alphaR * alphaR) / SQRT_PI));

        RealOpenMM sig       = atomParameters[ii][SigIndex] +  atomParameters[jj][SigIndex];
        RealOpenMM sig2      = inverseR*sig;
        sig2     *= sig2;
        RealOpenMM sig6      = sig2*sig2*sig2;
        RealOpenMM eps       = atomParameters[ii][EpsIndex]*atomParameters[jj][EpsIndex];
        RealOpenMM c6i = 8.0*pow(atomParameters[ii][SigIndex], 3.0) * atomParameters[ii][EpsIndex];
        RealOpenMM c6j = 8.0*pow(atomParameters[jj][SigIndex], 3.0) * atomParameters[jj][EpsIndex];
        // For the energies and forces, we first add the regular Lorentz−Berthelot terms.  The C12 term is treated as usual
        // but we then subtract out (remembering that the C6 term is negative) the multiplicative C6 term that has been
        // computed in real space.  Finally, we add a potential shift term to account for the difference between the LB
        // and multiplicative functional forms at the cutoff.
        dEdR     += eps*(twelve*sig6 - six)*sig6*inverseR*inverseR;
        vdwEnergy = eps*(sig6-one)*sig6;
        RealOpenMM emult = c6i*c6j*inverseR2*inverseR2*inverseR2*(1.0 - EXP(-dar2) * (1.0 + dar2 + 0.5*dar4));
        dEdR += 6.0*c6i*c6j*inverseR2*inverseR2*inverseR2*inverseR2*(1.0 - EXP(-dar2) * (1.0 + dar2 + 0.5*dar4 + dar6/6.0));
        // The additive part of the potential shift
        RealOpenMM inverseCut2 = 1.0/(cutoffDistance*cutoffDistance);
        RealOpenMM inverseCut6 = inverseCut2*inverseCut2*inverseCut2;
        sig2 = atomParameters[ii][SigIndex] +  atomParameters[jj][SigIndex];
        sig2 *= sig2;
        sig6 = sig2*sig2*sig2;
        RealOpenMM potentialshift = eps*sig6*inverseCut6;
        dalphaR   = alphaDispersionEwald * cutoffDistance;
        dar2 = dalphaR*dalphaR;
        dar4 = dar2*dar2;
        // The multiplicative part of the potential shift
        potentialshift -= c6i*c6j*inverseCut6*(1.0 - EXP(-dar2) * (1.0 + dar2 + 0.5*dar4));
        vdwEnergy += emult + potentialshift;

        // accumulate forces

        for (int kk = 0; kk < 3; kk++) {
            RealOpenMM force  = dEdR*deltaR[0][kk];
            forces[ii][kk]   += force;
            forces[jj][kk]   -= force;
        }

        // accumulate energies

        realSpaceEwaldEnergy        = (RealOpenMM) (ONE_4PI_EPS0*atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*inverseR*erfc(alphaR));

        totalVdwEnergy             += vdwEnergy;
        totalRealSpaceEwaldEnergy  += realSpaceEwaldEnergy;

        if (energyByAtom) {
            energyByAtom[ii] += realSpaceEwaldEnergy + vdwEnergy;
            energyByAtom[jj] += realSpaceEwaldEnergy + vdwEnergy;
        }

    }

    if (totalEnergy)
        *totalEnergy += totalRealSpaceEwaldEnergy + totalVdwEnergy;

    // Now subtract off the exclusions, since they were implicitly included in the reciprocal space sum.

    RealOpenMM totalExclusionEnergy = 0.0f;
    const double TWO_OVER_SQRT_PI = 2/sqrt(PI_M);
    for (int i = 0; i < numberOfAtoms; i++)
        for (set<int>::const_iterator iter = exclusions[i].begin(); iter != exclusions[i].end(); ++iter) {
            if (*iter > i) {
                int ii = i;
                int jj = *iter;

                RealOpenMM deltaR[2][ReferenceForce::LastDeltaRIndex];
                ReferenceForce::getDeltaR(atomCoordinates[jj], atomCoordinates[ii], deltaR[0]);
                RealOpenMM r         = deltaR[0][ReferenceForce::RIndex];
                RealOpenMM inverseR  = one/(deltaR[0][ReferenceForce::RIndex]);
                RealOpenMM alphaR    = alphaEwald * r;
                if (erf(alphaR) > 1e-6) {
                    RealOpenMM dEdR      = (RealOpenMM) (ONE_4PI_EPS0 * atomParameters[ii][QIndex] * atomParameters[jj][QIndex] * inverseR * inverseR * inverseR);
                               dEdR      = (RealOpenMM) (dEdR * (erf(alphaR) - 2 * alphaR * exp (- alphaR * alphaR) / SQRT_PI));

                    // accumulate forces

                    for (int kk = 0; kk < 3; kk++) {
                        RealOpenMM force  = dEdR*deltaR[0][kk];
                        forces[ii][kk]   -= force;
                        forces[jj][kk]   += force;
                    }

                    // accumulate energies

                    realSpaceEwaldEnergy = (RealOpenMM) (ONE_4PI_EPS0*atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*inverseR*erf(alphaR));
                }
                else {
                    realSpaceEwaldEnergy = (RealOpenMM) (alphaEwald*TWO_OVER_SQRT_PI*ONE_4PI_EPS0*atomParameters[ii][QIndex]*atomParameters[jj][QIndex]);
                }

                // Dispersion terms.  Here we just back out the reciprocal space terms, and don't add any extra real space terms.
                RealOpenMM dalphaR   = alphaDispersionEwald * r;
                RealOpenMM inverseR2 = inverseR*inverseR;
                RealOpenMM dar2 = dalphaR*dalphaR;
                RealOpenMM dar4 = dar2*dar2;
                RealOpenMM dar6 = dar4*dar2;
                RealOpenMM c6i = 8.0*pow(atomParameters[ii][SigIndex], 3.0) * atomParameters[ii][EpsIndex];
                RealOpenMM c6j = 8.0*pow(atomParameters[jj][SigIndex], 3.0) * atomParameters[jj][EpsIndex];
                realSpaceEwaldEnergy += c6i*c6j*inverseR2*inverseR2*inverseR2*(1.0 - EXP(-dar2) * (1.0 + dar2 + 0.5*dar4));
                RealOpenMM dEdR = 6.0*c6i*c6j*inverseR2*inverseR2*inverseR2*inverseR2*(1.0 - EXP(-dar2) * (1.0 + dar2 + 0.5*dar4 + dar6/6.0));
                for (int kk = 0; kk < 3; kk++) {
                    RealOpenMM force  = dEdR*deltaR[0][kk];
                    forces[ii][kk]   -= force;
                    forces[jj][kk]   += force;
                }

                totalExclusionEnergy += realSpaceEwaldEnergy;
                if (energyByAtom) {
                    energyByAtom[ii] -= realSpaceEwaldEnergy;
                    energyByAtom[jj] -= realSpaceEwaldEnergy;
                }
            }
        }
    if (totalEnergy)
        *totalEnergy -= totalExclusionEnergy;
}


/**---------------------------------------------------------------------------------------

   Calculate DPME Coulomb pair ixn

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param atomParameters   atom parameters                             atomParameters[atomIndex][paramterIndex]
   @param exclusions       atom exclusion indices
                           exclusions[atomIndex] contains the list of exclusions for that atom
   @param fixedParameters  non atom parameters (not currently used)
   @param forces           force array (forces added)
   @param energyByAtom     atom energy
   @param totalEnergy      total energy
   @param includeDirect      true if direct space interactions should be included
   @param includeReciprocal  true if reciprocal space interactions should be included

   --------------------------------------------------------------------------------------- */

void ReferenceDPMECoulombIxn::calculatePairIxn(int numberOfAtoms, vector<RealVec>& atomCoordinates,
                                               RealOpenMM** atomParameters, vector<set<int> >& exclusions,
                                               RealOpenMM* fixedParameters, vector<RealVec>& forces,
                                               RealOpenMM* energyByAtom, RealOpenMM* totalEnergy, bool includeDirect, bool includeReciprocal) const {

    calculatePMEIxn(numberOfAtoms, atomCoordinates, atomParameters, exclusions, fixedParameters, forces, energyByAtom,
                      totalEnergy, includeDirect, includeReciprocal);
    return;
}


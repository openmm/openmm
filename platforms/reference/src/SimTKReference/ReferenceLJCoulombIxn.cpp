
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

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "ReferenceLJCoulombIxn.h"
#include "ReferenceForce.h"
#include "ReferencePME.h"

// In case we're using some primitive version of Visual Studio this will
// make sure that erf() and erfc() are defined.
#include "openmm/internal/MSVC_erfc.h"

using std::set;
using std::vector;
using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

   ReferenceLJCoulombIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceLJCoulombIxn::ReferenceLJCoulombIxn( ) : cutoff(false), useSwitch(false), periodic(false), ewald(false), pme(false) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceLJCoulombIxn::ReferenceLJCoulombIxn";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   ReferenceLJCoulombIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceLJCoulombIxn::~ReferenceLJCoulombIxn( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceLJCoulombIxn::~ReferenceLJCoulombIxn";

   // ---------------------------------------------------------------------------------------

}

  /**---------------------------------------------------------------------------------------

     Set the force to use a cutoff.

     @param distance            the cutoff distance
     @param neighbors           the neighbor list to use
     @param solventDielectric   the dielectric constant of the bulk solvent

     --------------------------------------------------------------------------------------- */

  void ReferenceLJCoulombIxn::setUseCutoff( RealOpenMM distance, const OpenMM::NeighborList& neighbors, RealOpenMM solventDielectric ) {

    cutoff = true;
    cutoffDistance = distance;
    neighborList = &neighbors;
    krf = pow(cutoffDistance, -3.0)*(solventDielectric-1.0)/(2.0*solventDielectric+1.0);
    crf = (1.0/cutoffDistance)*(3.0*solventDielectric)/(2.0*solventDielectric+1.0);
  }

/**---------------------------------------------------------------------------------------

   Set the force to use a switching function on the Lennard-Jones interaction.

   @param distance            the switching distance

   --------------------------------------------------------------------------------------- */

void ReferenceLJCoulombIxn::setUseSwitchingFunction( RealOpenMM distance ) {
    useSwitch = true;
    switchingDistance = distance;
}

  /**---------------------------------------------------------------------------------------

     Set the force to use periodic boundary conditions.  This requires that a cutoff has
     also been set, and the smallest side of the periodic box is at least twice the cutoff
     distance.

     @param boxSize             the X, Y, and Z widths of the periodic box

     --------------------------------------------------------------------------------------- */

  void ReferenceLJCoulombIxn::setPeriodic( RealVec& boxSize ) {

    assert(cutoff);
    assert(boxSize[0] >= 2.0*cutoffDistance);
    assert(boxSize[1] >= 2.0*cutoffDistance);
    assert(boxSize[2] >= 2.0*cutoffDistance);
    periodic = true;
    periodicBoxSize[0] = boxSize[0];
    periodicBoxSize[1] = boxSize[1];
    periodicBoxSize[2] = boxSize[2];
  }

  /**---------------------------------------------------------------------------------------

     Set the force to use Ewald summation.

     @param alpha  the Ewald separation parameter
     @param kmaxx  the largest wave vector in the x direction
     @param kmaxy  the largest wave vector in the y direction
     @param kmaxz  the largest wave vector in the z direction

     --------------------------------------------------------------------------------------- */

  void ReferenceLJCoulombIxn::setUseEwald(RealOpenMM alpha, int kmaxx, int kmaxy, int kmaxz) {
      alphaEwald = alpha;
      numRx = kmaxx;
      numRy = kmaxy;
      numRz = kmaxz;
      ewald = true;
  }

  /**---------------------------------------------------------------------------------------

     Set the force to use Particle-Mesh Ewald (PME) summation.

     @param alpha  the Ewald separation parameter
     @param gridSize the dimensions of the mesh

     --------------------------------------------------------------------------------------- */

  void ReferenceLJCoulombIxn::setUsePME(RealOpenMM alpha, int meshSize[3]) {
      alphaEwald = alpha;
      meshDim[0] = meshSize[0];
      meshDim[1] = meshSize[1];
      meshDim[2] = meshSize[2];
      pme = true;
  }

/**---------------------------------------------------------------------------------------

   Calculate Ewald ixn

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

void ReferenceLJCoulombIxn::calculateEwaldIxn(int numberOfAtoms, vector<RealVec>& atomCoordinates,
                                             RealOpenMM** atomParameters, vector<set<int> >& exclusions,
                                             RealOpenMM* fixedParameters, vector<RealVec>& forces,
                                             RealOpenMM* energyByAtom, RealOpenMM* totalEnergy, bool includeDirect, bool includeReciprocal) const {
    typedef std::complex<RealOpenMM> d_complex;

    static const RealOpenMM epsilon     =  1.0;
    static const RealOpenMM one         =  1.0;
    static const RealOpenMM six         =  6.0;
    static const RealOpenMM twelve      = 12.0;

    int kmax                            = (ewald ? std::max(numRx, std::max(numRy,numRz)) : 0);
    RealOpenMM  factorEwald             = -1 / (4*alphaEwald*alphaEwald);
    RealOpenMM SQRT_PI                  = sqrt(PI_M);
    RealOpenMM TWO_PI                   = 2.0 * PI_M;
    RealOpenMM recipCoeff               = (RealOpenMM)(ONE_4PI_EPS0*4*PI_M/(periodicBoxSize[0] * periodicBoxSize[1] * periodicBoxSize[2]) /epsilon);

    RealOpenMM totalSelfEwaldEnergy     = 0.0;
    RealOpenMM realSpaceEwaldEnergy     = 0.0;
    RealOpenMM recipEnergy              = 0.0;
    RealOpenMM totalRecipEnergy         = 0.0;
    RealOpenMM vdwEnergy                = 0.0;

// **************************************************************************************
// SELF ENERGY
// **************************************************************************************

    if (includeReciprocal) {
        for( int atomID = 0; atomID < numberOfAtoms; atomID++ ){
            RealOpenMM selfEwaldEnergy       = (RealOpenMM) (ONE_4PI_EPS0*atomParameters[atomID][QIndex]*atomParameters[atomID][QIndex] * alphaEwald/SQRT_PI);
            totalSelfEwaldEnergy            -= selfEwaldEnergy;
            if( energyByAtom ){
                energyByAtom[atomID]        -= selfEwaldEnergy;
            }
        }
    }

    if( totalEnergy ){
        *totalEnergy += totalSelfEwaldEnergy;
    }

// **************************************************************************************
// RECIPROCAL SPACE EWALD ENERGY AND FORCES
// **************************************************************************************
    // PME

  if (pme && includeReciprocal) {
    pme_t          pmedata; /* abstract handle for PME data */
    RealOpenMM virial[3][3];

    pme_init(&pmedata,alphaEwald,numberOfAtoms,meshDim,5,1);

    vector<RealOpenMM> charges(numberOfAtoms);
    for (int i = 0; i < numberOfAtoms; i++)
        charges[i] = atomParameters[i][QIndex];
    pme_exec(pmedata,atomCoordinates,forces,charges,periodicBoxSize,&recipEnergy,virial);

    if( totalEnergy )
       *totalEnergy += recipEnergy;

    if( energyByAtom )
        for(int n = 0; n < numberOfAtoms; n++)
            energyByAtom[n] += recipEnergy;

        pme_destroy(pmedata);
  }

    // Ewald method

  else if (ewald && includeReciprocal) {

    // setup reciprocal box

           RealOpenMM recipBoxSize[3] = { TWO_PI / periodicBoxSize[0], TWO_PI / periodicBoxSize[1], TWO_PI / periodicBoxSize[2]};


    // setup K-vectors

  #define EIR(x, y, z) eir[(x)*numberOfAtoms*3+(y)*3+z]
  vector<d_complex> eir(kmax*numberOfAtoms*3);
  vector<d_complex> tab_xy(numberOfAtoms);
  vector<d_complex> tab_qxyz(numberOfAtoms);

  if (kmax < 1) {
      std::stringstream message;
      message << " kmax < 1 , Aborting" << std::endl;
      SimTKOpenMMLog::printError( message );
  }

  for(int i = 0; (i < numberOfAtoms); i++) {
    for(int m = 0; (m < 3); m++)
      EIR(0, i, m) = d_complex(1,0);

    for(int m=0; (m<3); m++)
      EIR(1, i, m) = d_complex(cos(atomCoordinates[i][m]*recipBoxSize[m]),
                               sin(atomCoordinates[i][m]*recipBoxSize[m]));

    for(int j=2; (j<kmax); j++)
      for(int m=0; (m<3); m++)
        EIR(j, i, m) = EIR(j-1, i, m) * EIR(1, i, m);
  }

    // calculate reciprocal space energy and forces

    int lowry = 0;
    int lowrz = 1;

    for(int rx = 0; rx < numRx; rx++) {

      RealOpenMM kx = rx * recipBoxSize[0];

      for(int ry = lowry; ry < numRy; ry++) {

        RealOpenMM ky = ry * recipBoxSize[1];

        if(ry >= 0) {
          for(int n = 0; n < numberOfAtoms; n++)
            tab_xy[n] = EIR(rx, n, 0) * EIR(ry, n, 1);
        }

        else {
          for(int n = 0; n < numberOfAtoms; n++)
            tab_xy[n]= EIR(rx, n, 0) * conj (EIR(-ry, n, 1));
        }

        for (int rz = lowrz; rz < numRz; rz++) {

          if( rz >= 0) {
           for( int n = 0; n < numberOfAtoms; n++)
             tab_qxyz[n] = atomParameters[n][QIndex] * (tab_xy[n] * EIR(rz, n, 2));
          }

          else {
            for( int n = 0; n < numberOfAtoms; n++)
              tab_qxyz[n] = atomParameters[n][QIndex] * (tab_xy[n] * conj(EIR(-rz, n, 2)));
          }

          RealOpenMM cs = 0.0f;
          RealOpenMM ss = 0.0f;

          for( int n = 0; n < numberOfAtoms; n++) {
            cs += tab_qxyz[n].real();
            ss += tab_qxyz[n].imag();
          }

          RealOpenMM kz = rz * recipBoxSize[2];
          RealOpenMM k2 = kx * kx + ky * ky + kz * kz;
          RealOpenMM ak = exp(k2*factorEwald) / k2;

          for(int n = 0; n < numberOfAtoms; n++) {
            RealOpenMM force = ak * (cs * tab_qxyz[n].imag() - ss * tab_qxyz[n].real());
            forces[n][0] += 2 * recipCoeff * force * kx ;
            forces[n][1] += 2 * recipCoeff * force * ky ;
            forces[n][2] += 2 * recipCoeff * force * kz ;
          }

          recipEnergy       = recipCoeff * ak * ( cs * cs + ss * ss);
          totalRecipEnergy += recipEnergy;

          if( totalEnergy )
             *totalEnergy += recipEnergy;

          if( energyByAtom )
             for(int n = 0; n < numberOfAtoms; n++)
               energyByAtom[n] += recipEnergy;

          lowrz = 1 - numRz;
        }
        lowry = 1 - numRy;
      }
    }
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
       ReferenceForce::getDeltaRPeriodic( atomCoordinates[jj], atomCoordinates[ii], periodicBoxSize, deltaR[0] );
       RealOpenMM r         = deltaR[0][ReferenceForce::RIndex];
       RealOpenMM inverseR  = one/(deltaR[0][ReferenceForce::RIndex]);
       RealOpenMM switchValue = 1, switchDeriv = 0;
       if (useSwitch && r > switchingDistance) {
           RealOpenMM t = (r-switchingDistance)/(cutoffDistance-switchingDistance);
           switchValue = 1+t*t*t*(-10+t*(15-t*6));
           switchDeriv = t*t*(-30+t*(60-t*30))/(cutoffDistance-switchingDistance);
       }
       RealOpenMM alphaR    = alphaEwald * r;


       RealOpenMM dEdR      = (RealOpenMM) (ONE_4PI_EPS0 * atomParameters[ii][QIndex] * atomParameters[jj][QIndex] * inverseR * inverseR * inverseR);
                  dEdR      = (RealOpenMM) (dEdR * (erfc(alphaR) + 2 * alphaR * exp ( - alphaR * alphaR) / SQRT_PI ));

       RealOpenMM sig       = atomParameters[ii][SigIndex] +  atomParameters[jj][SigIndex];
       RealOpenMM sig2      = inverseR*sig;
                  sig2     *= sig2;
       RealOpenMM sig6      = sig2*sig2*sig2;
       RealOpenMM eps       = atomParameters[ii][EpsIndex]*atomParameters[jj][EpsIndex];
                  dEdR     += switchValue*eps*( twelve*sig6 - six )*sig6*inverseR*inverseR;
       vdwEnergy = eps*(sig6-one)*sig6;
       if (useSwitch) {
           dEdR -= vdwEnergy*switchDeriv*inverseR;
           vdwEnergy *= switchValue;
       }

       // accumulate forces

       for( int kk = 0; kk < 3; kk++ ){
          RealOpenMM force  = dEdR*deltaR[0][kk];
          forces[ii][kk]   += force;
          forces[jj][kk]   -= force;
       }

       // accumulate energies

       realSpaceEwaldEnergy        = (RealOpenMM) (ONE_4PI_EPS0*atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*inverseR*erfc(alphaR));

       totalVdwEnergy             += vdwEnergy;
       totalRealSpaceEwaldEnergy  += realSpaceEwaldEnergy;

        if( energyByAtom ){
           energyByAtom[ii] += realSpaceEwaldEnergy + vdwEnergy;
           energyByAtom[jj] += realSpaceEwaldEnergy + vdwEnergy;
        }

    }

    if( totalEnergy )
        *totalEnergy += totalRealSpaceEwaldEnergy + totalVdwEnergy;

    // Now subtract off the exclusions, since they were implicitly included in the reciprocal space sum.

    RealOpenMM totalExclusionEnergy = 0.0f;
    for (int i = 0; i < numberOfAtoms; i++)
        for (set<int>::const_iterator iter = exclusions[i].begin(); iter != exclusions[i].end(); ++iter) {
            if (*iter > i) {
               int ii = i;
               int jj = *iter;

               RealOpenMM deltaR[2][ReferenceForce::LastDeltaRIndex];
               ReferenceForce::getDeltaR( atomCoordinates[jj], atomCoordinates[ii], deltaR[0] );
               RealOpenMM r         = deltaR[0][ReferenceForce::RIndex];
               RealOpenMM inverseR  = one/(deltaR[0][ReferenceForce::RIndex]);
               RealOpenMM alphaR    = alphaEwald * r;
               if (erf(alphaR) > 1e-6) {
                   RealOpenMM dEdR      = (RealOpenMM) (ONE_4PI_EPS0 * atomParameters[ii][QIndex] * atomParameters[jj][QIndex] * inverseR * inverseR * inverseR);
                              dEdR      = (RealOpenMM) (dEdR * (erf(alphaR) - 2 * alphaR * exp ( - alphaR * alphaR) / SQRT_PI ));

                   // accumulate forces

                   for( int kk = 0; kk < 3; kk++ ){
                      RealOpenMM force  = dEdR*deltaR[0][kk];
                      forces[ii][kk]   -= force;
                      forces[jj][kk]   += force;
                   }

                   // accumulate energies

                   realSpaceEwaldEnergy = (RealOpenMM) (ONE_4PI_EPS0*atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*inverseR*erf(alphaR));

                   totalExclusionEnergy += realSpaceEwaldEnergy;
                   if( energyByAtom ){
                       energyByAtom[ii] -= realSpaceEwaldEnergy;
                       energyByAtom[jj] -= realSpaceEwaldEnergy;
                   }
               }
            }
        }

    if( totalEnergy )
        *totalEnergy -= totalExclusionEnergy;
}


/**---------------------------------------------------------------------------------------

   Calculate LJ Coulomb pair ixn

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

void ReferenceLJCoulombIxn::calculatePairIxn(int numberOfAtoms, vector<RealVec>& atomCoordinates,
                                             RealOpenMM** atomParameters, vector<set<int> >& exclusions,
                                             RealOpenMM* fixedParameters, vector<RealVec>& forces,
                                             RealOpenMM* energyByAtom, RealOpenMM* totalEnergy, bool includeDirect, bool includeReciprocal) const {

   if (ewald || pme) {
       calculateEwaldIxn(numberOfAtoms, atomCoordinates, atomParameters, exclusions, fixedParameters, forces, energyByAtom,
               totalEnergy, includeDirect, includeReciprocal);
       return;
   }
   if (!includeDirect)
       return;
   if (cutoff) {
       for (int i = 0; i < (int) neighborList->size(); i++) {
           OpenMM::AtomPair pair = (*neighborList)[i];
           calculateOneIxn(pair.first, pair.second, atomCoordinates, atomParameters, forces, energyByAtom, totalEnergy);
       }
   }
   else {
       for( int ii = 0; ii < numberOfAtoms; ii++ ){
          // loop over atom pairs

          for( int jj = ii+1; jj < numberOfAtoms; jj++ )
              if (exclusions[jj].find(ii) == exclusions[jj].end())
                  calculateOneIxn(ii, jj, atomCoordinates, atomParameters, forces, energyByAtom, totalEnergy);
       }
   }
}

  /**---------------------------------------------------------------------------------------

     Calculate LJ Coulomb pair ixn between two atoms

     @param ii               the index of the first atom
     @param jj               the index of the second atom
     @param atomCoordinates  atom coordinates
     @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
     @param forces           force array (forces added)
     @param energyByAtom     atom energy
     @param totalEnergy      total energy

     --------------------------------------------------------------------------------------- */

void ReferenceLJCoulombIxn::calculateOneIxn( int ii, int jj, vector<RealVec>& atomCoordinates,
                        RealOpenMM** atomParameters, vector<RealVec>& forces,
                        RealOpenMM* energyByAtom, RealOpenMM* totalEnergy ) const {

    // ---------------------------------------------------------------------------------------

    static const std::string methodName = "\nReferenceLJCoulombIxn::calculateOneIxn";

    // ---------------------------------------------------------------------------------------

    // constants -- reduce Visual Studio warnings regarding conversions between float & double

    static const RealOpenMM zero        =  0.0;
    static const RealOpenMM one         =  1.0;
    static const RealOpenMM two         =  2.0;
    static const RealOpenMM three       =  3.0;
    static const RealOpenMM six         =  6.0;
    static const RealOpenMM twelve      = 12.0;
    static const RealOpenMM oneM        = -1.0;

    static const int threeI             = 3;

    static const int LastAtomIndex      = 2;

    RealOpenMM deltaR[2][ReferenceForce::LastDeltaRIndex];

    // get deltaR, R2, and R between 2 atoms

    if (periodic)
        ReferenceForce::getDeltaRPeriodic( atomCoordinates[jj], atomCoordinates[ii], periodicBoxSize, deltaR[0] );
    else
        ReferenceForce::getDeltaR( atomCoordinates[jj], atomCoordinates[ii], deltaR[0] );

    RealOpenMM r2        = deltaR[0][ReferenceForce::R2Index];
    RealOpenMM inverseR  = one/(deltaR[0][ReferenceForce::RIndex]);
    RealOpenMM switchValue = 1, switchDeriv = 0;
    if (useSwitch) {
        RealOpenMM r = deltaR[0][ReferenceForce::RIndex];
        if (r > switchingDistance) {
            RealOpenMM t = (r-switchingDistance)/(cutoffDistance-switchingDistance);
            switchValue = 1+t*t*t*(-10+t*(15-t*6));
            switchDeriv = t*t*(-30+t*(60-t*30))/(cutoffDistance-switchingDistance);
        }
    }
    RealOpenMM sig       = atomParameters[ii][SigIndex] +  atomParameters[jj][SigIndex];
    RealOpenMM sig2      = inverseR*sig;
               sig2     *= sig2;
    RealOpenMM sig6      = sig2*sig2*sig2;

    RealOpenMM eps       = atomParameters[ii][EpsIndex]*atomParameters[jj][EpsIndex];
    RealOpenMM dEdR      = switchValue*eps*( twelve*sig6 - six )*sig6;
    if (cutoff)
        dEdR += (RealOpenMM) (ONE_4PI_EPS0*atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*(inverseR-2.0f*krf*r2));
    else
        dEdR += (RealOpenMM) (ONE_4PI_EPS0*atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*inverseR);
    dEdR     *= inverseR*inverseR;
    RealOpenMM energy = eps*(sig6-one)*sig6;
    if (useSwitch) {
        dEdR -= energy*switchDeriv*inverseR;
        energy *= switchValue;
    }
    if (cutoff)
        energy += (RealOpenMM) (ONE_4PI_EPS0*atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*(inverseR+krf*r2-crf));
    else
        energy += (RealOpenMM) (ONE_4PI_EPS0*atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*inverseR);

    // accumulate forces

    for( int kk = 0; kk < 3; kk++ ){
       RealOpenMM force  = dEdR*deltaR[0][kk];
       forces[ii][kk]   += force;
       forces[jj][kk]   -= force;
    }

    // accumulate energies

    if( totalEnergy )
       *totalEnergy += energy;
    if( energyByAtom ){
       energyByAtom[ii] += energy;
       energyByAtom[jj] += energy;
    }
  }



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

#include <string.h>
#include <sstream>
#include <complex>

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "ReferenceFreeEnergyLJCoulombSoftcoreIxn.h"
#include "ReferenceForce.h"
#include "PME.h"

// In case we're using some primitive version of Visual Studio this will
// make sure that erf() and erfc() are defined.
//#include "MSVC_erfc.h"
#include "openmm/internal/MSVC_erfc.h"

using std::vector;

/**---------------------------------------------------------------------------------------

   ReferenceFreeEnergyLJCoulombSoftcoreIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceFreeEnergyLJCoulombSoftcoreIxn::ReferenceFreeEnergyLJCoulombSoftcoreIxn( ) : cutoff(false), periodic(false), ewald(false), pme(false), softCoreLJLambda(1.0) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceFreeEnergyLJCoulombSoftcoreIxn::ReferenceFreeEnergyLJCoulombSoftcoreIxn";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   ReferenceFreeEnergyLJCoulombSoftcoreIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceFreeEnergyLJCoulombSoftcoreIxn::~ReferenceFreeEnergyLJCoulombSoftcoreIxn( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceFreeEnergyLJCoulombSoftcoreIxn::~ReferenceFreeEnergyLJCoulombSoftcoreIxn";

   // ---------------------------------------------------------------------------------------

}

  /**---------------------------------------------------------------------------------------

     Set the force to use a cutoff.

     @param distance            the cutoff distance
     @param neighbors           the neighbor list to use
     @param solventDielectric   the dielectric constant of the bulk solvent

     @return ReferenceForce::DefaultReturn

     --------------------------------------------------------------------------------------- */

  int ReferenceFreeEnergyLJCoulombSoftcoreIxn::setUseCutoff( RealOpenMM distance, const OpenMM::NeighborList& neighbors, RealOpenMM solventDielectric ) {

    cutoff = true;
    cutoffDistance = distance;
    neighborList = &neighbors;
    krf = pow(cutoffDistance, (RealOpenMM)-3.0)*(solventDielectric-1.0f)/(2.0f*solventDielectric+1.0f);
    crf = (1.0f/cutoffDistance)*(3.0f*solventDielectric)/(2.0f*solventDielectric+1.0f);

    return ReferenceForce::DefaultReturn;
  }

  /**---------------------------------------------------------------------------------------

     Set the force to use periodic boundary conditions.  This requires that a cutoff has
     also been set, and the smallest side of the periodic box is at least twice the cutoff
     distance.

     @param boxSize             the X, Y, and Z widths of the periodic box

     @return ReferenceForce::DefaultReturn

     --------------------------------------------------------------------------------------- */

  int ReferenceFreeEnergyLJCoulombSoftcoreIxn::setPeriodic( RealOpenMM* boxSize ) {

    assert(cutoff);
    assert(boxSize[0] >= 2.0*cutoffDistance);
    assert(boxSize[1] >= 2.0*cutoffDistance);
    assert(boxSize[2] >= 2.0*cutoffDistance);
    periodic = true;
    periodicBoxSize[0] = boxSize[0];
    periodicBoxSize[1] = boxSize[1];
    periodicBoxSize[2] = boxSize[2];
    return ReferenceForce::DefaultReturn;

  }

  /**---------------------------------------------------------------------------------------

     Set the force to use Ewald summation.

     @param alpha  the Ewald separation parameter
     @param kmaxx  the largest wave vector in the x direction
     @param kmaxy  the largest wave vector in the y direction
     @param kmaxz  the largest wave vector in the z direction

     --------------------------------------------------------------------------------------- */

  void ReferenceFreeEnergyLJCoulombSoftcoreIxn::setUseEwald(RealOpenMM alpha, int kmaxx, int kmaxy, int kmaxz) {
      alphaEwald = alpha;
      numRx = kmaxx;
      numRy = kmaxy;
      numRz = kmaxz;
      ewald = true;
  }

  /**---------------------------------------------------------------------------------------

     Set the force to use Particle-Mesh Ewald (PME) summation.

     @param alpha  the Ewald separation parameter

     --------------------------------------------------------------------------------------- */

  void ReferenceFreeEnergyLJCoulombSoftcoreIxn::setUsePME(RealOpenMM alpha) {
      alphaEwald = alpha;
      pme = true;
  }


  /**---------------------------------------------------------------------------------------

     Set the soft core LJ lambda

     @param lambda the soft core LJ lambda

     --------------------------------------------------------------------------------------- */

  void ReferenceFreeEnergyLJCoulombSoftcoreIxn::setSoftCoreLJLambda(RealOpenMM lambda) {
      softCoreLJLambda = lambda;
  }

/**---------------------------------------------------------------------------------------

   Calculate parameters for LJ Coulomb ixn

   @param c6               c6
   @param c12              c12
   @param q1               q1 charge atom 1
   @param epsfac           epsfacSqrt ????????????
   @param parameters       output parameters:
										parameter[SigIndex]  = 0.5*( (c12/c6)**1/6 ) (sigma/2)
										parameter[EpsIndex]  = sqrt(c6*c6/c12)       (2*sqrt(epsilon))
										parameter[QIndex]    = epsfactorSqrt*q1

   @return ReferenceForce::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceFreeEnergyLJCoulombSoftcoreIxn::getDerivedParameters( RealOpenMM c6, RealOpenMM c12, RealOpenMM q1,
                                                                   RealOpenMM epsfacSqrt,
                                                                   RealOpenMM* parameters ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceFreeEnergyLJCoulombSoftcoreIxn::getDerivedParameters";

   static const RealOpenMM zero          =  0.0;
   static const RealOpenMM one           =  1.0;
   static const RealOpenMM six           =  6.0;
   static const RealOpenMM half          =  0.5;
   static const RealOpenMM oneSixth      =  one/six;
   static const RealOpenMM oneTweleth    =  half*oneSixth;

   // ---------------------------------------------------------------------------------------

   if( c12 <= 0.0 ){

      parameters[EpsIndex] = zero;
      parameters[SigIndex] = half;

   } else {

      parameters[EpsIndex]    = c6*SQRT( one/c12 );

      parameters[SigIndex]    = POW( (c12/c6), oneSixth );
      parameters[SigIndex]   *= half;
   }

   parameters[QIndex]   = epsfacSqrt*q1;

   return ReferenceForce::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Calculate Ewald ixn

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param atomParameters   atom parameters                             atomParameters[atomIndex][paramterIndex]
   @param exclusions       atom exclusion indices                      exclusions[atomIndex][atomToExcludeIndex]
                           exclusions[atomIndex][0] = number of exclusions
                           exclusions[atomIndex][1-no.] = atom indices of atoms to excluded from
                           interacting w/ atom atomIndex
   @param fixedParameters  non atom parameters (not currently used)
   @param forces           force array (forces added)
   @param energyByAtom     atom energy
   @param totalEnergy      total energy

   @return ReferenceForce::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceFreeEnergyLJCoulombSoftcoreIxn::calculateEwaldIxn( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                                                 RealOpenMM** atomParameters, int** exclusions,
                                                                 RealOpenMM* fixedParameters, RealOpenMM** forces,
                                                                 RealOpenMM* energyByAtom, RealOpenMM* totalEnergy) const {
#if 0
    typedef std::complex<RealOpenMM> d_complex;

    static const RealOpenMM epsilon     =  1.0;
    static const RealOpenMM one         =  1.0;
    static const RealOpenMM six         =  6.0;
    static const RealOpenMM twelve      = 12.0;

    int kmax                            = (ewald ? std::max(numRx, std::max(numRy,numRz)) : 0);
    RealOpenMM  factorEwald             = -1 / (4*alphaEwald*alphaEwald);
    RealOpenMM SQRT_PI                  = sqrt(PI);
    RealOpenMM TWO_PI                   = 2.0 * PI;
    RealOpenMM recipCoeff               = (RealOpenMM)(4*PI/(periodicBoxSize[0] * periodicBoxSize[1] * periodicBoxSize[2]) /epsilon);

    RealOpenMM totalSelfEwaldEnergy     = 0.0;
    RealOpenMM realSpaceEwaldEnergy     = 0.0;
    RealOpenMM recipEnergy              = 0.0;
    RealOpenMM totalRecipEnergy         = 0.0;
    RealOpenMM vdwEnergy                = 0.0;

// **************************************************************************************
// SELF ENERGY
// **************************************************************************************

    for( int atomID = 0; atomID < numberOfAtoms; atomID++ ){
        RealOpenMM selfEwaldEnergy       = atomParameters[atomID][QIndex]*atomParameters[atomID][QIndex] * alphaEwald/SQRT_PI;
        totalSelfEwaldEnergy            -= selfEwaldEnergy;

        if( energyByAtom ){
           energyByAtom[atomID]         -= selfEwaldEnergy;
        }
    }

    if( totalEnergy ){
        *totalEnergy += totalSelfEwaldEnergy;
    }

// **************************************************************************************
// RECIPROCAL SPACE EWALD ENERGY AND FORCES
// **************************************************************************************
    // PME

  if (pme) {
	pme_t          pmedata; /* abstract handle for PME data */
	int            ngrid[3];
	RealOpenMM virial[3][3];

	/* PME grid dimensions.
	 * We typically want to set this as the spacing rather than absolute dimensions, but
	 * to be able to reproduce results from other programs (e.g. Gromacs) we need to be
	 * able to set exact grid dimenisions occasionally.
	 */
	ngrid[0] = 16;
	ngrid[1] = 16;
	ngrid[2] = 16;

	pme_init(&pmedata,alphaEwald,numberOfAtoms,ngrid,4,1);

	pme_exec(pmedata,atomCoordinates,forces,atomParameters,periodicBoxSize,&recipEnergy,virial);

	if( totalEnergy )
       *totalEnergy += recipEnergy;

    if( energyByAtom )
        for(int n = 0; n < numberOfAtoms; n++)
            energyByAtom[n] += recipEnergy;

        pme_destroy(pmedata);
  }

    // Ewald method

  else if (ewald) {

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

  else {
      std::stringstream message;
      message << " Wrong method for Ewald summation, Aborting" << std::endl;
      SimTKOpenMMLog::printError( message );
  }


// **************************************************************************************
// SHORT-RANGE ENERGY AND FORCES
// **************************************************************************************

    RealOpenMM totalVdwEnergy            = 0.0f;
    RealOpenMM totalRealSpaceEwaldEnergy = 0.0f;

    for (int i = 0; i < (int) neighborList->size(); i++) {
       OpenMM::AtomPair pair = (*neighborList)[i];
       int ii = pair.first;
       int jj = pair.second;

       RealOpenMM deltaR[2][ReferenceForce::LastDeltaRIndex];
       ReferenceForce::getDeltaRPeriodic( atomCoordinates[jj], atomCoordinates[ii], periodicBoxSize, deltaR[0] );
       RealOpenMM r         = deltaR[0][ReferenceForce::RIndex];
       RealOpenMM r2        = deltaR[0][ReferenceForce::R2Index];
       RealOpenMM inverseR  = one/(deltaR[0][ReferenceForce::RIndex]);
       RealOpenMM alphaR    = alphaEwald * r;


       RealOpenMM dEdR      = atomParameters[ii][QIndex] * atomParameters[jj][QIndex] * inverseR * inverseR * inverseR;
                  dEdR      = (RealOpenMM)(dEdR * (erfc(alphaR) + 2 * alphaR * exp ( - alphaR * alphaR) / SQRT_PI ));

       RealOpenMM sig       = atomParameters[ii][SigIndex] +  atomParameters[jj][SigIndex];
       RealOpenMM sig2      = inverseR*sig;
                  sig2     *= sig2;
       RealOpenMM sig6      = sig2*sig2*sig2;
       RealOpenMM eps       = atomParameters[ii][EpsIndex]*atomParameters[jj][EpsIndex];
                  dEdR     += eps*( twelve*sig6 - six )*sig6*inverseR*inverseR;

       // accumulate forces

       for( int kk = 0; kk < 3; kk++ ){
          RealOpenMM force  = dEdR*deltaR[0][kk];
          forces[ii][kk]   += force;
          forces[jj][kk]   -= force;
       }

       // accumulate energies

       realSpaceEwaldEnergy        = (RealOpenMM) (atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*inverseR*erfc(alphaR));
       vdwEnergy                   = eps*(sig6-one)*sig6;

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
        for (int j = 1; j <= exclusions[i][0]; j++)
            if (exclusions[i][j] > i) {
               int ii = i;
               int jj = exclusions[i][j];

               RealOpenMM deltaR[2][ReferenceForce::LastDeltaRIndex];
               ReferenceForce::getDeltaRPeriodic( atomCoordinates[jj], atomCoordinates[ii], periodicBoxSize, deltaR[0] );
               RealOpenMM r         = deltaR[0][ReferenceForce::RIndex];
               RealOpenMM inverseR  = one/(deltaR[0][ReferenceForce::RIndex]);
               RealOpenMM alphaR    = alphaEwald * r;
               RealOpenMM dEdR      = atomParameters[ii][QIndex] * atomParameters[jj][QIndex] * inverseR * inverseR * inverseR;
                          dEdR      = (RealOpenMM)(dEdR * (erf(alphaR) - 2 * alphaR * exp ( - alphaR * alphaR) / SQRT_PI ));

               // accumulate forces

               for( int kk = 0; kk < 3; kk++ ){
                  RealOpenMM force  = dEdR*deltaR[0][kk];
                  forces[ii][kk]   -= force;
                  forces[jj][kk]   += force;
               }

               // accumulate energies

               realSpaceEwaldEnergy = (RealOpenMM) (atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*inverseR*erf(alphaR));

               totalExclusionEnergy += realSpaceEwaldEnergy;
               if( energyByAtom ){
                   energyByAtom[ii] -= realSpaceEwaldEnergy;
                   energyByAtom[jj] -= realSpaceEwaldEnergy;
               }
            }

    if( totalEnergy )
        *totalEnergy -= totalExclusionEnergy;

   // ***********************************************************************

#endif
   return ReferenceForce::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Calculate PME ixn

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param atomParameters   atom parameters                             atomParameters[atomIndex][paramterIndex]
   @param exclusions       atom exclusion indices                      exclusions[atomIndex][atomToExcludeIndex]
                           exclusions[atomIndex][0] = number of exclusions
                           exclusions[atomIndex][1-no.] = atom indices of atoms to excluded from
                           interacting w/ atom atomIndex
   @param fixedParameters  non atom parameters (not currently used)
   @param forces           force array (forces added)
   @param energyByAtom     atom energy
   @param totalEnergy      total energy

   @return ReferenceForce::DefaultReturn
      
   --------------------------------------------------------------------------------------- */
 
int ReferenceFreeEnergyLJCoulombSoftcoreIxn::calculatePMEIxn( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                                              RealOpenMM** atomParameters, int** exclusions,
                                                              RealOpenMM* fixedParameters, RealOpenMM** forces,
                                                              RealOpenMM* energyByAtom, RealOpenMM* totalEnergy) const {


#if 0
    RealOpenMM SQRT_PI = sqrt(PI);
    static const RealOpenMM one         =  1.0;

    RealOpenMM selfEwaldEnergy = 0.0;
    RealOpenMM realSpaceEwaldEnergy = 0.0;
    RealOpenMM recipEnergy = 0.0;


// **************************************************************************************
// SELF ENERGY
// **************************************************************************************

    for( int atomID = 0; atomID < numberOfAtoms; atomID++ ){
        selfEwaldEnergy = selfEwaldEnergy + atomParameters[atomID][QIndex]*atomParameters[atomID][QIndex];
    }
       selfEwaldEnergy = selfEwaldEnergy * alphaEwald/SQRT_PI ;

// **************************************************************************************
// RECIPROCAL SPACE EWALD ENERGY AND FORCES
// **************************************************************************************

	pme_t          pmedata; /* abstract handle for PME data */
	int            ngrid[3];
	RealOpenMM virial[3][3];
	
	/* PME grid dimensions.
	 * We typically want to set this as the spacing rather than absolute dimensions, but
	 * to be able to reproduce results from other programs (e.g. Gromacs) we need to be
	 * able to set exact grid dimenisions occasionally.
	 */
	ngrid[0] = 16;
	ngrid[1] = 16;
	ngrid[2] = 16;
	
	pme_init(&pmedata,alphaEwald,numberOfAtoms,ngrid,4,1);
	
	pme_exec(pmedata,atomCoordinates,forces,atomParameters,periodicBoxSize,&recipEnergy,virial);
	
// **************************************************************************************
// SHORT-RANGE ENERGY AND FORCES
// **************************************************************************************

       RealOpenMM deltaR[2][ReferenceForce::LastDeltaRIndex];

       for( int atomID1 = 0; atomID1 < numberOfAtoms; atomID1++ ){
        for( int atomID2 = atomID1 + 1; atomID2 < numberOfAtoms; atomID2++ ){

          ReferenceForce::getDeltaRPeriodic( atomCoordinates[atomID2], atomCoordinates[atomID1], periodicBoxSize, deltaR[0] );  
          RealOpenMM r         = deltaR[0][ReferenceForce::RIndex];
          RealOpenMM r2        = deltaR[0][ReferenceForce::R2Index];
          RealOpenMM inverseR  = one/(deltaR[0][ReferenceForce::RIndex]);

          realSpaceEwaldEnergy = 
              (RealOpenMM)(realSpaceEwaldEnergy + atomParameters[atomID1][QIndex]*atomParameters[atomID2][QIndex]*inverseR*erfc(alphaEwald*r)); 
        }
       }

       // allocate and initialize exclusion array

       vector<int> exclusionIndices(numberOfAtoms);
       for( int ii = 0; ii < numberOfAtoms; ii++ ){
          exclusionIndices[ii] = -1;
       }

       for( int ii = 0; ii < numberOfAtoms; ii++ ){

          // set exclusions

          for( int jj = 1; jj <= exclusions[ii][0]; jj++ ){
             exclusionIndices[exclusions[ii][jj]] = ii;
          }

          // loop over atom pairs

          for( int jj = ii+1; jj < numberOfAtoms; jj++ ){

             if( exclusionIndices[jj] != ii ){

       ReferenceForce::getDeltaRPeriodic( atomCoordinates[jj], atomCoordinates[ii], periodicBoxSize, deltaR[0] );  
       RealOpenMM r         = deltaR[0][ReferenceForce::RIndex];
       RealOpenMM r2        = deltaR[0][ReferenceForce::R2Index];
       RealOpenMM inverseR  = one/(deltaR[0][ReferenceForce::RIndex]);
       RealOpenMM alphaR    = alphaEwald * r;
 
       realSpaceEwaldEnergy = 
           (RealOpenMM)(realSpaceEwaldEnergy + atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*inverseR*erfc(alphaR)); 
       RealOpenMM dEdR = atomParameters[ii][QIndex] * atomParameters[jj][QIndex] * inverseR * inverseR * inverseR;
       dEdR = (RealOpenMM)(dEdR * (erfc(alphaR) + 2 * alphaR * exp ( - alphaR * alphaR) / SQRT_PI ));
 
                for( int kk = 0; kk < 3; kk++ ){
                   RealOpenMM force  = dEdR*deltaR[0][kk];
                   forces[ii][kk]   += force;
                   forces[jj][kk]   -= force;
                } 
             }
          }
       }

// ***********************************************************************

#endif
   return ReferenceForce::DefaultReturn;
}


/**---------------------------------------------------------------------------------------

   Calculate LJ Coulomb pair ixn

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param atomParameters   atom parameters                             atomParameters[atomIndex][paramterIndex]
   @param exclusions       atom exclusion indices                      exclusions[atomIndex][atomToExcludeIndex]
                           exclusions[atomIndex][0] = number of exclusions
                           exclusions[atomIndex][1-no.] = atom indices of atoms to excluded from
                           interacting w/ atom atomIndex
   @param fixedParameters  non atom parameters (not currently used)
   @param forces           force array (forces added)
   @param energyByAtom     atom energy
   @param totalEnergy      total energy

   @return ReferenceForce::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceFreeEnergyLJCoulombSoftcoreIxn::calculatePairIxn( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                             RealOpenMM** atomParameters, int** exclusions,
                                             RealOpenMM* fixedParameters, RealOpenMM** forces,
                                             RealOpenMM* energyByAtom, RealOpenMM* totalEnergy ) const {

   if (ewald || pme)
        return calculateEwaldIxn(numberOfAtoms, atomCoordinates, atomParameters, exclusions, fixedParameters, forces, energyByAtom, totalEnergy);
   if (cutoff) {
       for (int i = 0; i < (int) neighborList->size(); i++) {
           OpenMM::AtomPair pair = (*neighborList)[i];
           calculateOneIxn(pair.first, pair.second, atomCoordinates, atomParameters, forces, energyByAtom, totalEnergy);
       }
   }
   else {
       // allocate and initialize exclusion array

       int* exclusionIndices = new int[numberOfAtoms];
       for( int ii = 0; ii < numberOfAtoms; ii++ ){
          exclusionIndices[ii] = -1;
       }

       for( int ii = 0; ii < numberOfAtoms; ii++ ){

          // set exclusions

          for( int jj = 1; jj <= exclusions[ii][0]; jj++ ){
             exclusionIndices[exclusions[ii][jj]] = ii;
          }

          // loop over atom pairs

          for( int jj = ii+1; jj < numberOfAtoms; jj++ ){

             if( exclusionIndices[jj] != ii ){
                 calculateOneIxn(ii, jj, atomCoordinates, atomParameters, forces, energyByAtom, totalEnergy);
             }
          }
       }

       delete[] exclusionIndices;
   }

   return ReferenceForce::DefaultReturn;
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

     @return ReferenceForce::DefaultReturn

     --------------------------------------------------------------------------------------- */

int ReferenceFreeEnergyLJCoulombSoftcoreIxn::calculateOneIxn( int ii, int jj, RealOpenMM** atomCoordinates,
                        RealOpenMM** atomParameters, RealOpenMM** forces,
                        RealOpenMM* energyByAtom, RealOpenMM* totalEnergy ) const {

    // ---------------------------------------------------------------------------------------

    static const std::string methodName = "\nReferenceFreeEnergyLJCoulombSoftcoreIxn::calculateOneIxn";

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

    // debug flag

    static const int debug              = -1;

    static const int LastAtomIndex      = 2;

    RealOpenMM deltaR[2][ReferenceForce::LastDeltaRIndex];

    // get deltaR, R2, and R between 2 atoms

    if (periodic)
        ReferenceForce::getDeltaRPeriodic( atomCoordinates[jj], atomCoordinates[ii], periodicBoxSize, deltaR[0] );
    else
        ReferenceForce::getDeltaR( atomCoordinates[jj], atomCoordinates[ii], deltaR[0] );

    RealOpenMM r2                     = deltaR[0][ReferenceForce::R2Index];
    RealOpenMM inverseR               = one/(deltaR[0][ReferenceForce::RIndex]);
    RealOpenMM sig                    = atomParameters[ii][SigIndex] +  atomParameters[jj][SigIndex];
    RealOpenMM eps                    = atomParameters[ii][EpsIndex]*atomParameters[jj][EpsIndex];
    RealOpenMM minSoftCoreLJLambda    = atomParameters[ii][SoftCoreLJLambdaIndex] < atomParameters[jj][SoftCoreLJLambdaIndex] ?
                                        atomParameters[ii][SoftCoreLJLambdaIndex] : atomParameters[jj][SoftCoreLJLambdaIndex];

    // LJ: use soft core LJ if lambda < 1

    RealOpenMM energy          = zero;
    RealOpenMM dEdR;

    if( minSoftCoreLJLambda < one ){
       calculateOneSoftCoreLJIxn( deltaR[0][ReferenceForce::RIndex], sig, eps, minSoftCoreLJLambda, &dEdR, &energy );
    } else {
       calculateOneLJIxn( inverseR, sig, eps, &dEdR, &energy );
    }

    // Coulomb

    if (cutoff)
       dEdR       += atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*(inverseR-2.0f*krf*r2);
    else
       dEdR       += atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*inverseR;

    dEdR          *= inverseR*inverseR;

    // accumulate forces

    for( int kk = 0; kk < 3; kk++ ){
       RealOpenMM force  = dEdR*deltaR[0][kk];
       forces[ii][kk]   += force;
       forces[jj][kk]   -= force;
    }

    // accumulate energies

    if( totalEnergy || energyByAtom ) {
        if (cutoff)
            energy += atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*(inverseR+krf*r2-crf);
        else
            energy += atomParameters[ii][QIndex]*atomParameters[jj][QIndex]*inverseR;

        if( totalEnergy )
           *totalEnergy += energy;
        if( energyByAtom ){
           energyByAtom[ii] += energy;
           energyByAtom[jj] += energy;
        }
    }

    // debug

    if( debug == ii ){
       static bool printHeader = false;
       std::stringstream message;
       message << methodName;
       message << std::endl;
       int pairArray[2] = { ii, jj };
       if( !printHeader  ){
          printHeader = true;
          message << std::endl;
          message << methodName.c_str() << " a0 k [c q p s] r1 r2  angle dt rp p[] dot cosine angle dEdR*r F[]" << std::endl;
       }

       message << std::endl;
       for( int kk = 0; kk < 2; kk++ ){
          message << " Atm " << pairArray[kk] << " [" << atomCoordinates[pairArray[kk]][0] << " " << atomCoordinates[pairArray[kk]][1] << " " << atomCoordinates[pairArray[kk]][2] << "] ";
       }
       message << std::endl << " Delta:";
       for( int kk = 0; kk < (LastAtomIndex - 1); kk++ ){
          message << " [";
          for( int jj = 0; jj < ReferenceForce::LastDeltaRIndex; jj++ ){
             message << deltaR[kk][jj] << " ";
          }
          message << "]";
       }
       message << std::endl;

       for( int kk = 0; kk < 2; kk++ ){
          message << " p" << pairArray[kk] << " [";
          message << atomParameters[pairArray[kk]][0] << " " << atomParameters[pairArray[kk]][1] << " " << atomParameters[pairArray[kk]][2];
          message << "]";
       }
      message << std::endl;

       message << " dEdR=" << dEdR;
       message << " E=" << energy << " force factors: ";
       message << "F=compute force; f=cumulative force";

       message << std::endl << "  ";
       message << " f" << ii << "[";
       SimTKOpenMMUtilities::formatRealStringStream( message, deltaR[0], threeI, dEdR );
       message << "]";

       for( int kk = 0; kk < 2; kk++ ){
          message << " F" <<  pairArray[kk] << " [";
          SimTKOpenMMUtilities::formatRealStringStream( message, forces[pairArray[kk]], threeI );
          message << "]";
       }

       //SimTKOpenMMLog::printMessage( message );
    }
    return ReferenceForce::DefaultReturn;
  }

  /**---------------------------------------------------------------------------------------

     Calculate LJ pair ixn between two atoms

     @param inverseR         1/r
     @param sig              sigma
     @param eps              epsilon
     @param dEdR             output force factor
     @param energy           LJ energy

     @return ReferenceForce::DefaultReturn

     --------------------------------------------------------------------------------------- */

int ReferenceFreeEnergyLJCoulombSoftcoreIxn::calculateOneLJIxn( RealOpenMM inverseR, RealOpenMM sig, RealOpenMM eps,
                        RealOpenMM* dEdR, RealOpenMM* energy ) const {

    // ---------------------------------------------------------------------------------------

    //static const std::string methodName = "\nReferenceFreeEnergyLJCoulombSoftcoreIxn::calculateOneLJIxn";

    // ---------------------------------------------------------------------------------------

    // constants -- reduce Visual Studio warnings regarding conversions between float & double

    static const RealOpenMM zero        =  0.0;
    static const RealOpenMM one         =  1.0;
    static const RealOpenMM six         =  6.0;
    static const RealOpenMM twelve      = 12.0;

    RealOpenMM sig2                     = inverseR*sig;
               sig2                    *= sig2;
    RealOpenMM sig6                     = sig2*sig2*sig2;
              *dEdR                     = eps*( twelve*sig6 - six )*sig6;
               
        *energy                        += eps*(sig6-one)*sig6;

    return ReferenceForce::DefaultReturn;
}

  /**---------------------------------------------------------------------------------------

     Calculate softcore LJ pair ixn between two atoms

     @param r                r
     @param sig              sigma
     @param eps              epsilon
     @param lambda           lambda
     @param dEdR             output force factor
     @param energy           LJ energy

     @return ReferenceForce::DefaultReturn

     --------------------------------------------------------------------------------------- */

int ReferenceFreeEnergyLJCoulombSoftcoreIxn::calculateOneSoftCoreLJIxn( RealOpenMM r, RealOpenMM sig, RealOpenMM eps,
                                                                        RealOpenMM lambda,
                                                                        RealOpenMM* dEdR, RealOpenMM* energy ) const {

    // ---------------------------------------------------------------------------------------

    //static const std::string methodName = "\nReferenceFreeEnergyLJCoulombSoftcoreIxn::calculateOneSoftCoreLJIxn";

    // ---------------------------------------------------------------------------------------

    // constants -- reduce Visual Studio warnings regarding conversions between float & double

    static const RealOpenMM zero        =  0.0;
    static const RealOpenMM one         =  1.0;
    static const RealOpenMM four        =  4.0;
    static const RealOpenMM six         =  6.0;
    static const RealOpenMM twelve      = 12.0;
    static const RealOpenMM alphaLJ     = 0.5;

#if 0
RealOpenMM dEdROrig = 0.0;
RealOpenMM E_Orig   = 0.0;
static int maxPrint = 0;
calculateOneLJIxn( one/r, sig, eps, &dEdROrig, &E_Orig );
#endif

    // soft-core LJ energy = lambda*4*eps*[ 1/{alphaLJ*(1-lambda) + (r/sig)**6}**2 - 1/{alphaLJ*(1-lambda) + (r/sig)**6} ]

    eps                                *= lambda;
    RealOpenMM sig2                     = r/sig;
               sig2                    *= sig2;
    RealOpenMM sig6                     = sig2*sig2*sig2;

    RealOpenMM softcoreLJTerm           = alphaLJ*(one -  lambda) + sig6;
    RealOpenMM softcoreLJInv            = one/softcoreLJTerm;
    RealOpenMM softcoreLJInv2           = softcoreLJInv*softcoreLJInv;

    *dEdR                               = eps*softcoreLJInv2*( twelve*softcoreLJInv - six )*sig6;
               
    *energy                            += eps*softcoreLJInv*( softcoreLJInv - one );

#if 0
if( maxPrint++ < 5 ){
   printf( "r=%14.6e sig=%14.6e eps=%14.6e lambda=%14.6e de[%14.6e %14.6e] e[%14.6e %14.6e] %14.6e %14.6e\n",
           r, sig, eps/lambda, lambda, dEdROrig, *dEdR, E_Orig, *energy, softcoreLJInv, sig6 ); 
}
#endif
    return ReferenceForce::DefaultReturn;
}

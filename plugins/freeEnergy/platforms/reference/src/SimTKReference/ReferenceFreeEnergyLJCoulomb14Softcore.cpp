
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

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "ReferenceFreeEnergyLJCoulomb14Softcore.h"
#include "ReferenceForce.h"

using std::vector;
using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

   ReferenceFreeEnergyLJCoulomb14Softcore constructor

   --------------------------------------------------------------------------------------- */

ReferenceFreeEnergyLJCoulomb14Softcore::ReferenceFreeEnergyLJCoulomb14Softcore( ) : cutoff(false) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceFreeEnergyLJCoulomb14Softcore::ReferenceFreeEnergyLJCoulomb14Softcore";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   ReferenceFreeEnergyLJCoulomb14Softcore destructor

   --------------------------------------------------------------------------------------- */

ReferenceFreeEnergyLJCoulomb14Softcore::~ReferenceFreeEnergyLJCoulomb14Softcore( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceFreeEnergyLJCoulomb14Softcore::~ReferenceFreeEnergyLJCoulomb14Softcore";

   // ---------------------------------------------------------------------------------------

}

  /**---------------------------------------------------------------------------------------

     Set the force to use a cutoff.

     @param distance            the cutoff distance
     @param solventDielectric   the dielectric constant of the bulk solvent

     --------------------------------------------------------------------------------------- */

  void ReferenceFreeEnergyLJCoulomb14Softcore::setUseCutoff( RealOpenMM distance, RealOpenMM solventDielectric ) {
    
    cutoff = true;
    cutoffDistance = distance;
    krf = pow(cutoffDistance, (RealOpenMM)-3.0)*(solventDielectric-1.0f)/(2.0f*solventDielectric+1.0f);
    crf = (1.0f/cutoffDistance)*(3.0f*solventDielectric)/(2.0f*solventDielectric+1.0f);
  }
  
/**---------------------------------------------------------------------------------------

   Calculate parameters for LJ 1-4 ixn

   @param c6               c6
   @param c12              c12
   @param q1               q1 charge atom 1
   @param q2               q2 charge atom 2
   @param epsfac           epsfac ????????????
   @param parameters       output parameters:
										parameter[0]= c6*c6/c12
										parameter[1]= (c12/c6)**1/6
										parameter[2]= epsfactor*q1*q2

   --------------------------------------------------------------------------------------- */

void ReferenceFreeEnergyLJCoulomb14Softcore::getDerivedParameters( RealOpenMM c6, RealOpenMM c12, RealOpenMM q1,
                                                                   RealOpenMM q2, RealOpenMM epsfac,
                                                                   RealOpenMM* parameters ) const {

   // ---------------------------------------------------------------------------------------

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;
   static const RealOpenMM six        =  6.0;
   static const RealOpenMM oneSixth   =  one/six;

   // ---------------------------------------------------------------------------------------

   if( c12 <= zero ){
      parameters[0] = one;
      parameters[1] = zero;
   } else {
      parameters[0] = (c6*c6)/c12;
      parameters[1] = POW( (c12/c6), oneSixth );
   }
   parameters[2] = epsfac*q1*q2;
}

/**---------------------------------------------------------------------------------------

   Calculate LJ 1-4 ixn

   @param atomIndices      atom indices of 4 atoms in bond
   @param atomCoordinates  atom coordinates
   @param parameters       three parameters:
                                        parameters[0]= (c12/c6)**1/6  (sigma)
										parameters[1]= c6*c6/c12      (4*epsilon)
										parameters[2]= epsfac*q1*q2
   @param forces           force array (forces added to current values)
   @param totalEnergy      if not null, the energy will be added to this

   --------------------------------------------------------------------------------------- */

void ReferenceFreeEnergyLJCoulomb14Softcore::calculateBondIxn( int* atomIndices, vector<RealVec>& atomCoordinates,
                                                               RealOpenMM* parameters, vector<RealVec>& forces,
                                                               RealOpenMM* totalEnergy ) const {

   static const std::string methodName = "\nReferenceFreeEnergyLJCoulomb14Softcore::calculateBondIxn";

   // constants -- reduce Visual Studio warnings regarding conversions between float & double

   static const RealOpenMM zero        =  0.0;
   static const RealOpenMM one         =  1.0;
   static const RealOpenMM two         =  2.0;
   static const RealOpenMM three       =  3.0;
   static const RealOpenMM six         =  6.0;
   static const RealOpenMM twelve      = 12.0;
   static const RealOpenMM oneM        = -1.0;

   static const int threeI             = 3;

   // number of parameters

   static const int numberOfParameters = 3;

   // debug flag

   static const int debug              = 0;

   static const int LastAtomIndex      = 2;

   RealOpenMM deltaR[2][ReferenceForce::LastDeltaRIndex];

   // ---------------------------------------------------------------------------------------

   // get deltaR, R2, and R between 2 atoms

   int atomAIndex = atomIndices[0];
   int atomBIndex = atomIndices[1];
   ReferenceForce::getDeltaR( atomCoordinates[atomBIndex], atomCoordinates[atomAIndex], deltaR[0] );  

   if (cutoff && deltaR[0][ReferenceForce::RIndex] > cutoffDistance)
       return;
   RealOpenMM r2        = deltaR[0][ReferenceForce::R2Index];
   RealOpenMM inverseR  = one/(deltaR[0][ReferenceForce::RIndex]);

   RealOpenMM sig                    = parameters[0];
   RealOpenMM eps                    = parameters[1];
   RealOpenMM minSoftCoreLJLambda    = parameters[3];
   RealOpenMM energy                 = zero;
   RealOpenMM dEdR                   = zero;

   if( minSoftCoreLJLambda < one ){
       calculateOneSoftCoreLJ14Ixn( deltaR[0][ReferenceForce::RIndex], sig, eps, minSoftCoreLJLambda, &dEdR, &energy );
    } else {
       calculateOneLJ14Ixn( inverseR, sig, eps, &dEdR, &energy );
    }
    if (cutoff)
       dEdR += parameters[2]*(inverseR-2.0f*krf*r2);
    else
       dEdR += parameters[2]*inverseR;
    dEdR     *= inverseR*inverseR;

   // accumulate forces

   for( int ii = 0; ii < 3; ii++ ){
      RealOpenMM force        = dEdR*deltaR[0][ii];
      forces[atomAIndex][ii] += force;
      forces[atomBIndex][ii] -= force;
   }

   if (cutoff)
       energy += parameters[2]*(inverseR+krf*r2-crf);
   else
       energy += parameters[2]*inverseR;

   // accumulate energies

   if (totalEnergy != NULL)
       *totalEnergy += energy;
}

/**---------------------------------------------------------------------------------------

     Calculate LJ pair ixn between two atoms

     @param inverseR         1/r
     @param sig              sigma
     @param eps              epsilon
     @param dEdR             output force factor
     @param energy           LJ energy

     --------------------------------------------------------------------------------------- */

void ReferenceFreeEnergyLJCoulomb14Softcore::calculateOneLJ14Ixn( RealOpenMM inverseR, RealOpenMM sig, RealOpenMM eps,
                                                                  RealOpenMM* dEdR, RealOpenMM* energy ) const {

    // ---------------------------------------------------------------------------------------

    //static const std::string methodName = "\nReferenceLJ14CoulombIxn::calculateOneLJIxn";

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
}

  /**---------------------------------------------------------------------------------------

     Calculate softcore LJ pair ixn between two atoms

     @param r                r
     @param sig              sigma
     @param eps              epsilon
     @param lambda           lambda
     @param dEdR             output force factor
     @param energy           LJ energy

     --------------------------------------------------------------------------------------- */

void ReferenceFreeEnergyLJCoulomb14Softcore::calculateOneSoftCoreLJ14Ixn( RealOpenMM r, RealOpenMM sig, RealOpenMM eps,
                                                                         RealOpenMM lambda,
                                                                         RealOpenMM* dEdR, RealOpenMM* energy ) const {

    // ---------------------------------------------------------------------------------------

    //static const std::string methodName = "\nReferenceFreeEnergyLJCoulomb14Softcore::calculateOneSoftCoreLJ14Ixn";

    // ---------------------------------------------------------------------------------------

    // constants -- reduce Visual Studio warnings regarding conversions between float & double

    static const RealOpenMM zero        =  0.0;
    static const RealOpenMM one         =  1.0;
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
}

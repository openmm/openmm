
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
#include "ReferenceLJCoulomb14.h"
#include "ReferenceForce.h"

/**---------------------------------------------------------------------------------------

   ReferenceLJCoulomb14 constructor

   --------------------------------------------------------------------------------------- */

ReferenceLJCoulomb14::ReferenceLJCoulomb14( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceLJCoulomb14::ReferenceLJCoulomb14";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   ReferenceLJCoulomb14 destructor

   --------------------------------------------------------------------------------------- */

ReferenceLJCoulomb14::~ReferenceLJCoulomb14( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceLJCoulomb14::~ReferenceLJCoulomb14";

   // ---------------------------------------------------------------------------------------

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

   @return ReferenceForce::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceLJCoulomb14::getDerivedParameters( RealOpenMM c6, RealOpenMM c12, RealOpenMM q1,
                                         RealOpenMM q2, RealOpenMM epsfac,
                                         RealOpenMM* parameters ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceLJCoulomb14::getDerivedParameters";

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

   return ReferenceForce::DefaultReturn;
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
   @param energiesByBond   energies by bond: energiesByBond[bondIndex]
   @param energiesByAtom   energies by atom: energiesByAtom[atomIndex]

   @return ReferenceForce::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceLJCoulomb14::calculateBondIxn( int* atomIndices, RealOpenMM** atomCoordinates,
                                     RealOpenMM* parameters, RealOpenMM** forces,
                                     RealOpenMM* energiesByBond,
                                     RealOpenMM* energiesByAtom ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceLJCoulomb14::calculateBondIxn";

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nReferenceLJCoulomb14::calculateBondIxn";

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

   RealOpenMM inverseR  = one/(deltaR[0][ReferenceForce::RIndex]);
   RealOpenMM sig2      = inverseR*parameters[0];
              sig2     *= sig2;
   RealOpenMM sig6      = sig2*sig2*sig2;

   RealOpenMM dEdR      = parameters[1]*( twelve*sig6 - six )*sig6;
              dEdR     += parameters[2]*inverseR;
              dEdR     *= inverseR*inverseR;

   // accumulate forces

   for( int ii = 0; ii < 3; ii++ ){
      RealOpenMM force        = dEdR*deltaR[0][ii];
      forces[atomAIndex][ii] += force;
      forces[atomBIndex][ii] -= force;
   }

   RealOpenMM energy = parameters[1]*( sig6 - one )*sig6 + parameters[2]*inverseR;

   // accumulate energies

   updateEnergy( energy, energiesByBond, LastAtomIndex, atomIndices, energiesByAtom );

   // debug 

   if( debug ){
      static bool printHeader = false;
      std::stringstream message;
      message << methodName;
      message << std::endl;
      if( !printHeader  ){  
         printHeader = true;
         message << std::endl;
         message << methodName.c_str() << " a0 k [c q p s] r1 r2  angle dt rp p[] dot cosine angle dEdR*r F[]" << std::endl;
      }   

      message << std::endl;
      for( int ii = 0; ii < LastAtomIndex; ii++ ){
         message << " Atm " << atomIndices[ii] << " [" << atomCoordinates[atomIndices[ii]][0] << " " << atomCoordinates[atomIndices[ii]][1] << "] ";
      }
      message << std::endl << " Delta:";
      for( int ii = 0; ii < (LastAtomIndex - 1); ii++ ){
         message << " [";
         for( int jj = 0; jj < ReferenceForce::LastDeltaRIndex; jj++ ){
            message << deltaR[ii][jj] << " ";
         }
         message << "]";
      }
      message << std::endl;

      message << " p1="     << parameters[0];
      message << " p2="     << parameters[1];
      message << " p3="     << parameters[2];
      message << std::endl << "  ";

      message << " dEdR=" << dEdR;
      message << " E=" << energy << " force factors: ";
      message << "F=compute force; f=cumulative force";

      message << std::endl << "  ";
      for( int ii = 0; ii < LastAtomIndex; ii++ ){
         message << " F" << (ii+1) << "[";
         SimTKOpenMMUtilities::formatRealStringStream( message, deltaR[0], threeI, dEdR );
         message << "]";
      }   
      message << std::endl << "  ";

      for( int ii = 0; ii < LastAtomIndex; ii++ ){
         message << " f" << (ii+1) << "[";
         SimTKOpenMMUtilities::formatRealStringStream( message, forces[atomIndices[ii]], threeI );
         message << "]";
      }

      SimTKOpenMMLog::printMessage( message );
   }   

   return ReferenceForce::DefaultReturn;
}

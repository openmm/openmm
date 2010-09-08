
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

#include "AmoebaReferenceHarmonicBondForce.h"
#include <vector>

/**---------------------------------------------------------------------------------------

   Calculate Amoeba harmonic bond ixn (force and energy)

   @param positionAtomA           Cartesian coordinates of atom A
   @param positionAtomB           Cartesian coordinates of atom B
   @param bondLength              bond length
   @param bondK                   bond force
   @param bondCubic               cubic bond force parameter
   @param bondQuartic             quartic bond force parameter
   @param forces                  force vector

   @return energy

   --------------------------------------------------------------------------------------- */

RealOpenMM AmoebaReferenceHarmonicBondForce::calculateForceAndEnergy( RealOpenMM* positionAtomA, RealOpenMM* positionAtomB,
                                                                      RealOpenMM bondLength, RealOpenMM bondK,
                                                                      RealOpenMM bondCubic, RealOpenMM bondQuartic,
                                                                      RealOpenMM** forces ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "AmoebaReferenceHarmonicBondForce::calculateHarmonicForce";

   static const RealOpenMM zero          = 0.0;
   static const RealOpenMM one           = 1.0;
   static const RealOpenMM onePt5        = 1.5;
   static const RealOpenMM two           = 2.0;

   // ---------------------------------------------------------------------------------------

   // get deltaR, R2, and R between 2 atoms

   std::vector<RealOpenMM> deltaR;
   deltaR.push_back( positionAtomB[0] - positionAtomA[0] );
   deltaR.push_back( positionAtomB[1] - positionAtomA[1] );
   deltaR.push_back( positionAtomB[2] - positionAtomA[2] );

   RealOpenMM r               = SQRT( deltaR[0]*deltaR[0] + deltaR[1]*deltaR[1] + deltaR[2]*deltaR[2] );

   // deltaIdeal = r - r_0

   RealOpenMM deltaIdeal      = r - bondLength;
   RealOpenMM deltaIdeal2     = deltaIdeal*deltaIdeal;

   RealOpenMM dEdR            = (one + onePt5*bondCubic*deltaIdeal + two*bondQuartic*deltaIdeal2);
   dEdR                      *= two*bondK*deltaIdeal;
   dEdR                       = r > zero ? (dEdR/r) : zero;

   forces[0][0]              += dEdR*deltaR[0];
   forces[0][1]              += dEdR*deltaR[1];
   forces[0][2]              += dEdR*deltaR[2];

   forces[1][0]              -= dEdR*deltaR[0];
   forces[1][1]              -= dEdR*deltaR[1];
   forces[1][2]              -= dEdR*deltaR[2];

   RealOpenMM energy          = bondK*deltaIdeal2*( one + bondCubic*deltaIdeal + bondQuartic*deltaIdeal2 );
   return energy;
}


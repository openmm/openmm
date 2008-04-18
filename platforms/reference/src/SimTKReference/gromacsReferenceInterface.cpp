
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

// #include "SimTKOpenMMGromacsUtilities.h"

#include "../SimTKUtilities/SimTKOpenMMGromacsUtilities.h"
#include "gromacsReferenceInterface.h"
#include "gromacsReferenceInterfaceCpp.h"

#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"

#include "ReferenceBondForce.h"
#include "ReferenceHarmonicBondIxn.h"
#include "ReferenceAngleBondIxn.h"
#include "ReferenceProperDihedralBond.h"
#include "ReferenceRbDihedralBond.h"
#include "ReferenceLJ14.h"
#include "ReferenceLJCoulombIxn.h"

#include "ReferenceStochasticDynamics.h"
#include "ReferenceConstraint.h"

/**---------------------------------------------------------------------------------------

   Allocate memory for atom indices participating in bonds and the bond parameters

   @param numberOfBonds            number of bonds
   @param numberOfAtomIndices      (number of atoms)/bond
   @param outputBondIndices        allocated memory for atom indices outputBondIndices[bondIndex][atomIndex]
   @param numberParameters         number of parameters/bond
   @param outputBondParameters     allocated memory for parameters outputBondParameters[bondIndex][parameterIndex]

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsAllocateBondIxnArrays( int numberOfBonds, int numberOfAtomIndices,
                                  int*** outputBondIndices, int numberParameters,
                                      RealOpenMM*** outputBondParameters ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "\ngromacsAllocateBondIxnArrays: ";

   // ---------------------------------------------------------------------------------------

   int** atomIndices                   = (int**) malloc( sizeof( int* )*numberOfBonds );
   int* atomIndicesBlock               = (int*)  malloc( sizeof( int  )*numberOfBonds*numberOfAtomIndices );
   memset( atomIndicesBlock, 0, sizeof( int  )*numberOfBonds*numberOfAtomIndices );

   RealOpenMM** bondParameters         = (RealOpenMM**) malloc( sizeof( RealOpenMM* )*numberOfBonds );
   RealOpenMM* bondParametersBlock     = (RealOpenMM*)  malloc( sizeof( RealOpenMM )*numberOfBonds*numberParameters );
   memset( bondParametersBlock, 0, sizeof( RealOpenMM )*numberOfBonds*numberParameters );

   for( int ii = 0; ii < numberOfBonds; ii++ ){ 

      atomIndices[ii]       = atomIndicesBlock;
      atomIndicesBlock     += numberOfAtomIndices;

      bondParameters[ii]    = bondParametersBlock;
      bondParametersBlock  += numberParameters;

   }

   *outputBondIndices    = atomIndices;
   *outputBondParameters = bondParameters;

   return 0;
}

/**---------------------------------------------------------------------------------------

   Free memory for arrays

   @param indexOneDArraysToFree    number of 1D RealOpenMM arrays
   @param indexTwoDArraysToFree    number of 2D RealOpenMM arrays
   @param indexTwoDIntArraysToFree number of 2D int arrays
   @param MaxFreeArrayIndex        max index on arrays
   @param oneDArraysToFree         1D RealOpenMM arrays to free
   @param twoDArraysToFree         2D RealOpenMM arrays to free
   @param twoDIntArraysToFree      2D int arrays to free
   @param whoIsCalling             callee

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsFreeSundryArrays( int indexOneDArraysToFree, int indexTwoDArraysToFree,
                             int indexTwoDIntArraysToFree, int MaxFreeArrayIndex,
                             RealOpenMM** oneDArraysToFree,
                             RealOpenMM*** twoDArraysToFree, int*** twoDIntArraysToFree,
                             const char* whoIsCalling ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "\ngromacsFreeSundryArrays ";

   // ---------------------------------------------------------------------------------------

   // ---------------------------------------------------------------------------------------

   // free memory

   if( indexTwoDArraysToFree    >= MaxFreeArrayIndex ||
       indexTwoDIntArraysToFree >= MaxFreeArrayIndex ||
       indexOneDArraysToFree    >= MaxFreeArrayIndex ){

      std::stringstream message;
      message << methodName << " " << whoIsCalling;
      message << " Free array index too large: [" << indexTwoDArraysToFree << " " << indexTwoDIntArraysToFree << " " << indexOneDArraysToFree << "]";
      message << " should be less than " << MaxFreeArrayIndex;

      SimTKOpenMMLog::printMessage( message );

      indexTwoDArraysToFree     = indexTwoDArraysToFree    >= MaxFreeArrayIndex ? MaxFreeArrayIndex : indexTwoDArraysToFree;
      indexTwoDIntArraysToFree  = indexTwoDIntArraysToFree >= MaxFreeArrayIndex ? MaxFreeArrayIndex : indexTwoDIntArraysToFree;
      indexOneDArraysToFree     = indexOneDArraysToFree    >= MaxFreeArrayIndex ? MaxFreeArrayIndex : indexOneDArraysToFree;
   }

   for( int ii = 0; ii < indexTwoDArraysToFree; ii++ ){
      free( twoDArraysToFree[ii][0] );
      free(  twoDArraysToFree[ii] );
   }

   for( int ii = 0; ii < indexTwoDIntArraysToFree; ii++ ){
      free( twoDIntArraysToFree[ii][0] );
      free(  twoDIntArraysToFree[ii] );
   }

   for( int ii = 0; ii < indexOneDArraysToFree; ii++ ){
      free( oneDArraysToFree[ii] );
   }

   return 0;
}

/**---------------------------------------------------------------------------------------

   Get bond atom indices and parameters for harmonic bonds

   @param top                      Gromacs t_topology data struct
   @param atomIndices              array of atomIndices[bondIndex][atomIndex]
   @param bondParameters           array of bond parameters[bondIndex][parameterIndex]

   The arrays atomIndices & bondParameters are allocated off the heap;
   The memory for the second index is one block of size (number of bonds)*(number of atoms in bond)

   The memory should be freed by the callee via calls:
		free( atomIndices[0] ); free( atomIndices );

   @return number of bonds

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsGetHarmonicBondParameters( const t_topology* top, int*** atomIndices,
                                      RealOpenMM*** bondParameters ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "\ngromacsGetHarmonicBondParameters: ";
   static const int twoI           = 2;
   static const RealOpenMM two     = 2.0;

   // ---------------------------------------------------------------------------------------

   const t_idef* idef   = &top->idef;
   int numberOfBonds    = idef->il[ F_BONDS ].nr/3; 
   int* atoms           = (int*) idef->il[ F_BONDS ].iatoms;

   gromacsAllocateBondIxnArrays( numberOfBonds, twoI, atomIndices, twoI, bondParameters );

   int offsetIndex = 0;
   for( int ii = 0; ii < numberOfBonds; ii++ ) { 

      int type                 = atoms[ offsetIndex++ ];  
      (*atomIndices)[ii][0]    = atoms[ offsetIndex++ ];
      (*atomIndices)[ii][1]    = atoms[ offsetIndex++ ];

      (*bondParameters)[ii][0] = idef->iparams[ type ].harmonic.rA;
      (*bondParameters)[ii][1] = idef->iparams[ type ].harmonic.krA;
   }

   return numberOfBonds;
}

/**---------------------------------------------------------------------------------------

   Calculate harmonic bond forces/energies for given configuration and topology

   @param top                      Gromacs t_topology data struct
   @param atomCoordinates          atom configuration
   @param harmonicBondForces       output array of forces: harmonicBondForces[atomIndex][3]
   @param harmonicBondEnergies     output array of bond energies[atomIndex]
   @param totalEnergy              output total energy

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsHarmonicBondReferenceForce( const t_topology* top, RealOpenMM** atomCoordinates,
                                       RealOpenMM** harmonicBondForces,
                                       RealOpenMM* harmonicBondAtomEnergies,
                                       RealOpenMM* totalEnergy ){ 

   // ---------------------------------------------------------------------------------------

   static const char* methodName         = "\ngromacsHarmonicReferenceForce: ";
   static int debug                      = true;
   static int startAtom                  = 0;
   static int stopAtom                   = 3;

   static const RealOpenMM zero          = 0.0;
   static const RealOpenMM one           = 1.0;

   static const int  MaxFreeArrayIndex   = 20;

   int indexOneDArraysToFree             = 0;
   RealOpenMM* oneDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDArraysToFree             = 0;
   RealOpenMM** twoDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDIntArraysToFree          = 0;
   int** twoDIntArraysToFree[MaxFreeArrayIndex];

   // ---------------------------------------------------------------------------------------

   int numberOfAtoms = top->atoms.nr;

   // ---------------------------------------------------------------------------------------

   // get bond atom indices & parameters

   int** harmonicBondAtomIndices;
   RealOpenMM** harmonicBondParameters;
   int numberOfBonds                                 = gromacsGetHarmonicBondParameters( top, &harmonicBondAtomIndices, &harmonicBondParameters );
   twoDIntArraysToFree[indexTwoDIntArraysToFree++]   = harmonicBondAtomIndices;
   twoDArraysToFree[indexTwoDArraysToFree++]         = harmonicBondParameters;

   // bond energy array

   RealOpenMM* harmonicBondBondEnergies              = (RealOpenMM*) malloc( sizeof( RealOpenMM )*numberOfBonds );
   memset( harmonicBondBondEnergies, 0, sizeof( RealOpenMM )*numberOfBonds );
   oneDArraysToFree[indexOneDArraysToFree++]         = harmonicBondBondEnergies;

   ReferenceHarmonicBondIxn referenceHarmonicBondIxn;
   ReferenceBondForce referenceBondForce;

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << " harmonic bonds=" << numberOfBonds;
      SimTKOpenMMLog::printMessage( message );
   }

   // harmonic bond forces

   referenceBondForce.calculateForce( numberOfBonds, harmonicBondAtomIndices, atomCoordinates,
                                      harmonicBondParameters, harmonicBondForces,
                                      harmonicBondBondEnergies, harmonicBondAtomEnergies, 
                                      totalEnergy, referenceHarmonicBondIxn );

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << " Harmonic Total energy=" << totalEnergy << std::endl;
      for( int ii = startAtom; ii < stopAtom; ii++ ){
         message << ii << " F[";
         SimTKOpenMMUtilities::formatRealStringStream( message, harmonicBondForces[ii], 3, one );
         message << "] E=" << harmonicBondAtomEnergies[ii] << std::endl;
      } 
      SimTKOpenMMLog::printMessage( message );
   }

   // ---------------------------------------------------------------------------------------

   // free memory

   gromacsFreeSundryArrays( indexOneDArraysToFree, indexTwoDArraysToFree, indexTwoDIntArraysToFree,
                            MaxFreeArrayIndex, oneDArraysToFree, twoDArraysToFree, twoDIntArraysToFree,
                            methodName );

   return 0;
}

/**---------------------------------------------------------------------------------------

   Get bond atom indices and parameters for angle bonds

   @param top                      Gromacs t_topology data struct
   @param atomIndices              array of atomIndices[bondIndex][atomIndex]
   @param bondParameters           array of bond parameters[bondIndex][parameterIndex]

   The arrays atomIndices & bondParameters are allocated off the heap;
   The memory for the second index is one block of size (number of bonds)*(number of atoms in bond)

   The memory should be freed by the callee via calls:
		free( atomIndices[0] ); free( atomIndices );

   @return number of bonds

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsGetAngleBondParameters( const t_topology* top, int*** atomIndices,
                                   RealOpenMM*** bondParameters ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "\ngromacsGetAngleBondParameters: ";
   static const int twoI           = 2;
   static const int threeI         = 3;
   static const RealOpenMM two     = 2.0;

   // ---------------------------------------------------------------------------------------

   const t_idef* idef   = &top->idef;
   int numberOfBonds    = idef->il[ F_ANGLES ].nr/4; 
   int* atoms           = (int*) idef->il[ F_ANGLES ].iatoms;

   gromacsAllocateBondIxnArrays( numberOfBonds, threeI, atomIndices, twoI, bondParameters );

   int offsetIndex = 0;
   for( int ii = 0; ii < numberOfBonds; ii++ ) { 

      int type                 = atoms[ offsetIndex++ ];  
      (*atomIndices)[ii][0]    = atoms[ offsetIndex++ ];
      (*atomIndices)[ii][1]    = atoms[ offsetIndex++ ];
      (*atomIndices)[ii][2]    = atoms[ offsetIndex++ ];

      (*bondParameters)[ii][0] = idef->iparams[ type ].harmonic.rA;
      (*bondParameters)[ii][1] = idef->iparams[ type ].harmonic.krA;
   }

   return numberOfBonds;
}

/**---------------------------------------------------------------------------------------

   Calculate angle bond forces/energies for given configuration and topology

   @param top                      Gromacs t_topology data struct
   @param atomCoordinates          atom configuration
   @param angleBondForces       output array of forces: angleBondForces[atomIndex][3]
   @param angleBondEnergies     output array of bond energies[atomIndex]
   @param totalEnergy              output total energy

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsAngleBondReferenceForce( const t_topology* top, RealOpenMM** atomCoordinates,
                                    RealOpenMM** angleBondForces,
                                    RealOpenMM* angleBondAtomEnergies,
                                    RealOpenMM* totalEnergy ){ 

   // ---------------------------------------------------------------------------------------

   static const char* methodName         = "\ngromacsAngleReferenceForce: ";
   static int debug                      = true;
   static int startAtom                  = 0;
   static int stopAtom                   = 3;

   static const RealOpenMM zero          = 0.0;
   static const RealOpenMM one           = 1.0;

   static const int  MaxFreeArrayIndex   = 20;

   int indexOneDArraysToFree             = 0;
   RealOpenMM* oneDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDArraysToFree             = 0;
   RealOpenMM** twoDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDIntArraysToFree          = 0;
   int** twoDIntArraysToFree[MaxFreeArrayIndex];

   // ---------------------------------------------------------------------------------------

   int numberOfAtoms = top->atoms.nr;

   // ---------------------------------------------------------------------------------------

   // get bond atom indices & parameters

   int** angleBondAtomIndices;
   RealOpenMM** angleBondParameters;
   int numberOfBonds                                 = gromacsGetAngleBondParameters( top, &angleBondAtomIndices, &angleBondParameters );
   twoDIntArraysToFree[indexTwoDIntArraysToFree++]   = angleBondAtomIndices;
   twoDArraysToFree[indexTwoDArraysToFree++]         = angleBondParameters;

   // bond energy array

   RealOpenMM* angleBondBondEnergies                 = (RealOpenMM*) malloc( sizeof( RealOpenMM )*numberOfBonds );
   memset( angleBondBondEnergies, 0, sizeof( RealOpenMM )*numberOfBonds );
   oneDArraysToFree[indexOneDArraysToFree++]         = angleBondBondEnergies;

   ReferenceAngleBondIxn referenceAngleBondIxn;
   ReferenceBondForce referenceBondForce;

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << " angle bonds=" << numberOfBonds;
      SimTKOpenMMLog::printMessage( message );
   }

   // angle bond forces

   referenceBondForce.calculateForce( numberOfBonds, angleBondAtomIndices, atomCoordinates,
                                      angleBondParameters, angleBondForces, 
                                      angleBondBondEnergies, angleBondAtomEnergies, 
                                      totalEnergy, referenceAngleBondIxn );

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << "Angle  Total energy=" << *totalEnergy << std::endl;
      for( int ii = startAtom; ii < stopAtom; ii++ ){
         message << ii << " F[";
         SimTKOpenMMUtilities::formatRealStringStream( message, angleBondForces[ii], 3, one );
         message << "] E=" << angleBondAtomEnergies[ii] << std::endl;
      } 
      SimTKOpenMMLog::printMessage( message );
   }

   // ---------------------------------------------------------------------------------------

   // free memory

   gromacsFreeSundryArrays( indexOneDArraysToFree, indexTwoDArraysToFree, indexTwoDIntArraysToFree,
                            MaxFreeArrayIndex, oneDArraysToFree, twoDArraysToFree, twoDIntArraysToFree,
                            methodName );

   // ---------------------------------------------------------------------------------------

   return 0;
}

/**---------------------------------------------------------------------------------------

   Get bond atom indices and parameters for properDihedral bonds

   @param top                      Gromacs t_topology data struct
   @param atomIndices              array of atomIndices[bondIndex][atomIndex]
   @param bondParameters           array of bond parameters[bondIndex][parameterIndex]

   The arrays atomIndices & bondParameters are allocated off the heap;
   The memory for the second index is one block of size (number of bonds)*(number of atoms in bond)

   The memory should be freed by the callee via calls:
		free( atomIndices[0] ); free( atomIndices );

   @return number of bonds

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsGetProperDihedralBondParameters( const t_topology* top, int*** atomIndices,
                                            RealOpenMM*** bondParameters ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "\ngromacsGetProperDihedralBondParameters: ";
   static const int twoI           = 2;
   static const int threeI         = 3;
   static const int fourI          = 4;
   static const RealOpenMM two     = 2.0;

   // ---------------------------------------------------------------------------------------

   const t_idef* idef   = &top->idef;
   int numberOfBonds    = idef->il[ F_PDIHS ].nr/5; 
   int* atoms           = (int*) idef->il[ F_PDIHS ].iatoms;

   gromacsAllocateBondIxnArrays( numberOfBonds, fourI, atomIndices, threeI, bondParameters );

   int offsetIndex = 0;
   for( int ii = 0; ii < numberOfBonds; ii++ ) { 

      int type                 = atoms[ offsetIndex++ ];  
      (*atomIndices)[ii][0]    = atoms[ offsetIndex++ ];
      (*atomIndices)[ii][1]    = atoms[ offsetIndex++ ];
      (*atomIndices)[ii][2]    = atoms[ offsetIndex++ ];
      (*atomIndices)[ii][3]    = atoms[ offsetIndex++ ];

      (*bondParameters)[ii][0] = idef->iparams[ type ].pdihs.cpA;
      (*bondParameters)[ii][1] = idef->iparams[ type ].pdihs.phiA;
      (*bondParameters)[ii][2] = (RealOpenMM) idef->iparams[ type ].pdihs.mult;
   }

   return numberOfBonds;
}

/**---------------------------------------------------------------------------------------

   Calculate properDihedral bond forces/energies for given configuration and topology

   @param top                      Gromacs t_topology data struct
   @param atomCoordinates          atom configuration
   @param properDihedralBondForces       output array of forces: angleBondForces[atomIndex][3]
   @param properDihedralBondEnergies     output array of bond energies[atomIndex]
   @param totalEnergy              output total energy

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsProperDihedralBondReferenceForce( const t_topology* top, RealOpenMM** atomCoordinates,
                                             RealOpenMM** properDihedralBondForces,
                                             RealOpenMM* properDihedralBondAtomEnergies,
                                             RealOpenMM* totalEnergy ){ 

   // ---------------------------------------------------------------------------------------

   static const char* methodName         = "\ngromacsProperDihedralReferenceForce: ";
   static int debug                      = true;
   static int startAtom                  = 0;
   static int stopAtom                   = 3;

   static const RealOpenMM zero          = 0.0;
   static const RealOpenMM one           = 1.0;

   static const int  MaxFreeArrayIndex   = 20;

   int indexOneDArraysToFree             = 0;
   RealOpenMM* oneDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDArraysToFree             = 0;
   RealOpenMM** twoDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDIntArraysToFree          = 0;
   int** twoDIntArraysToFree[MaxFreeArrayIndex];

   // ---------------------------------------------------------------------------------------

   int numberOfAtoms = top->atoms.nr;

   // ---------------------------------------------------------------------------------------

   // get bond atom indices & parameters

   int** properDihedralBondAtomIndices;
   RealOpenMM** properDihedralBondParameters;
   int numberOfBonds                                 = gromacsGetProperDihedralBondParameters( top, &properDihedralBondAtomIndices, &properDihedralBondParameters );
   twoDIntArraysToFree[indexTwoDIntArraysToFree++]   = properDihedralBondAtomIndices;
   twoDArraysToFree[indexTwoDArraysToFree++]         = properDihedralBondParameters;

   // bond energy array

   RealOpenMM* properDihedralBondBondEnergies        = (RealOpenMM*) malloc( sizeof( RealOpenMM )*numberOfBonds );
   memset( properDihedralBondBondEnergies, 0, sizeof( RealOpenMM )*numberOfBonds );
   oneDArraysToFree[indexOneDArraysToFree++]         = properDihedralBondBondEnergies;

   ReferenceProperDihedralBond referenceProperDihedralBond;
   ReferenceBondForce referenceBondForce;

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << " properDihedral bonds=" << numberOfBonds << std::endl;
/*
      for( int ii  = 0; ii < numberOfBonds; ii++ ){
         message << "   bnd=" << ii << " a[";
         for( int jj  = 0; jj < 4; jj++ ){
            message << properDihedralBondAtomIndices[ii][jj] << " ";
         }
         message << "] p[";
         for( int jj  = 0; jj < 3; jj++ ){
            message << properDihedralBondParameters[ii][jj] << " ";
         }
         message << "]" << std::endl; 
      }
*/
      SimTKOpenMMLog::printMessage( message );
   }

   // properDihedral bond forces

   referenceBondForce.calculateForce( numberOfBonds, properDihedralBondAtomIndices, atomCoordinates,
                                      properDihedralBondParameters, properDihedralBondForces, 
                                      properDihedralBondBondEnergies, properDihedralBondAtomEnergies, 
                                      totalEnergy, referenceProperDihedralBond);

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << "ProperDihedral  Total energy=" << *totalEnergy << std::endl;
      for( int ii = startAtom; ii < stopAtom; ii++ ){
         message << ii << " F[";
         SimTKOpenMMUtilities::formatRealStringStream( message, properDihedralBondForces[ii], 3, one );
         message << "] E=" << properDihedralBondAtomEnergies[ii] << std::endl;
      } 
      SimTKOpenMMLog::printMessage( message );
   }

   // ---------------------------------------------------------------------------------------

   // free memory

   gromacsFreeSundryArrays( indexOneDArraysToFree, indexTwoDArraysToFree, indexTwoDIntArraysToFree,
                            MaxFreeArrayIndex, oneDArraysToFree, twoDArraysToFree, twoDIntArraysToFree,
                            methodName );

   // ---------------------------------------------------------------------------------------

   return 0;
}

/**---------------------------------------------------------------------------------------

   Get bond atom indices and parameters for RB Dihedral bonds

   @param top                      Gromacs t_topology data struct
   @param atomIndices              array of atomIndices[bondIndex][atomIndex]
   @param bondParameters           array of bond parameters[bondIndex][parameterIndex]

   The arrays atomIndices & bondParameters are allocated off the heap;
   The memory for the second index is one block of size (number of bonds)*(number of atoms in bond)

   The memory should be freed by the callee via calls:
		free( atomIndices[0] ); free( atomIndices );

   @return number of bonds

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsGetRbDihedralBondParameters( const t_topology* top, int*** atomIndices,
                                            RealOpenMM*** bondParameters ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "\ngromacsGetRbDihedralBondParameters: ";
   static const int twoI           = 2;
   static const int threeI         = 3;
   static const int fourI          = 4;
   static const int sixI           = 6;

   static const RealOpenMM two     = 2.0;

   // ---------------------------------------------------------------------------------------

   const t_idef* idef   = &top->idef;
   int numberOfBonds    = idef->il[ F_RBDIHS ].nr/5; 
   int* atoms           = (int*) idef->il[ F_RBDIHS ].iatoms;

   gromacsAllocateBondIxnArrays( numberOfBonds, fourI, atomIndices, sixI, bondParameters );

   int offsetIndex = 0;
   for( int ii = 0; ii < numberOfBonds; ii++ ) { 

      int type                 = atoms[ offsetIndex++ ];  
      (*atomIndices)[ii][0]    = atoms[ offsetIndex++ ];
      (*atomIndices)[ii][1]    = atoms[ offsetIndex++ ];
      (*atomIndices)[ii][2]    = atoms[ offsetIndex++ ];
      (*atomIndices)[ii][3]    = atoms[ offsetIndex++ ];

      for( int jj = 0; jj < sixI; jj++ ) { 
         (*bondParameters)[ii][jj] = idef->iparams[ type ].rbdihs.rbc[jj];
      }
   }

   return numberOfBonds;
}

/**---------------------------------------------------------------------------------------

   Calculate RB dihedral bond forces/energies for given configuration and topology

   @param top                      Gromacs t_topology data struct
   @param atomCoordinates          atom configuration
   @param rbDihedralBondForces     output array of forces: angleBondForces[atomIndex][3]
   @param rbDihedralBondEnergies   output array of bond energies[atomIndex]
   @param totalEnergy              output total energy

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsRbDihedralBondReferenceForce( const t_topology* top, RealOpenMM** atomCoordinates,
                                         RealOpenMM** rbDihedralBondForces,
                                         RealOpenMM* rbDihedralBondAtomEnergies,
                                         RealOpenMM* totalEnergy ){ 

   // ---------------------------------------------------------------------------------------

   static const char* methodName         = "\ngromacsRbDihedralReferenceForce: ";
   static int debug                      = true;
   static int startAtom                  = 0;
   static int stopAtom                   = 3;

   static const RealOpenMM zero          = 0.0;
   static const RealOpenMM one           = 1.0;

   static const int  MaxFreeArrayIndex   = 20;

   int indexOneDArraysToFree             = 0;
   RealOpenMM* oneDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDArraysToFree             = 0;
   RealOpenMM** twoDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDIntArraysToFree          = 0;
   int** twoDIntArraysToFree[MaxFreeArrayIndex];

   // ---------------------------------------------------------------------------------------

   int numberOfAtoms = top->atoms.nr;

   // ---------------------------------------------------------------------------------------

   // get bond atom indices & parameters

   int** rbDihedralBondAtomIndices;
   RealOpenMM** rbDihedralBondParameters;
   int numberOfBonds                                 = gromacsGetRbDihedralBondParameters( top, &rbDihedralBondAtomIndices, &rbDihedralBondParameters );
   twoDIntArraysToFree[indexTwoDIntArraysToFree++]   = rbDihedralBondAtomIndices;
   twoDArraysToFree[indexTwoDArraysToFree++]         = rbDihedralBondParameters;

   // bond energy array

   RealOpenMM* rbDihedralBondBondEnergies            = (RealOpenMM*) malloc( sizeof( RealOpenMM )*numberOfBonds );
   memset( rbDihedralBondBondEnergies, 0, sizeof( RealOpenMM )*numberOfBonds );
   oneDArraysToFree[indexOneDArraysToFree++]         = rbDihedralBondBondEnergies;

   ReferenceRbDihedralBond referenceRbDihedralBond;
   ReferenceBondForce referenceBondForce;

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << " rbDihedral bonds=" << numberOfBonds << std::endl;
/*
      for( int ii  = 0; ii < numberOfBonds; ii++ ){
         message << "   bnd=" << ii << " a[";
         for( int jj  = 0; jj < 4; jj++ ){
            message << rbDihedralBondAtomIndices[ii][jj] << " ";
         }
         message << "] p[";
         for( int jj  = 0; jj < 3; jj++ ){
            message << rbDihedralBondParameters[ii][jj] << " ";
         }
         message << "]" << std::endl; 
      }
*/
      SimTKOpenMMLog::printMessage( message );
   }

   // rbDihedral bond forces

   referenceBondForce.calculateForce( numberOfBonds, rbDihedralBondAtomIndices, atomCoordinates,
                                      rbDihedralBondParameters, rbDihedralBondForces, 
                                      rbDihedralBondBondEnergies, rbDihedralBondAtomEnergies, 
                                      totalEnergy, referenceRbDihedralBond);

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << "RbDihedral  Total energy=" << *totalEnergy << std::endl;
      for( int ii = startAtom; ii < stopAtom; ii++ ){
         message << ii << " F[";
         SimTKOpenMMUtilities::formatRealStringStream( message, rbDihedralBondForces[ii], 3, one );
         message << "] E=" << rbDihedralBondAtomEnergies[ii] << std::endl;
      } 
      SimTKOpenMMLog::printMessage( message );
   }

   // ---------------------------------------------------------------------------------------

   // free memory

   gromacsFreeSundryArrays( indexOneDArraysToFree, indexTwoDArraysToFree, indexTwoDIntArraysToFree,
                            MaxFreeArrayIndex, oneDArraysToFree, twoDArraysToFree, twoDIntArraysToFree,
                            methodName );

   // ---------------------------------------------------------------------------------------

   return 0;
}

/**---------------------------------------------------------------------------------------

   Get bond atom indices and parameters for LJ 14 bonds

   @param top                      Gromacs t_topology data struct
   @param atomIndices              array of atomIndices[bondIndex][atomIndex]
   @param bondParameters           array of bond parameters[bondIndex][parameterIndex]

   The arrays atomIndices & bondParameters are allocated off the heap;
   The memory for the second index is one block of size (number of bonds)*(number of atoms in bond)

   The memory should be freed by the callee via calls:
		free( atomIndices[0] ); free( atomIndices );

   @return number of bonds

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsGetLJ14Parameters( const t_topology* top, int*** atomIndices,
                              RealOpenMM*** bondParameters ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "\ngromacsGetLJ14Parameters: ";
   static const int twoI           = 2;
   static const int threeI         = 3;
   static const int fourI          = 4;
   static const int sixI           = 6;

   static const RealOpenMM two     = 2.0;

   // ---------------------------------------------------------------------------------------

   const t_idef* idef   = &top->idef;
   int numberOfBonds    = idef->il[ F_LJ14 ].nr/3; 
   int* atoms           = (int*) idef->il[ F_LJ14 ].iatoms;

   gromacsAllocateBondIxnArrays( numberOfBonds, twoI, atomIndices, threeI, bondParameters );

   int offsetIndex = 0;
   for( int ii = 0; ii < numberOfBonds; ii++ ) { 

      int type                 = atoms[ offsetIndex++ ];  
      (*atomIndices)[ii][0]    = atoms[ offsetIndex++ ];
      (*atomIndices)[ii][1]    = atoms[ offsetIndex++ ];

      (*bondParameters)[ii][0] = idef->iparams[ type ].lj14.c6A;
      (*bondParameters)[ii][1] = idef->iparams[ type ].lj14.c12A;
   }

   return numberOfBonds;
}

/**---------------------------------------------------------------------------------------

   Calculate LJ 14 forces/energies for given configuration and topology

   @param top                      Gromacs t_topology data struct
   @param atomCoordinates          atom configuration
   @param fr                       Gromacs t_forcerec struct
   @param lj14Forces               output array of forces: lj14Forces[atomIndex][3]
   @param lj14AtomEnergies         output array of LJ 14 energies[atomIndex]
   @param totalEnergy              output total energy

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsLJ14ReferenceForce( const t_topology* top, RealOpenMM** atomCoordinates,
                               t_forcerec *fr, const RealOpenMM* charges, RealOpenMM** lj14Forces,
                               RealOpenMM* lj14AtomEnergies,
                               RealOpenMM* totalEnergy ){ 

   // ---------------------------------------------------------------------------------------

   static const char* methodName         = "\ngromacsLJ14ReferenceForce: ";
   static int debug                      = true;
   static int startAtom                  = 0;
   static int stopAtom                   = 3;

   static const RealOpenMM zero          = 0.0;
   static const RealOpenMM one           = 1.0;

   static const int  MaxFreeArrayIndex   = 20;

   int indexOneDArraysToFree             = 0;
   RealOpenMM* oneDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDArraysToFree             = 0;
   RealOpenMM** twoDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDIntArraysToFree          = 0;
   int** twoDIntArraysToFree[MaxFreeArrayIndex];

   // ---------------------------------------------------------------------------------------

   int numberOfAtoms = top->atoms.nr;

   // ---------------------------------------------------------------------------------------

   // get bond atom indices & parameters

   int** lj14AtomIndices;
   RealOpenMM** lj14Parameters;
   int numberOfBonds                                 = gromacsGetLJ14Parameters( top, &lj14AtomIndices, &lj14Parameters );
   twoDIntArraysToFree[indexTwoDIntArraysToFree++]   = lj14AtomIndices;
   twoDArraysToFree[indexTwoDArraysToFree++]         = lj14Parameters;

   // bond energy array

   RealOpenMM* lj14BondEnergies                      = (RealOpenMM*) malloc( sizeof( RealOpenMM )*numberOfBonds );
   memset( lj14BondEnergies, 0, sizeof( RealOpenMM )*numberOfBonds );
   oneDArraysToFree[indexOneDArraysToFree++]         = lj14BondEnergies;

   ReferenceLJ14 referenceLJ14;
   ReferenceBondForce referenceBondForce;

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << " lj 14 bonds=" << numberOfBonds << std::endl;
/*
      for( int ii  = 0; ii < numberOfBonds; ii++ ){
         message << "   bnd=" << ii << " a[";
         for( int jj  = 0; jj < 4; jj++ ){
            message << lj14AtomIndices[ii][jj] << " ";
         }
         message << "] p[";
         for( int jj  = 0; jj < 3; jj++ ){
            message << lj14Parameters[ii][jj] << " ";
         }
         message << "]" << std::endl; 
      }
*/
      SimTKOpenMMLog::printMessage( message );
   }

   // calculate derived parameters

   for( int ii  = 0; ii < numberOfBonds; ii++ ){
      referenceLJ14.getDerivedParameters( lj14Parameters[ii][0], lj14Parameters[ii][1],
                                          charges[lj14AtomIndices[ii][0]],
                                          charges[lj14AtomIndices[ii][1]], fr->epsfac*fr->fudgeQQ, lj14Parameters[ii] );
   }

   // lj14 bond forces

   referenceBondForce.calculateForce( numberOfBonds, lj14AtomIndices, atomCoordinates,
                                      lj14Parameters, lj14Forces, 
                                      lj14BondEnergies, lj14AtomEnergies, 
                                      totalEnergy, referenceLJ14);

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << "RbDihedral  Total energy=" << *totalEnergy << std::endl;
      for( int ii = startAtom; ii < stopAtom; ii++ ){
         message << ii << " F[";
         SimTKOpenMMUtilities::formatRealStringStream( message, lj14Forces[ii], 3, one );
         message << "] E=" << lj14AtomEnergies[ii] << std::endl;
      } 
      SimTKOpenMMLog::printMessage( message );
   }

   // ---------------------------------------------------------------------------------------

   // free memory

   gromacsFreeSundryArrays( indexOneDArraysToFree, indexTwoDArraysToFree, indexTwoDIntArraysToFree,
                            MaxFreeArrayIndex, oneDArraysToFree, twoDArraysToFree, twoDIntArraysToFree,
                            methodName );

   // ---------------------------------------------------------------------------------------

   return 0;
}

/**---------------------------------------------------------------------------------------

   Get bond atom indices and parameters for LJ Coulomb bonds

   @param top                      Gromacs t_topology data struct
   @param atomIndices              array of atomIndices[bondIndex][atomIndex]
   @param bondParameters           array of bond parameters[bondIndex][parameterIndex]

   The arrays atomIndices & bondParameters are allocated off the heap;
   The memory for the second index is one block of size (number of bonds)*(number of atoms in bond)

   The memory should be freed by the callee via calls:
		free( atomIndices[0] ); free( atomIndices );

   @return number of bonds

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsGetLJCoulombParameters( const t_topology* top, t_forcerec *fr, const t_mdatoms *md,
                                   int*** ljCoulombExclusions, RealOpenMM*** outputAtomParameters ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "\ngromacsGetLJCoulombParameters: ";
   static const int twoI           = 2;
   static const int threeI         = 3;
   static const int fourI          = 4;
   static const int sixI           = 6;

   static const RealOpenMM two     = 2.0;

   int numberOfAtoms, exclusionBlockSize, numberOfParameters;
   int** exclusions;
   int* exclusionBlock;

   RealOpenMM** atomParameters;
   RealOpenMM* atomParametersBlock;

   // ---------------------------------------------------------------------------------------

   const t_block *gmx_excl             = &top->atoms.excl;

   int ntypes                          = fr->ntype;
   int* types                          = md->typeA;
   real* nbfp                          = fr->nbfp;
   
   // ---------------------------------------------------------------------------------------

   numberOfAtoms                       = top->atoms.nr;
   numberOfParameters                  = 4;

   // allocate & initialize memory for exclusions and parameters

   exclusionBlockSize                  = gmx_excl->index[numberOfAtoms] - gmx_excl->index[0] + numberOfAtoms;
   exclusions                          = (int**) malloc( sizeof( int* )*numberOfAtoms );
   exclusionBlock                      = (int*)  malloc( sizeof( int  )*exclusionBlockSize );
   memset( exclusionBlock, 0, sizeof( int  )*exclusionBlockSize );
   *ljCoulombExclusions                = exclusions;

   atomParameters                      = (RealOpenMM**) malloc( sizeof( RealOpenMM* )*numberOfAtoms );
   atomParametersBlock                 = (RealOpenMM*)  malloc( sizeof( RealOpenMM )*numberOfAtoms*numberOfParameters );
   memset( atomParametersBlock, 0, sizeof( RealOpenMM )*numberOfAtoms*numberOfParameters );
   *outputAtomParameters               = atomParameters;

   for( int ii = 0; ii < numberOfAtoms; ii++ ){ 

      exclusions[ii]        = exclusionBlock;
      exclusionBlock       += (gmx_excl->index[ii+1] - gmx_excl->index[ii] + 1);

      atomParameters[ii]    = atomParametersBlock;
      atomParametersBlock  += numberOfParameters;

   }
   
   // ---------------------------------------------------------------------------------------

   // load parameters

   RealOpenMM* charges = md->chargeA;
   for( int ii = 0; ii < numberOfAtoms; ii++ ) { 

      exclusions[ii][0]        = gmx_excl->index[ii+1] - gmx_excl->index[ii];
      int offsetIndex          = 1;
      for ( int jj = gmx_excl->index[ii]; jj < gmx_excl->index[ii+1]; jj++ ) {
         exclusions[ii][offsetIndex++] = gmx_excl->a[jj];
      }    
 
      atomParameters[ii][0]     = nbfp[ types[ii]*2*ntypes + types[ii]*2 ];
      atomParameters[ii][1]     = nbfp[ types[ii]*2*ntypes + types[ii]*2 + 1 ]; 
      atomParameters[ii][2]     = charges[ii]; 

   }

   return 0;
}

/**---------------------------------------------------------------------------------------

   Calculate LJ Coulomb forces/energies for given configuration and topology

   @param top                      Gromacs t_topology data struct
   @param atomCoordinates          atom configuration
   @param fr                       Gromacs t_forcerec struct
   @param md                       Gromacs t_mdatoms struct
   @param ljCoulombForces          output array of forces: ljCoulombForces[atomIndex][3]
   @param ljCoulombAtomEnergies    output array of LJ Coulomb energies[atomIndex]
   @param totalEnergy              output total energy

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsLJCoulombReferenceForce( const t_topology* top, RealOpenMM** atomCoordinates,
                                    t_forcerec *fr, const t_mdatoms *md, RealOpenMM** ljCoulombForces,
                                    RealOpenMM* ljCoulombAtomEnergies,
                                    RealOpenMM* totalEnergy ){ 

   // ---------------------------------------------------------------------------------------

   static const char* methodName         = "\ngromacsLJCoulombReferenceForce: ";
   static int debug                      = true;
   static int startAtom                  = 0;
   static int stopAtom                   = 3;

   static const RealOpenMM zero          = 0.0;
   static const RealOpenMM one           = 1.0;

   static const int  MaxFreeArrayIndex   = 20;

   int indexOneDArraysToFree             = 0;
   RealOpenMM* oneDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDArraysToFree             = 0;
   RealOpenMM** twoDArraysToFree[MaxFreeArrayIndex];

   // currently not used!

   RealOpenMM fixedParameters[MaxFreeArrayIndex];

   int indexTwoDIntArraysToFree          = 0;
   int** twoDIntArraysToFree[MaxFreeArrayIndex];

   // ---------------------------------------------------------------------------------------

   int numberOfAtoms = top->atoms.nr;

   // ---------------------------------------------------------------------------------------

   // get exclusions & parameters

   int** ljCoulombExclusions;
   RealOpenMM** ljCoulombParameters;
   int numberOfBonds                                 = gromacsGetLJCoulombParameters( top, fr, md, &ljCoulombExclusions, &ljCoulombParameters );
   twoDIntArraysToFree[indexTwoDIntArraysToFree++]   = ljCoulombExclusions;
   twoDArraysToFree[indexTwoDArraysToFree++]         = ljCoulombParameters;

   ReferenceLJCoulombIxn referenceLJCoulomb;

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << " lj Coulomb atoms=" << numberOfAtoms << std::endl;
/*
      for( int ii  = 0; ii < numberOfAtoms; ii++ ){

         message << "   " << ii << " P[";
         for( int jj  = 0; jj < 3; jj++ ){
            message << ljCoulombParameters[ii][jj] << " ";
         }
         message << "]";

         message << "  X " << ljCoulombExclusions[ii][0] << " [";
         for( int jj  = 1; jj <= ljCoulombExclusions[ii][0]; jj++ ){
            message << ljCoulombExclusions[ii][jj] << " ";
         }
         message << "]" << std::endl; 
      }
*/
      SimTKOpenMMLog::printMessage( message );
   }

   // calculate derived parameters

   FILE* paramFile = fopen( "LJCoulombDerived.txt", "w" );
   RealOpenMM epsFacSqrt = SQRT( fr->epsfac );
   for( int ii  = 0; ii < numberOfAtoms; ii++ ){
      (void) fprintf( paramFile, "%d %12.5e %12.5e %12.5e", ii, ljCoulombParameters[ii][0], ljCoulombParameters[ii][1], ljCoulombParameters[ii][2] );
      referenceLJCoulomb.getDerivedParameters( ljCoulombParameters[ii][0], ljCoulombParameters[ii][1],
                                               ljCoulombParameters[ii][2],
                                               epsFacSqrt, ljCoulombParameters[ii] );
      (void) fprintf( paramFile, "   %12.5e %12.5e %12.5e\n", ljCoulombParameters[ii][0], ljCoulombParameters[ii][1], ljCoulombParameters[ii][2] );
   }
   (void) fclose( paramFile );

   // ljCoulomb bond forces

   referenceLJCoulomb.calculatePairIxn( numberOfAtoms, atomCoordinates,
                                        ljCoulombParameters, ljCoulombExclusions,
                                        fixedParameters, ljCoulombForces, 
                                        ljCoulombAtomEnergies, totalEnergy );

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << "LJCoulomb Total energy=" << *totalEnergy << std::endl;
      for( int ii = startAtom; ii < stopAtom; ii++ ){
         message << ii << " F[";
         SimTKOpenMMUtilities::formatRealStringStream( message, ljCoulombForces[ii], 3, one );
         message << "] E=" << ljCoulombAtomEnergies[ii] << std::endl;
      } 
      SimTKOpenMMLog::printMessage( message );
   }

   // ---------------------------------------------------------------------------------------

   // free memory

   gromacsFreeSundryArrays( indexOneDArraysToFree, indexTwoDArraysToFree, indexTwoDIntArraysToFree,
                            MaxFreeArrayIndex, oneDArraysToFree, twoDArraysToFree, twoDIntArraysToFree,
                            methodName );

   // ---------------------------------------------------------------------------------------

   return 0;
}

/**---------------------------------------------------------------------------------------

   Calculate forces for given configuration and topology

   @param top                      Gromacs t_topology data struct
   @param gromacAtomCoordinates    atom configuration
   @param md                       Gromacs t_mdatoms data struct
   @param partialChargesIn         array of partial charges
   @param log                      log reference (stdlog in md.c)

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsReferenceForce( const t_topology* top, const rvec* gromacAtomCoordinates,
                           const t_mdatoms *md, t_forcerec* fr,
                           char* baseFileName, FILE* log ){ 

   // ---------------------------------------------------------------------------------------

   static const char* methodName               = "\ngromacsReferenceForce: ";
   static int debug                            = true;
   static int startAtom                        = 0;
   static int stopAtom                         = 3;

   static const RealOpenMM zero                = 0.0;
   static const RealOpenMM one                 = 1.0;

   static const int  zeroI                     = 0;
   static const int  twoI                      = 2;
   static const int  threeI                    = 3;
   static const int  fourI                     = 4;

   static const int  MaxFreeArrayIndex         = 20;

   static const int  includeHarmonic           = 1;
   static const int  includeAngle              = 1;
   static const int  includeProperDihedral     = 1;
   static const int  includeRbDihedral         = 1;
   static const int  includeLJ14               = 1;
   static const int  includeLJCoulomb          = 1;

   static const int  outputParameterFile       = 0;
   static const int  outputForceFile           = 1;

   RealOpenMM totalEnergy                      = zero;

   int indexOneDArraysToFree                   = 0;
   RealOpenMM* oneDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDArraysToFree                   = 0;
   RealOpenMM** twoDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDIntArraysToFree                = 0;
   int** twoDIntArraysToFree[MaxFreeArrayIndex];

   // ---------------------------------------------------------------------------------------

   // set log file, if available

   if( log ){
      SimTKOpenMMLog::setSimTKOpenMMLog( log );
   }

   int numberOfAtoms = top->atoms.nr;

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << " atoms=" << numberOfAtoms;
      SimTKOpenMMLog::printMessage( message );
   }

   RealOpenMM** atomCoordinates              = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
   twoDArraysToFree[indexTwoDArraysToFree++] = atomCoordinates;

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      atomCoordinates[ii][0] = gromacAtomCoordinates[ii][0];
      atomCoordinates[ii][1] = gromacAtomCoordinates[ii][1];
      atomCoordinates[ii][2] = gromacAtomCoordinates[ii][2];
   }

   // accumulate bonded force contributions

   RealOpenMM** bondedForces                 = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
   twoDArraysToFree[indexTwoDArraysToFree++] = bondedForces;

   // ---------------------------------------------------------------------------------------

   // harmonic bond forces

   if( includeHarmonic ){
      RealOpenMM** harmonicBondForces                   = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
      twoDArraysToFree[indexTwoDArraysToFree++]         = harmonicBondForces;
   
      RealOpenMM* harmonicBondEnergies                  = (RealOpenMM*) malloc( sizeof( RealOpenMM )*numberOfAtoms );
      memset( harmonicBondEnergies, 0, sizeof( RealOpenMM )*numberOfAtoms );
      oneDArraysToFree[indexOneDArraysToFree++]         = harmonicBondEnergies;
   
      gromacsHarmonicBondReferenceForce( top, atomCoordinates, harmonicBondForces, harmonicBondEnergies, &totalEnergy ); 

      if( debug ){
         std::stringstream message;
         message << methodName;
         message << " Bond Total energy=" << totalEnergy << std::endl;
         for( int ii = startAtom; ii < stopAtom; ii++ ){
            message << ii << " F[";
            SimTKOpenMMUtilities::formatRealStringStream( message, harmonicBondForces[ii], 3, one );
            message << "] E=" << harmonicBondEnergies[ii] << std::endl;
         } 
         SimTKOpenMMLog::printMessage( message );
      }

      // output file
   
      if( outputForceFile ){
         std::stringstream bondForceFileName;
         bondForceFileName << baseFileName;
         bondForceFileName << "HarmonicBond.txt";
         ReferenceForce::writeForces( numberOfAtoms, twoI, atomCoordinates, harmonicBondForces,
                                      harmonicBondEnergies, bondForceFileName.str(), one, one, one );
      }

      // accumulate bonded forces

      if( bondedForces ){
         SimTKOpenMMUtilities::addTwoDimArray( numberOfAtoms, 3, harmonicBondForces, bondedForces );
      }
   }

   // ---------------------------------------------------------------------------------------

   // angle bond forces

   if( includeAngle ){

      RealOpenMM** angleBondForces                      = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
      twoDArraysToFree[indexTwoDArraysToFree++]         = angleBondForces;
   
      RealOpenMM* angleBondAtomEnergies                 = (RealOpenMM*) malloc( sizeof( RealOpenMM )*numberOfAtoms );
      memset( angleBondAtomEnergies, 0, sizeof( RealOpenMM )*numberOfAtoms );
      oneDArraysToFree[indexOneDArraysToFree++]         = angleBondAtomEnergies;
   
      gromacsAngleBondReferenceForce( top, atomCoordinates, angleBondForces, angleBondAtomEnergies, &totalEnergy ); 
   
      if( debug ){
         std::stringstream message;
         message << methodName;
         message << " Angle Total energy=" << totalEnergy << std::endl;
         for( int ii = startAtom; ii < stopAtom; ii++ ){
            message << ii << " F[";
            SimTKOpenMMUtilities::formatRealStringStream( message, angleBondForces[ii], 3, one );
            message << "] E=" << angleBondAtomEnergies[ii] << std::endl;
         } 
         SimTKOpenMMLog::printMessage( message );
      }

      // output file
   
      if( outputForceFile ){
         std::stringstream angleForceFileName;
         angleForceFileName << baseFileName;
         angleForceFileName << "AngleBond.txt";
         ReferenceForce::writeForces( numberOfAtoms, threeI, atomCoordinates, angleBondForces,
                                      angleBondAtomEnergies, angleForceFileName.str(), one, one, one );
      }

      // accumulate bonded forces

      if( bondedForces ){
         SimTKOpenMMUtilities::addTwoDimArray( numberOfAtoms, 3, angleBondForces, bondedForces );
      }

   }

   // ---------------------------------------------------------------------------------------

   // properDihedral bond forces

   if( includeProperDihedral ){

      RealOpenMM** properDihedralBondForces             = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
      twoDArraysToFree[indexTwoDArraysToFree++]         = properDihedralBondForces;
   
      RealOpenMM* properDihedralBondAtomEnergies        = (RealOpenMM*) malloc( sizeof( RealOpenMM )*numberOfAtoms );
      memset( properDihedralBondAtomEnergies, 0, sizeof( RealOpenMM )*numberOfAtoms );
      oneDArraysToFree[indexOneDArraysToFree++]         = properDihedralBondAtomEnergies;
   
      gromacsProperDihedralBondReferenceForce( top, atomCoordinates, properDihedralBondForces, properDihedralBondAtomEnergies, &totalEnergy ); 
   
      if( debug ){
         std::stringstream message;
         message << methodName;
         message << " ProperDihedral Total energy=" << totalEnergy << std::endl;
         for( int ii = startAtom; ii < stopAtom; ii++ ){
            message << ii << " F[";
            SimTKOpenMMUtilities::formatRealStringStream( message, properDihedralBondForces[ii], 3, one );
            message << "] E=" << properDihedralBondAtomEnergies[ii] << std::endl;
         } 
         SimTKOpenMMLog::printMessage( message );
      }

      // output file
   
      if( outputForceFile ){
         std::stringstream properDihedralForceFileName;
         properDihedralForceFileName << baseFileName;
         properDihedralForceFileName << "ProperDihedral.txt";
         ReferenceForce::writeForces( numberOfAtoms, fourI, atomCoordinates, properDihedralBondForces,
                                      properDihedralBondAtomEnergies, properDihedralForceFileName.str(), one, one, one );
      }

      // accumulate bonded forces

      if( bondedForces ){
         SimTKOpenMMUtilities::addTwoDimArray( numberOfAtoms, 3, properDihedralBondForces, bondedForces );
      }
   }

   // ---------------------------------------------------------------------------------------

   // rbDihedral bond forces

   if( includeRbDihedral ){

      RealOpenMM** rbDihedralBondForces             = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
      twoDArraysToFree[indexTwoDArraysToFree++]     = rbDihedralBondForces;
   
      RealOpenMM* rbDihedralBondAtomEnergies        = (RealOpenMM*) malloc( sizeof( RealOpenMM )*numberOfAtoms );
      memset( rbDihedralBondAtomEnergies, 0, sizeof( RealOpenMM )*numberOfAtoms );
      oneDArraysToFree[indexOneDArraysToFree++]     = rbDihedralBondAtomEnergies;
   
      gromacsRbDihedralBondReferenceForce( top, atomCoordinates, rbDihedralBondForces, rbDihedralBondAtomEnergies, &totalEnergy ); 
  
      if( debug ){
         std::stringstream message;
         message << methodName;
         message << " RbDihedral Total energy=" << totalEnergy << std::endl;
         for( int ii = startAtom; ii < stopAtom; ii++ ){
            message << ii << " F[";
            SimTKOpenMMUtilities::formatRealStringStream( message, rbDihedralBondForces[ii], 3, one );
            message << "] E=" << rbDihedralBondAtomEnergies[ii] << std::endl;
         } 
         SimTKOpenMMLog::printMessage( message );
      }

      // output file
   
      if( outputForceFile ){
         std::stringstream rbDihedralForceFileName;
         rbDihedralForceFileName << baseFileName;
         rbDihedralForceFileName << "RbDihedral.txt";
         ReferenceForce::writeForces( numberOfAtoms, fourI, atomCoordinates, rbDihedralBondForces,
                                      rbDihedralBondAtomEnergies, rbDihedralForceFileName.str(), one, one, one );
      }

      // accumulate bonded forces

      if( bondedForces ){
         SimTKOpenMMUtilities::addTwoDimArray( numberOfAtoms, 3, rbDihedralBondForces, bondedForces );
      }
   }

   // ---------------------------------------------------------------------------------------

   // lj14 bond forces

   if( includeLJ14 ){

      RealOpenMM** lj14BondForces                   = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
      twoDArraysToFree[indexTwoDArraysToFree++]     = lj14BondForces;
   
      RealOpenMM* lj14BondAtomEnergies              = (RealOpenMM*) malloc( sizeof( RealOpenMM )*numberOfAtoms );
      memset( lj14BondAtomEnergies, 0, sizeof( RealOpenMM )*numberOfAtoms );
      oneDArraysToFree[indexOneDArraysToFree++]     = lj14BondAtomEnergies;
   
      gromacsLJ14ReferenceForce( top, atomCoordinates, fr, md->chargeA, lj14BondForces, lj14BondAtomEnergies, &totalEnergy ); 
  
      if( debug ){
         std::stringstream message;
         message << methodName;
         message << " LJ14 Total energy=" << totalEnergy << std::endl;
         for( int ii = startAtom; ii < stopAtom; ii++ ){
            message << ii << " F[";
            SimTKOpenMMUtilities::formatRealStringStream( message, lj14BondForces[ii], 3, one );
            message << "] E=" << lj14BondAtomEnergies[ii] << std::endl;
         } 
         SimTKOpenMMLog::printMessage( message );
      }

      // output file
   
      if( outputForceFile ){
         std::stringstream lj14ForceFileName;
         lj14ForceFileName << baseFileName;
         lj14ForceFileName << "LJ14.txt";
         ReferenceForce::writeForces( numberOfAtoms, fourI, atomCoordinates, lj14BondForces,
                                      lj14BondAtomEnergies, lj14ForceFileName.str(), one, one, one );
      }

      // accumulate bonded forces

      if( bondedForces ){
         SimTKOpenMMUtilities::addTwoDimArray( numberOfAtoms, 3, lj14BondForces, bondedForces );
      }
   }

   // output file

   if( outputForceFile && bondedForces ){
      std::stringstream bondedForceFileName;
      bondedForceFileName << baseFileName;
      bondedForceFileName << "Bonded.txt";
      ReferenceForce::writeForces( numberOfAtoms, zeroI, atomCoordinates, bondedForces,
                                   NULL, bondedForceFileName.str(), one, one, one );
   }

// ---------------------------------------------------------------------------------------

   // ljCoulomb bond forces

   if( includeLJCoulomb ){

      RealOpenMM** ljCoulombBondForces              = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
      twoDArraysToFree[indexTwoDArraysToFree++]     = ljCoulombBondForces;
   
      RealOpenMM* ljCoulombBondAtomEnergies         = (RealOpenMM*) malloc( sizeof( RealOpenMM )*numberOfAtoms );
      memset( ljCoulombBondAtomEnergies, 0, sizeof( RealOpenMM )*numberOfAtoms );
      oneDArraysToFree[indexOneDArraysToFree++]     = ljCoulombBondAtomEnergies;
   
      gromacsLJCoulombReferenceForce( top, atomCoordinates, fr, md, ljCoulombBondForces, ljCoulombBondAtomEnergies, &totalEnergy ); 
  
      if( debug ){
         std::stringstream message;
         message << methodName;
         message << " LJCoulomb Total energy=" << totalEnergy << std::endl;
         for( int ii = startAtom; ii < stopAtom; ii++ ){
            message << ii << " F[";
            SimTKOpenMMUtilities::formatRealStringStream( message, ljCoulombBondForces[ii], 3, one );
            message << "] E=" << ljCoulombBondAtomEnergies[ii] << std::endl;
         } 
         SimTKOpenMMLog::printMessage( message );
      }

      // output file
   
      if( outputForceFile ){
         std::stringstream ljCoulombForceFileName;
         ljCoulombForceFileName << baseFileName;
         ljCoulombForceFileName << "LJCoulomb.txt";
         ReferenceForce::writeForces( numberOfAtoms, fourI, atomCoordinates, ljCoulombBondForces,
                                      ljCoulombBondAtomEnergies, ljCoulombForceFileName.str(), one, one, one );
      }
   }

   // ---------------------------------------------------------------------------------------

   // free memory

   gromacsFreeSundryArrays( indexOneDArraysToFree, indexTwoDArraysToFree, indexTwoDIntArraysToFree,
                            MaxFreeArrayIndex, oneDArraysToFree, twoDArraysToFree, twoDIntArraysToFree,
                            methodName );

   return 0;
}

/**---------------------------------------------------------------------------------------

   Get Shake parameters and atom indices (simple version)

   @param top                      Gromacs t_topology data struct
   @param md                       Gromacs t_mdatoms data struct
   @param outputAtomIndices        array of atomIndices -- see below
   @param outputShakeParameters    array of Shake parameters -- see below

   The arrays atomIndices & outputShakeParameters are allocated off the heap;
   The memory for the second index is one block of size (number of blocks)*(second dimension)

   The memory should be freed by the callee via calls:
		free( outputAtomIndices[0] ); free( outputAtomIndices );
		free( outputShakeParameters[0] ); free( outputShakeParameters );


   outputAtomIndices has the following structure:

   outputAtomIndices[blockIndex][i], where

			i = 0:   index of atom1
			i = 1:   index of atom2

   outputShakeParameters has the following structure:

   outputShakeParameters[blockIndex][i], where

			i = 0:  constraint distance

   @return number of constraints

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsGetShakeParameters( const t_topology* top, const t_mdatoms *md,
                               int*** outputAtomIndices, RealOpenMM*** outputShakeParameters ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "\ngromacsGetShakeParameters: ";
   static const int twoI           = 2;
   static const int threeI         = 3;
   static const int fourI          = 4;
   static const int fiveI          = 5;
   static const int sixI           = 6;

   static const int debug          = 1;

   static const RealOpenMM zero    = 0.0;
   static const RealOpenMM one     = 1.0;
   static const RealOpenMM two     = 2.0;
   static const RealOpenMM epsilon = (RealOpenMM) 1.0e-05;

   // ---------------------------------------------------------------------------------------

   const t_idef* idef           = &top->idef;
   const t_iatom* atoms         = idef->il[ F_SHAKE ].iatoms;
   int numberOfConstraints      = idef->il[ F_SHAKE ].nr/3;

   gromacsAllocateBondIxnArrays( numberOfConstraints, twoI, outputAtomIndices, twoI, outputShakeParameters );
   int** atomIndices            = *outputAtomIndices;
   RealOpenMM** shakeParameters = *outputShakeParameters;

   // loop over constraints 

   int offset     = 0;
   for( int ii = 0; ii < numberOfConstraints; ii++ ) { 

      int type                    = atoms[ offset     ];  
      int atomI                   = atoms[ offset + 1 ];  
      int atomJ                   = atoms[ offset + 2 ];  

      offset                     += 3;

      atomIndices[ii][0]          = atomI;
      atomIndices[ii][1]          = atomJ;

      shakeParameters[ii][0]      = idef->iparams[type].shake.dA;
      shakeParameters[ii][1]      = zero; // unused
   }

   return numberOfConstraints;
}

/**---------------------------------------------------------------------------------------

   Get Shake parameters and atom indices in blocks

   @param top                      Gromacs t_topology data struct
   @param md                       Gromacs t_mdatoms data struct
   @param outputAtomIndices        array of atomIndices -- see below
   @param outputShakeParameters    array of Shake parameters -- see below

   The arrays atomIndices & outputShakeParameters are allocated off the heap;
   The memory for the second index is one block of size (number of blocks)*(second dimension [2 or 5] )

   The memory should be freed by the callee via calls:
		free( outputAtomIndices[0] ); free( outputAtomIndices );
		free( outputShakeParameters[0] ); free( outputShakeParameters );



   outputAtomIndices has the following structure:

   outputAtomIndices[blockIndex][i], where

			i = 0:   number of atoms in block (max is 4)
			i = 1:   index of heavy atom in constraint block
			i = 2-4: index of light atoms in constraint block

   If max > 4, then program aborts



   outputShakeParameters has the following structure:

   outputShakeParameters[blockIndex][i], where

			i = 0:  inverse mass of heavy atom
			i = 2:  1.0/2*[ inverse mass of heavy atom + inverse mass of light atom ]
			i = 3:  distance constraint**2

   The following assumptions are made:

			all light atoms have same mass
			constraint distance is the same for all constraints in a block

   If these are violated, the program should abort


   @return number of constraint blocks

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsGetBlockShakeParameters( const t_topology* top, const t_mdatoms *md,
                                    int*** outputAtomIndices, RealOpenMM*** outputShakeParameters ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "\ngromacsGetBlockShakeParameters: ";
   static const int twoI           = 2;
   static const int threeI         = 3;
   static const int fourI          = 4;
   static const int fiveI          = 5;
   static const int sixI           = 6;

   static const int debug          = 1;

   static const RealOpenMM zero    = 0.0;
   static const RealOpenMM one     = 1.0;
   static const RealOpenMM two     = 2.0;
   static const RealOpenMM epsilon = (RealOpenMM) 1.0e-05;

   // ---------------------------------------------------------------------------------------

   const t_idef* idef          = &top->idef;
   const t_iatom* atoms        = idef->il[ F_SHAKE ].iatoms;
   int numberOfConstraints     = idef->il[ F_SHAKE ].nr/3;

   // loop over constraints 
   // create a ReferenceShakeConstraint object w/ info for that constraint
   // add the constraint to a map of constraints based on the index of the heavy atom
   //      shakeMap[heavyAtomIndex] = Vector{ ReferenceShakeConstraints } containing that atom

   int offset = 0;
   IntShakeMap shakeMap;
   for( int ii = 0; ii < numberOfConstraints; ii++ ) { 

      int type       = atoms[ offset     ];  
      int atomI      = atoms[ offset + 1 ];  
      int atomJ      = atoms[ offset + 2 ];  
      offset        += 3;
      ReferenceShakeConstraint* referenceShakeConstraint = new ReferenceShakeConstraint( atomI, atomJ, idef->iparams[type].shake.dA,
                                                                                         md->invmass[atomI], md->invmass[atomJ] );
      int heavyAtomIndex = referenceShakeConstraint->getHeavyAtomIndex();

      IntShakeMapI shakeMapI= shakeMap.find( heavyAtomIndex );
      ShakeVector* shakeVector;
      if( shakeMapI == shakeMap.end() ){
         shakeVector              = new ShakeVector();
         shakeMap[heavyAtomIndex] = shakeVector;
      } else {
         shakeVector = (*shakeMapI).second;
      }
      shakeVector->push_back( referenceShakeConstraint );
   }

   // validate blocks
   //    (1) at most 3 elements/block
   //    (2) distances the same
   //    (3) masses the same for light atoms

   int numberOfConstraintBlocks = (int) shakeMap.size();
   int blockCount               = 1;
   int errors                   = 0;
   std::stringstream message;
   message << "Number of Shake constraint blocks=" << numberOfConstraintBlocks;
   for( IntShakeMapI ii = shakeMap.begin(); ii != shakeMap.end(); ii++ ){ 
      ShakeVector* shakeVector = (*ii).second;
      message << "\nBlock " << blockCount++;
      if( shakeVector->size() > 3 ){
         message << "\n      ERROR: Block has too many elements(=" << shakeVector->size() << std::endl;
         errors++;
      }
      RealOpenMM blockDistance     = -one;
      RealOpenMM blockInverseMass  = -one;
      for( ShakeVectorI jj = shakeVector->begin(); jj != shakeVector->end(); jj++ ) { 

         (*jj)->printState( message );
         message << "   ";

         // distance check

         if( blockDistance < zero ){
            blockDistance = (*jj)->getConstraintDistance();
         } else {
            RealOpenMM constraintDistance = (*jj)->getConstraintDistance();
            RealOpenMM diff               = FABS( blockDistance - constraintDistance );
            if( diff > epsilon ){
               errors++;
               message << "\n      ERROR: distances not equal for block[" << blockDistance << " " << constraintDistance << "] diff=" << diff << std::endl;
            }
         }

         // mass check

         if( blockInverseMass < zero ){
            blockInverseMass = (*jj)->getLightAtomInverseMass();
         } else {
            RealOpenMM lightMass          = (*jj)->getLightAtomInverseMass();
            RealOpenMM diff               = FABS( blockInverseMass - lightMass );
            if( diff > epsilon ){
               errors++;
               message << "\n      ERROR: light atom masses not equal for block[" << blockInverseMass << " " << lightMass << "] diff=" << diff << std::endl;
            }
         }
      }
   }

   // abort if errors

   if( errors ){
      message << "\n\nErrors=" << errors << " -- aborting.";
      SimTKOpenMMLog::printError( message );
   } else {
      SimTKOpenMMLog::printMessage( message );
   }

   // allocate & initialize memory for Shake atom indices & parameters
   // and load data into arrays

   gromacsAllocateBondIxnArrays( numberOfConstraintBlocks, fiveI, outputAtomIndices, threeI, outputShakeParameters );
   int** atomIndices            = *outputAtomIndices;
   RealOpenMM** shakeParameters = *outputShakeParameters;

   blockCount = 0;
   for( IntShakeMapI ii = shakeMap.begin(); ii != shakeMap.end(); ii++, blockCount++ ){ 
      ShakeVector* shakeVector   = (*ii).second;
      atomIndices[blockCount][0] = (int) (shakeVector->size() + 1);
      int lightAtomCount         = 2;
      for( ShakeVectorI jj = shakeVector->begin(); jj != shakeVector->end(); jj++ ) { 

         if( lightAtomCount == 2 ){
            atomIndices[blockCount][1]                = (*jj)->getHeavyAtomIndex();

            shakeParameters[blockCount][0]            = (*jj)->getHeavyAtomInverseMass();
   
            shakeParameters[blockCount][1]            = two*( (*jj)->getHeavyAtomInverseMass() +  (*jj)->getLightAtomInverseMass() );
            shakeParameters[blockCount][1]            = one/shakeParameters[blockCount][1];
   
            shakeParameters[blockCount][2]            = (*jj)->getConstraintDistance()* (*jj)->getConstraintDistance();
         }
         atomIndices[blockCount][lightAtomCount++]    = (*jj)->getLightAtomIndex();
      }
   }
         
   // free allocated memory

   for( IntShakeMapI ii = shakeMap.begin(); ii != shakeMap.end(); ii++ ) { 
      ShakeVector* shakeVector = (*ii).second;
      for( ShakeVectorI jj = shakeVector->begin(); jj != shakeVector->end(); jj++ ) { 
         delete *jj;
      }
     delete shakeVector;
   }

   return numberOfConstraintBlocks;
}

/**---------------------------------------------------------------------------------------

   Update coordinates using stochastic dynamics integrator

   @param top                      Gromacs t_topology data struct
   @param gromacAtomCoordinates    atom configuration
   @param md                       Gromacs t_mdatoms data struct
   @param fr                       Gromacs t_forcerec data struct
   @param gromacsVelocities        velocities
   @param gromacsForces            forces
   @param masses                   masses
   @param deltaT                   delta t
   @param tau                      viscosity(?)
   @param temperature              temperature
   @param baseFileName             base file name
   @param log                      log reference (stdlog in md.c)

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsReferenceSdUpdate( int numberOfAtoms, const t_topology* top, const rvec* gromacAtomCoordinates,
                              const t_mdatoms *md, t_forcerec* fr,
                              const rvec* gromacsVelocities, const rvec* gromacsForces,
                              real* masses, real deltaT, real tau, real temperature,
                              char* baseFileName, FILE* log ){ 

   // ---------------------------------------------------------------------------------------

   static const char* methodName               = "\ngromacsReferenceSdUpdate: ";
   static int debug                            = true;

   static const RealOpenMM zero                = 0.0;
   static const RealOpenMM one                 = 1.0;

   static const int  zeroI                     = 0;
   static const int  twoI                      = 2;
   static const int  threeI                    = 3;
   static const int  fourI                     = 4;

   static const int  MaxFreeArrayIndex         = 20;

   static const int  outputParameterFile       = 0;
   static const int  outputForceFile           = 1;

   RealOpenMM totalEnergy                      = zero;

   int indexOneDArraysToFree                   = 0;
   RealOpenMM* oneDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDArraysToFree                   = 0;
   RealOpenMM** twoDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDIntArraysToFree                = 0;
   int** twoDIntArraysToFree[MaxFreeArrayIndex];

   // ---------------------------------------------------------------------------------------

   // set log file, if available

   if( log ){
      SimTKOpenMMLog::setSimTKOpenMMLog( log );
   }

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << " atoms=" << numberOfAtoms;
      SimTKOpenMMLog::printMessage( message );
   }

   // ---------------------------------------------------------------------------------------

   // load coordinates, velocities, and forces

   RealOpenMM** atomCoordinates              = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
   twoDArraysToFree[indexTwoDArraysToFree++] = atomCoordinates;

   RealOpenMM** velocities                   = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
   twoDArraysToFree[indexTwoDArraysToFree++] = velocities;

   RealOpenMM** forces                       = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
   twoDArraysToFree[indexTwoDArraysToFree++] = forces;

   for( int ii = 0; ii < numberOfAtoms; ii++ ){

      atomCoordinates[ii][0] = gromacAtomCoordinates[ii][0];
      atomCoordinates[ii][1] = gromacAtomCoordinates[ii][1];
      atomCoordinates[ii][2] = gromacAtomCoordinates[ii][2];

      velocities[ii][0]      = gromacsVelocities[ii][0];
      velocities[ii][1]      = gromacsVelocities[ii][1];
      velocities[ii][2]      = gromacsVelocities[ii][2];

      forces[ii][0]          = gromacsForces[ii][0];
      forces[ii][1]          = gromacsForces[ii][1];
      forces[ii][2]          = gromacsForces[ii][2];

   }

   // ---------------------------------------------------------------------------------------

   // perform update

   // create referenceStochasticDynamics object
   // get Shake parameters

   ReferenceStochasticDynamics referenceStochasticDynamics( numberOfAtoms, deltaT, tau, temperature ); 

   if( debug ){
      std::stringstream message;
      message << methodName << " Parameters:\n   ";
      referenceStochasticDynamics.printParameters( message );
      message << std::endl;
      SimTKOpenMMLog::printMessage( message );
   }

   int** shakeAtomIndices;
   RealOpenMM** shakeParameters;
   //int numberOfConstraintBlocks                      = gromacsGetBlockShakeParameters( top, md, &shakeAtomIndices, &shakeParameters );
   int numberOfConstraints                           = gromacsGetShakeParameters( top, md, &shakeAtomIndices, &shakeParameters );
   twoDIntArraysToFree[indexTwoDIntArraysToFree++]   = shakeAtomIndices;
   twoDArraysToFree[indexTwoDArraysToFree++]         = shakeParameters;

   if( debug ){
/*
      // block constraints

      std::stringstream message;
      message << methodName << " Shake Parameters: " << numberOfConstraintBlocks << "\n";
      for( int ii = 0; ii < numberOfConstraintBlocks; ii++ ){
         message << " " << (ii+1) << " " << shakeParameters[ii][0] << " " << shakeParameters[ii][1] << " " << shakeParameters[ii][2] << " ";
         message << shakeAtomIndices[ii][0] << " [";
         for( int jj = 0; jj < shakeAtomIndices[ii][0]; jj++ ){
            message << shakeAtomIndices[ii][jj+1] << " ";
         }
         message << "]\n";
      }
      message << std::endl;
      SimTKOpenMMLog::printMessage( message );
*/

      
      std::stringstream message;
      message << methodName << " Shake Parameters: " << numberOfConstraints << "\n";
      for( int ii = 0; ii < numberOfConstraints && ii < 10; ii++ ){
         message << " " << (ii+1) << " " << shakeParameters[ii][0] << " [" << shakeAtomIndices[ii][0] << " " << shakeAtomIndices[ii][1] << "]\n";
      }
      message << std::endl;
      SimTKOpenMMLog::printMessage( message );
   }

   // add Shake to dynamics

   ReferenceShakeAlgorithm referenceShakeAlgorithm( numberOfConstraints, shakeAtomIndices, shakeParameters );
   referenceStochasticDynamics.setReferenceShakeAlgorithm( &referenceShakeAlgorithm );

	referenceStochasticDynamics.writeState( numberOfAtoms, atomCoordinates,
                                           velocities, forces, masses,
                                           zeroI, baseFileName );

   referenceStochasticDynamics.update( numberOfAtoms, atomCoordinates, velocities, forces, masses );

   if( debug ){
      std::stringstream message;
      message << methodName << " Called update: " << referenceStochasticDynamics.getTimeStep() << "\n";
      SimTKOpenMMLog::printMessage( message );
   }

	referenceStochasticDynamics.writeState( numberOfAtoms, atomCoordinates,
                                           velocities, forces, masses,
                                           zeroI, baseFileName );

   // ---------------------------------------------------------------------------------------

   // free memory

   gromacsFreeSundryArrays( indexOneDArraysToFree, indexTwoDArraysToFree, indexTwoDIntArraysToFree,
                            MaxFreeArrayIndex, oneDArraysToFree, twoDArraysToFree, twoDIntArraysToFree,
                            methodName );

   return 0;
}

/**---------------------------------------------------------------------------------------

   Update coordinates using stochastic dynamics integrator

   @param top                      Gromacs t_topology data struct
   @param gromacAtomCoordinates    atom configuration
   @param md                       Gromacs t_mdatoms data struct
   @param fr                       Gromacs t_forcerec data struct
   @param gromacsVelocities        velocities
   @param gromacsForces            forces
   @param masses                   masses
   @param timesteps                number of timesteps
   @param deltaT                   delta t
   @param tau                      viscosity(?)
   @param temperature              temperature
   @param baseFileName             base file name
   @param log                      log reference (stdlog in md.c)

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsReferenceSimulation( int numberOfAtoms, const t_topology* top, const rvec* gromacAtomCoordinates,
                                const t_mdatoms *md, t_forcerec* fr,
                                const rvec* gromacsVelocities, const rvec* gromacsForces,
                                real* masses, int timesteps, real deltaT, real tau, real temperature,
                                char* baseFileName, FILE* log ){ 

   // ---------------------------------------------------------------------------------------

   static const char* methodName               = "\ngromacsReferenceSimulation: ";
   static int debug                            = false;

   static const RealOpenMM zero                = 0.0;
   static const RealOpenMM one                 = 1.0;

   static const int  zeroI                     = 0;
   static const int  twoI                      = 2;
   static const int  threeI                    = 3;
   static const int  fourI                     = 4;

   static const int  MaxFreeArrayIndex         = 100;

   static const int  outputParameterFile       = 0;
   static const int  outputForceFile           = 1;

   RealOpenMM fixedParameters[2]; // unused!

   RealOpenMM totalEnergy                      = zero;

   int indexOneDArraysToFree                   = 0;
   RealOpenMM* oneDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDArraysToFree                   = 0;
   RealOpenMM** twoDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDIntArraysToFree                = 0;
   int** twoDIntArraysToFree[MaxFreeArrayIndex];

   // ---------------------------------------------------------------------------------------

   // set log file, if available

   if( log ){
      SimTKOpenMMLog::setSimTKOpenMMLog( log );
   }

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << " atoms=" << numberOfAtoms;
      SimTKOpenMMLog::printMessage( message );
   }

   // ---------------------------------------------------------------------------------------

   // load coordinates, velocities, and forces

   RealOpenMM** atomCoordinates              = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
   twoDArraysToFree[indexTwoDArraysToFree++] = atomCoordinates;

   RealOpenMM** velocities                   = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
   twoDArraysToFree[indexTwoDArraysToFree++] = velocities;

   RealOpenMM** forces                       = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
   twoDArraysToFree[indexTwoDArraysToFree++] = forces;

   RealOpenMM* energies                      = SimTKOpenMMUtilities::allocateOneDRealOpenMMArray( numberOfAtoms, NULL, true, (RealOpenMM) 0.0 );
   oneDArraysToFree[indexOneDArraysToFree++] = energies;

   for( int ii = 0; ii < numberOfAtoms; ii++ ){

      atomCoordinates[ii][0] = gromacAtomCoordinates[ii][0];
      atomCoordinates[ii][1] = gromacAtomCoordinates[ii][1];
      atomCoordinates[ii][2] = gromacAtomCoordinates[ii][2];

      velocities[ii][0]      = gromacsVelocities[ii][0];
      velocities[ii][1]      = gromacsVelocities[ii][1];
      velocities[ii][2]      = gromacsVelocities[ii][2];
   }

   // ---------------------------------------------------------------------------------------

   // force setup

   enum ForceTypes { HarmonicBond, AngleBond, ProperDihedral, RbDihedral, LJ14, LJCoulomb, MaxForceTypes };
 
   int** bondAtomIndices[MaxForceTypes];
   RealOpenMM** bondParameters[MaxForceTypes];
   int numberOfBonds[MaxForceTypes];
   ReferenceBondIxn* referenceBondIxn[MaxForceTypes];

   ReferenceBondForce referenceBondForce;

   for( int ii = 0; ii < MaxForceTypes; ii++ ){
      bondAtomIndices[ii]     = NULL;
      bondParameters[ii]      = NULL;
      numberOfBonds[ii]       = 0;
      referenceBondIxn[ii]    = NULL;
   }

   // ---------------------------------------------------------------------------------------

   // harmonic bonds

   numberOfBonds[HarmonicBond]                       = gromacsGetHarmonicBondParameters( top, &bondAtomIndices[HarmonicBond], &bondParameters[HarmonicBond] );
   referenceBondIxn[HarmonicBond]                    = new ReferenceHarmonicBondIxn();

   // angle bonds

   numberOfBonds[AngleBond]                          = gromacsGetAngleBondParameters( top, &bondAtomIndices[AngleBond], &bondParameters[AngleBond] );
   referenceBondIxn[AngleBond]                       = new ReferenceAngleBondIxn();
   
   // proper dihedral

   numberOfBonds[ProperDihedral]                     = gromacsGetProperDihedralBondParameters( top, &bondAtomIndices[ProperDihedral], &bondParameters[ProperDihedral] );
   referenceBondIxn[ProperDihedral]                  = new ReferenceProperDihedralBond();

   // rb dihedral

   numberOfBonds[RbDihedral]                         = gromacsGetRbDihedralBondParameters( top, &bondAtomIndices[RbDihedral], &bondParameters[RbDihedral] );
   referenceBondIxn[RbDihedral]                      = new ReferenceRbDihedralBond();

   // 1-4 ixns

   numberOfBonds[LJ14]                               = gromacsGetLJ14Parameters( top, &bondAtomIndices[LJ14], &bondParameters[LJ14] );
   ReferenceLJ14* referenceLJ14                      = new ReferenceLJ14();
   referenceBondIxn[LJ14]                            = referenceLJ14;

   // calculate derived parameters for 1-4 LJ

   RealOpenMM* charges = md->chargeA;
   for( int ii  = 0; ii < numberOfBonds[LJ14]; ii++ ){
      referenceLJ14->getDerivedParameters( bondParameters[LJ14][ii][0],
                                           bondParameters[LJ14][ii][1],
                                           charges[bondAtomIndices[LJ14][ii][0]],
                                           charges[bondAtomIndices[LJ14][ii][1]], fr->epsfac*fr->fudgeQQ, bondParameters[LJ14][ii] );
   }

   // LJ Coulomb

   gromacsGetLJCoulombParameters( top, fr, md, &bondAtomIndices[LJCoulomb], &bondParameters[LJCoulomb] );
   ReferenceLJCoulombIxn referenceLJCoulomb;

   // calculate derived parameters

   RealOpenMM epsFacSqrt = SQRT( fr->epsfac );
   for( int ii  = 0; ii < numberOfAtoms; ii++ ){
      referenceLJCoulomb.getDerivedParameters( bondParameters[LJCoulomb][ii][0], bondParameters[LJCoulomb][ii][1],
                                               bondParameters[LJCoulomb][ii][2],
                                               epsFacSqrt, bondParameters[LJCoulomb][ii] );
   }

   // ---------------------------------------------------------------------------------------

   // record arrays to be freed

   for( int ii = 0; ii < MaxForceTypes; ii++ ){
      if( bondAtomIndices[ii] != NULL ){
         twoDIntArraysToFree[indexTwoDIntArraysToFree++] = bondAtomIndices[ii];
      }
      if( bondParameters[ii] != NULL ){
         twoDArraysToFree[indexTwoDArraysToFree++]       = bondParameters[ii];
      }
   }

   // ---------------------------------------------------------------------------------------

   // update() setup

   // create referenceStochasticDynamics object
   // get Shake parameters

   ReferenceStochasticDynamics referenceStochasticDynamics( numberOfAtoms, deltaT, tau, temperature ); 

   if( debug ){
      std::stringstream message;
      message << methodName << " Parameters:\n   ";
      referenceStochasticDynamics.printParameters( message );
      message << std::endl;
      SimTKOpenMMLog::printMessage( message );
   }

   int** shakeAtomIndices;
   RealOpenMM** shakeParameters;
   int numberOfConstraints                           = gromacsGetShakeParameters( top, md, &shakeAtomIndices, &shakeParameters );
   twoDIntArraysToFree[indexTwoDIntArraysToFree++]   = shakeAtomIndices;
   twoDArraysToFree[indexTwoDArraysToFree++]         = shakeParameters;

   if( debug ){
      std::stringstream message;
      message << methodName << " Shake Parameters: " << numberOfConstraints << "\n";
      for( int ii = 0; ii < numberOfConstraints && ii < 10; ii++ ){
         message << " " << (ii+1) << " " << shakeParameters[ii][0] << " [" << shakeAtomIndices[ii][0] << " " << shakeAtomIndices[ii][1] << "]\n";
      }
      message << std::endl;
      SimTKOpenMMLog::printMessage( message );
   }

   // add Shake to dynamics

   ReferenceShakeAlgorithm referenceShakeAlgorithm( numberOfConstraints, shakeAtomIndices, shakeParameters );
   referenceStochasticDynamics.setReferenceShakeAlgorithm( &referenceShakeAlgorithm );

	referenceStochasticDynamics.writeState( numberOfAtoms, atomCoordinates,
                                           velocities, forces, masses,
                                           zeroI, baseFileName );

   // ---------------------------------------------------------------------------------------

   // time steps

   for( int step = 0; step < timesteps; step++ ){

      // zero forces

      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         forces[ii][0]   = zero;
         forces[ii][1]   = zero;
         forces[ii][2]   = zero;
         energies[ii]    = zero;
      }

      // calculate bond forces

      for( int ii = 0; ii < MaxForceTypes; ii++ ){
         if( numberOfBonds[ii] > 0 ){
            referenceBondForce.calculateForce( numberOfBonds[ii], bondAtomIndices[ii], atomCoordinates,
                                               bondParameters[ii], forces,
                                               NULL, energies, &totalEnergy, *(referenceBondIxn[ii]) );
         }
      }

      if( debug ){
         std::stringstream message;
         int maxAtom = 5;
         message << methodName << " Bonded forces:\n";
         message << " Bond counts: [";
         for( int ii = 0; ii < MaxForceTypes; ii++ ){
            message <<  numberOfBonds[ii] << " ";
         }
         message << "]";
         message << std::endl;

         for( int ii = 0; ii < maxAtom; ii++ ){
            message << ii;
            message << " x[";
            SimTKOpenMMUtilities::formatRealStringStream( message, atomCoordinates[ii], 3, one );
            message << "] f[";
            SimTKOpenMMUtilities::formatRealStringStream( message, forces[ii], 3, one );
            message << "]";
            message << std::endl;
         }
         SimTKOpenMMLog::printMessage( message );
      }
   
      // LJ Coulomb forces
      // bondAtomIndices contains exclusions

      referenceLJCoulomb.calculatePairIxn( numberOfAtoms, atomCoordinates,
                                           bondParameters[LJCoulomb], bondAtomIndices[LJCoulomb],
                                           fixedParameters, forces, NULL, &totalEnergy );

      if( debug ){
         std::stringstream message;
         int maxAtom = 5;
         message << methodName << " Bonded + LJ/Coulomb forces:\n   ";
         for( int ii = 0; ii < maxAtom; ii++ ){
            message << ii;
            message << " x[";
            SimTKOpenMMUtilities::formatRealStringStream( message, atomCoordinates[ii], 3, one );
            message << "] f[";
            SimTKOpenMMUtilities::formatRealStringStream( message, forces[ii], 3, one );
            message << "]";
            message << std::endl;
         }
         SimTKOpenMMLog::printMessage( message );
      }
   
      // do update

      referenceStochasticDynamics.update( numberOfAtoms, atomCoordinates, velocities, forces, masses );

      if( debug ){
         std::stringstream message;
         message << methodName << " Called update: " << referenceStochasticDynamics.getTimeStep() << "\n";
         SimTKOpenMMLog::printMessage( message );
      }
   }

	referenceStochasticDynamics.writeState( numberOfAtoms, atomCoordinates,
                                           velocities, forces, masses,
                                           zeroI, baseFileName );

   // ---------------------------------------------------------------------------------------

   // free memory

   gromacsFreeSundryArrays( indexOneDArraysToFree, indexTwoDArraysToFree, indexTwoDIntArraysToFree,
                            MaxFreeArrayIndex, oneDArraysToFree, twoDArraysToFree, twoDIntArraysToFree,
                            methodName );

   for( int ii = 0; ii < MaxForceTypes; ii++ ){
      if( referenceBondIxn[ii] ){
         delete referenceBondIxn[ii];
      }
   }

   return 0;
}

/**---------------------------------------------------------------------------------------

   Remove linear momentum from velocities

   @param numberOfAtoms            number of atoms
   @param gromacsMasses            masses
   @param gromacsVelocities        velocities
   @param baseFileName             base file name
   @param log                      log reference (stdlog in md.c)

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsReferenceRemoveLinearMomentum( int numberOfAtoms, const real* gromacsMasses,
                                          rvec* gromacsVelocities,
                                          char* baseFileName, FILE* log ){ 

   // ---------------------------------------------------------------------------------------

   static const char* methodName               = "\ngromacsReferenceRemoveLinearMomentum: ";
   static int debug                            = false;

   static const RealOpenMM zero                = 0.0;
   static const RealOpenMM one                 = 1.0;

   static const int  zeroI                     = 0;
   static const int  twoI                      = 2;
   static const int  threeI                    = 3;
   static const int  fourI                     = 4;

   static const int  MaxFreeArrayIndex         = 100;

   int indexOneDArraysToFree                   = 0;
   RealOpenMM* oneDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDArraysToFree                   = 0;
   RealOpenMM** twoDArraysToFree[MaxFreeArrayIndex];

   int indexTwoDIntArraysToFree                = 0;
   int** twoDIntArraysToFree[MaxFreeArrayIndex];

   // ---------------------------------------------------------------------------------------

   // set log file, if available

   if( log ){
      SimTKOpenMMLog::setSimTKOpenMMLog( log );
   }

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << " atoms=" << numberOfAtoms;
      SimTKOpenMMLog::printMessage( message );
   }

   // ---------------------------------------------------------------------------------------

   // load masses and velocities

   RealOpenMM*  masses                       = SimTKOpenMMUtilities::allocateOneDRealOpenMMArray( numberOfAtoms, NULL, true, (RealOpenMM) 0.0 );
   oneDArraysToFree[indexOneDArraysToFree++] = masses;

   RealOpenMM** velocities                   = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
   twoDArraysToFree[indexTwoDArraysToFree++] = velocities;

   for( int ii = 0; ii < numberOfAtoms; ii++ ){

      masses[ii]             = gromacsMasses[ii];

      velocities[ii][0]      = gromacsVelocities[ii][0];
      velocities[ii][1]      = gromacsVelocities[ii][1];
      velocities[ii][2]      = gromacsVelocities[ii][2];
   }

   // ---------------------------------------------------------------------------------------

   // remove linear momentum, writing out velocities before and after

   ReferenceDynamics referenceDynamics( numberOfAtoms, zero, zero );
	referenceDynamics.writeState( numberOfAtoms, NULL, velocities, NULL, masses, zeroI, baseFileName );
   referenceDynamics.removeTotalLinearMomentum( numberOfAtoms, masses, velocities );
   referenceDynamics.incrementTimeStep();
	referenceDynamics.writeState( numberOfAtoms, NULL, velocities, NULL, masses, zeroI, baseFileName );

   // ---------------------------------------------------------------------------------------

   // free memory

   gromacsFreeSundryArrays( indexOneDArraysToFree, indexTwoDArraysToFree, indexTwoDIntArraysToFree,
                            MaxFreeArrayIndex, oneDArraysToFree, twoDArraysToFree, twoDIntArraysToFree,
                            methodName );

   // ---------------------------------------------------------------------------------------

   return 0;
}

/**---------------------------------------------------------------------------------------

   Calculate forces for given configuration and topology
   Used in self-test; parameter arrays, ... only allocated once (static); freed
   if debug is -1

   @param top                      Gromacs t_topology data struct
   @param gromacAtomCoordinates    atom configuration
   @param md                       Gromacs t_mdatoms data struct
   @param fr                       Gromacs t_forcerec data struct
   @param forces                   output forces
   @param includeNonBonded         include nonbonded ixn
   @param baseFileName             Base file name
   @param debug                    debug flag (== -1, free arrays and return)
   @param log                      log reference (stdlog in md.c)

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsReferenceForceSelfTest( 
                           const t_topology* top, const rvec* gromacAtomCoordinates,
                           const t_mdatoms *md, t_forcerec* fr, rvec* forces,
                           int includeNonBonded, char* baseFileName, int debug, FILE* log ){ 

   // ---------------------------------------------------------------------------------------

   static const char* methodName                     = "\ngromacsReferenceForceSelfTest: ";

   // ---------------------------------------------------------------------------------------

   static const RealOpenMM zero                      = 0.0;
   static const RealOpenMM one                       = 1.0;

   static const int  zeroI                           = 0;
   static const int  twoI                            = 2;
   static const int  threeI                          = 3;
   static const int  fourI                           = 4;

   // ---------------------------------------------------------------------------------------

   static const int  includeHarmonic                 = 1;
   static const int  includeAngle                    = 1;
   static const int  includeProperDihedral           = 1;
   static const int  includeRbDihedral               = 1;
   static const int  includeLJ14                     = 1;
   static const int  includeLJCoulomb                = 1;

   // ---------------------------------------------------------------------------------------

   static const int  HarmonicIndex                   = 0;
   static const int  AngleIndex                      = 1;
   static const int  ProperDihedralIndex             = 2;
   static const int  RbDihedralIndex                 = 3;
   static const int  LJ14Index                       = 4;
   static const int  LJCoulombIndex                  = 5;
   static const int  MaxForceIndex                   = 6;
   static int bondCount[MaxForceIndex];

   static const int  outputParameterFile             = 0;
   static int        outputForceFile                 = 1;

   RealOpenMM totalEnergy                            = zero;

   // ---------------------------------------------------------------------------------------

   static const int  MaxFreeArrayIndex               = 20;
   static int indexOneDArraysToFree                  = 0;
   static RealOpenMM* oneDArraysToFree[MaxFreeArrayIndex];

   static int indexTwoDArraysToFree                  = 0;
   static RealOpenMM** twoDArraysToFree[MaxFreeArrayIndex];

   static int indexTwoDIntArraysToFree               = 0;
   static int** twoDIntArraysToFree[MaxFreeArrayIndex];

   // ---------------------------------------------------------------------------------------

   static RealOpenMM** atomCoordinates               = NULL;
   static RealOpenMM** atomForces                    = NULL;
   static RealOpenMM** tempForces                    = NULL;
   static RealOpenMM*  tempEnergies                  = NULL;

   static int** harmonicBondAtomIndices              = NULL;
   static RealOpenMM** harmonicBondParameters        = NULL;

   static int** angleBondAtomIndices                 = NULL;
   static RealOpenMM** angleBondParameters           = NULL;

   static int** properDihedralBondAtomIndices        = NULL;
   static RealOpenMM** properDihedralBondParameters  = NULL;

   static int** rbDihedralBondAtomIndices            = NULL;
   static RealOpenMM** rbDihedralBondParameters      = NULL;

   static int** lj14AtomIndices                      = NULL;
   static RealOpenMM** lj14Parameters                = NULL;

   static int** ljCoulombExclusions                  = NULL;
   static RealOpenMM** ljCoulombParameters           = NULL;

   // ---------------------------------------------------------------------------------------

   // set log file, if available

   if( log ){
      SimTKOpenMMLog::setSimTKOpenMMLog( log );
   }

   if( !baseFileName ){
      outputForceFile = 0;
   }

   // ---------------------------------------------------------------------------------------

   // free arrays and return

   if( debug == -1 ){

      gromacsFreeSundryArrays( indexOneDArraysToFree, indexTwoDArraysToFree, indexTwoDIntArraysToFree,
                               MaxFreeArrayIndex, oneDArraysToFree, twoDArraysToFree, twoDIntArraysToFree,
                               methodName );

      atomCoordinates               = NULL;
      atomForces                    = NULL;
      tempForces                    = NULL;
      tempEnergies                  = NULL;
   
      harmonicBondAtomIndices       = NULL;
      harmonicBondParameters        = NULL;
   
      angleBondAtomIndices          = NULL;
      angleBondParameters           = NULL;
   
      properDihedralBondAtomIndices = NULL;
      properDihedralBondParameters  = NULL;
   
      rbDihedralBondAtomIndices     = NULL;
      rbDihedralBondParameters      = NULL;
   
      lj14AtomIndices               = NULL;
      lj14Parameters                = NULL;
   
      ljCoulombExclusions           = NULL;
      ljCoulombParameters           = NULL;
   
      return 0;
   }

   // ---------------------------------------------------------------------------------------

   int numberOfAtoms = top->atoms.nr;

   if( debug ){
      std::stringstream message;
      message << methodName;
      message << " atoms=" << numberOfAtoms;
      SimTKOpenMMLog::printMessage( message );
   }

   // setup

   if( atomCoordinates == NULL ){

      atomCoordinates                                     = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
      twoDArraysToFree[indexTwoDArraysToFree++]           = atomCoordinates;

      atomForces                                          = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
      twoDArraysToFree[indexTwoDArraysToFree++]           = atomForces;


      tempForces                                          = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
      twoDArraysToFree[indexTwoDArraysToFree++]           = tempForces;
   
      tempEnergies                                        = (RealOpenMM*) malloc( sizeof( RealOpenMM )*numberOfAtoms );
      oneDArraysToFree[indexOneDArraysToFree++]           = tempEnergies;
   
      // harmonic 

      bondCount[HarmonicIndex]                            = gromacsGetHarmonicBondParameters( top, &harmonicBondAtomIndices, &harmonicBondParameters );
      twoDIntArraysToFree[indexTwoDIntArraysToFree++]     = harmonicBondAtomIndices;
      twoDArraysToFree[indexTwoDArraysToFree++]           = harmonicBondParameters;

      // angle

      bondCount[AngleIndex]                               = gromacsGetAngleBondParameters( top, &angleBondAtomIndices, &angleBondParameters );
      twoDIntArraysToFree[indexTwoDIntArraysToFree++]     = angleBondAtomIndices;
      twoDArraysToFree[indexTwoDArraysToFree++]           = angleBondParameters;

      // proper dihedral

      bondCount[ProperDihedralIndex]                      = gromacsGetProperDihedralBondParameters( top, &properDihedralBondAtomIndices, &properDihedralBondParameters );
      twoDIntArraysToFree[indexTwoDIntArraysToFree++]     = properDihedralBondAtomIndices;
      twoDArraysToFree[indexTwoDArraysToFree++]           = properDihedralBondParameters;

      // rb dihedral

      bondCount[RbDihedralIndex]                          = gromacsGetRbDihedralBondParameters( top, &rbDihedralBondAtomIndices, &rbDihedralBondParameters );
      twoDIntArraysToFree[indexTwoDIntArraysToFree++]     = rbDihedralBondAtomIndices;
      twoDArraysToFree[indexTwoDArraysToFree++]           = rbDihedralBondParameters;
   
      // LJ14

      bondCount[LJ14Index]                                = gromacsGetLJ14Parameters( top, &lj14AtomIndices, &lj14Parameters );
      twoDIntArraysToFree[indexTwoDIntArraysToFree++]     = lj14AtomIndices;
      twoDArraysToFree[indexTwoDArraysToFree++]           = lj14Parameters;

      // calculate derived parameters

      ReferenceLJ14 referenceLJ14;
      RealOpenMM* charges = md->chargeA;
      for( int ii  = 0; ii < bondCount[LJ14Index]; ii++ ){
         referenceLJ14.getDerivedParameters( lj14Parameters[ii][0], lj14Parameters[ii][1],
                                             charges[lj14AtomIndices[ii][0]],
                                             charges[lj14AtomIndices[ii][1]], fr->epsfac*fr->fudgeQQ, lj14Parameters[ii] );
      }

      // LJ Coulomb

      if( includeNonBonded ){
         bondCount[LJCoulombIndex]                           = gromacsGetLJCoulombParameters( top, fr, md, &ljCoulombExclusions, &ljCoulombParameters );
         twoDIntArraysToFree[indexTwoDIntArraysToFree++]     = ljCoulombExclusions;
         twoDArraysToFree[indexTwoDArraysToFree++]           = ljCoulombParameters;
      
         ReferenceLJCoulombIxn referenceLJCoulomb;
      
         // calculate derived parameters
      
         RealOpenMM epsFacSqrt = SQRT( fr->epsfac );
         for( int ii  = 0; ii < numberOfAtoms; ii++ ){
            referenceLJCoulomb.getDerivedParameters( ljCoulombParameters[ii][0], ljCoulombParameters[ii][1],
                                                     ljCoulombParameters[ii][2],
                                                     epsFacSqrt, ljCoulombParameters[ii] );
         }
      }
   
   }
   memset( atomForces[0], 0, sizeof( RealOpenMM )*numberOfAtoms*3 );
   memset( tempEnergies, 0, sizeof( RealOpenMM )*numberOfAtoms );

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      atomCoordinates[ii][0] = gromacAtomCoordinates[ii][0];
      atomCoordinates[ii][1] = gromacAtomCoordinates[ii][1];
      atomCoordinates[ii][2] = gromacAtomCoordinates[ii][2];
   }

   // accumulate bonded force contributions

   // ---------------------------------------------------------------------------------------

   ReferenceBondForce referenceBondForce;

   // harmonic bond forces

   if( includeHarmonic ){

      ReferenceHarmonicBondIxn referenceHarmonicBondIxn;

      // harmonic bond forces

      memset( tempForces[0], 0, sizeof( RealOpenMM )*numberOfAtoms*3 );
      memset( tempEnergies, 0, sizeof( RealOpenMM )*numberOfAtoms );

      referenceBondForce.calculateForce( bondCount[HarmonicIndex], harmonicBondAtomIndices, atomCoordinates,
                                         harmonicBondParameters, tempForces,
                                         NULL, tempEnergies, 
                                         &totalEnergy, referenceHarmonicBondIxn );

      if( atomForces ){
         SimTKOpenMMUtilities::addTwoDimArray( numberOfAtoms, 3, tempForces, atomForces );
      }

      // output file
   
      if( outputForceFile ){
         std::stringstream harmonicBondForceFileName;
         harmonicBondForceFileName << baseFileName;
         harmonicBondForceFileName << "SelfHarmonicBond.txt";
         ReferenceForce::writeForces( numberOfAtoms, fourI, atomCoordinates, tempForces,
                                      tempEnergies, harmonicBondForceFileName.str(), one, one, one );
      }

   }

   // ---------------------------------------------------------------------------------------

   // angle bond forces

   if( includeAngle ){

      ReferenceAngleBondIxn referenceAngleBondIxn;

      // angle bond forces

      memset( tempForces[0], 0, sizeof( RealOpenMM )*numberOfAtoms*3 );
      memset( tempEnergies, 0, sizeof( RealOpenMM )*numberOfAtoms );

      referenceBondForce.calculateForce( bondCount[AngleIndex], angleBondAtomIndices, atomCoordinates,
                                         angleBondParameters, tempForces, 
                                         NULL, tempEnergies, 
                                         &totalEnergy, referenceAngleBondIxn );

      if( atomForces ){
         SimTKOpenMMUtilities::addTwoDimArray( numberOfAtoms, 3, tempForces, atomForces );
      }

      // output file
   
      if( outputForceFile ){
         std::stringstream angleForceFileName;
         angleForceFileName << baseFileName;
         angleForceFileName << "SelfAngleBond.txt";
         ReferenceForce::writeForces( numberOfAtoms, fourI, atomCoordinates, tempForces,
                                      tempEnergies, angleForceFileName.str(), one, one, one );
      }

   }

   // ---------------------------------------------------------------------------------------

   // properDihedral bond forces

   if( includeProperDihedral ){

      ReferenceProperDihedralBond referenceProperDihedralBond;

      memset( tempForces[0], 0, sizeof( RealOpenMM )*numberOfAtoms*3 );
      memset( tempEnergies, 0, sizeof( RealOpenMM )*numberOfAtoms );

      referenceBondForce.calculateForce( bondCount[ProperDihedralIndex], properDihedralBondAtomIndices, atomCoordinates,
                                         properDihedralBondParameters, tempForces, 
                                         NULL, tempEnergies, 
                                         &totalEnergy, referenceProperDihedralBond );

      // output file
   
      if( outputForceFile ){
         std::stringstream properDihedralForceFileName;
         properDihedralForceFileName << baseFileName;
         properDihedralForceFileName << "SelfProperDihedral.txt";
         ReferenceForce::writeForces( numberOfAtoms, fourI, atomCoordinates, tempForces,
                                      tempEnergies, properDihedralForceFileName.str(), one, one, one );
      }

      // accumulate bonded forces

      if( atomForces ){
         SimTKOpenMMUtilities::addTwoDimArray( numberOfAtoms, 3, tempForces, atomForces );
      }
   }

   // ---------------------------------------------------------------------------------------

   // rbDihedral bond forces

   if( includeRbDihedral ){

      ReferenceRbDihedralBond referenceRbDihedralBond;

      memset( tempForces[0], 0, sizeof( RealOpenMM )*numberOfAtoms*3 );
      memset( tempEnergies, 0, sizeof( RealOpenMM )*numberOfAtoms );

      referenceBondForce.calculateForce( bondCount[RbDihedralIndex], rbDihedralBondAtomIndices, atomCoordinates,
                                         rbDihedralBondParameters, tempForces, 
                                         NULL, tempEnergies, 
                                         &totalEnergy, referenceRbDihedralBond );

      // output file
   
      if( outputForceFile ){
         std::stringstream rbDihedralForceFileName;
         rbDihedralForceFileName << baseFileName;
         rbDihedralForceFileName << "SelfRbDihedral.txt";
         ReferenceForce::writeForces( numberOfAtoms, fourI, atomCoordinates, tempForces,
                                      tempEnergies, rbDihedralForceFileName.str(), one, one, one );
      }

      // accumulate bonded forces

      if( atomForces ){
         SimTKOpenMMUtilities::addTwoDimArray( numberOfAtoms, 3, tempForces, atomForces );
      }
   }

   // ---------------------------------------------------------------------------------------

   // lj14 bond forces

   if( includeLJ14 ){

      ReferenceLJ14 referenceLJ14;

      memset( tempForces[0], 0, sizeof( RealOpenMM )*numberOfAtoms*3 );
      memset( tempEnergies, 0, sizeof( RealOpenMM )*numberOfAtoms );

      referenceBondForce.calculateForce( bondCount[LJ14Index], lj14AtomIndices, atomCoordinates,
                                         lj14Parameters, tempForces, 
                                         NULL, tempEnergies, 
                                         &totalEnergy, referenceLJ14 );

      // output file
   
      if( outputForceFile ){
         std::stringstream lj14ForceFileName;
         lj14ForceFileName << baseFileName;
         lj14ForceFileName << "SelfLJ14.txt";
         ReferenceForce::writeForces( numberOfAtoms, fourI, atomCoordinates, tempForces,
                                      tempEnergies, lj14ForceFileName.str(), one, one, one );
      }

      // accumulate bonded forces

      if( atomForces ){
         SimTKOpenMMUtilities::addTwoDimArray( numberOfAtoms, 3, tempForces, atomForces );
      }
   }

   // output file

   if( outputForceFile && atomForces ){
      std::stringstream bondedForceFileName;
      bondedForceFileName << baseFileName;
      bondedForceFileName << "SelfBonded.txt";
      ReferenceForce::writeForces( numberOfAtoms, zeroI, atomCoordinates, atomForces,
                                   NULL, bondedForceFileName.str(), one, one, one );
   }

// ---------------------------------------------------------------------------------------

   // ljCoulomb forces

   if( includeNonBonded && includeLJCoulomb ){

      ReferenceLJCoulombIxn referenceLJCoulomb;

      // currently not used!

      RealOpenMM fixedParameters[MaxFreeArrayIndex];

      memset( tempForces[0], 0, sizeof( RealOpenMM )*numberOfAtoms*3 );
      memset( tempEnergies, 0, sizeof( RealOpenMM )*numberOfAtoms );

      // ljCoulomb bond forces
   
      referenceLJCoulomb.calculatePairIxn( numberOfAtoms, atomCoordinates,
                                           ljCoulombParameters, ljCoulombExclusions,
                                           fixedParameters, tempForces, 
                                           tempEnergies, &totalEnergy );
   
      if( outputForceFile ){
         std::stringstream lJCoulombForceFileName;
         lJCoulombForceFileName << baseFileName;
         lJCoulombForceFileName << "SelfLJCoulomb.txt";
         ReferenceForce::writeForces( numberOfAtoms, fourI, atomCoordinates, tempForces,
                                      tempEnergies, lJCoulombForceFileName.str(), one, one, one );
      }

      // accumulate forces

      if( atomForces ){
         SimTKOpenMMUtilities::addTwoDimArray( numberOfAtoms, 3, tempForces, atomForces );
      }
   }

   // output file

   if( outputForceFile && atomForces ){
      std::stringstream forceFileName;
      forceFileName << baseFileName;
      forceFileName << "SelfTotalF.txt";
      ReferenceForce::writeForces( numberOfAtoms, zeroI, atomCoordinates, atomForces,
                                   NULL, forceFileName.str(), one, one, one );
   }

   // ---------------------------------------------------------------------------------------

   // accumulate forces

   if( forces ){
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         forces[ii][0] = atomForces[ii][0];
         forces[ii][1] = atomForces[ii][1];
         forces[ii][2] = atomForces[ii][2];
      }
   }

   return 0;
}

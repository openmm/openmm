
/* Portions copyright (c) 2006-2009 Stanford University and Simbios.
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

#ifndef __ObcParameters_H__
#define __ObcParameters_H__

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "ImplicitSolventParameters.h"

// ---------------------------------------------------------------------------------------

class ObcParameters : public ImplicitSolventParameters {

   public:

       // OBC types

       enum ObcType { ObcTypeI, ObcTypeII };

       static const std::string ParameterFileName;

   private:

      // OBC constants & parameters
   
      RealOpenMM _dielectricOffset;
      RealOpenMM _alphaObc;
      RealOpenMM _betaObc;
      RealOpenMM _gammaObc;
      ObcType _obcType;

      // scaled radius factors (S_kk in HCT paper)

      int _ownScaledRadiusFactors;
      RealOpenMM* _scaledRadiusFactors;

      // cutoff and periodic boundary conditions
      
      bool cutoff;
      bool periodic;
      RealOpenMM periodicBoxSize[3];
      RealOpenMM cutoffDistance;


      /**---------------------------------------------------------------------------------------
      
         Set solvent dielectric (Simbios) 
      
         @param dielectricOffset         solvent dielectric

         --------------------------------------------------------------------------------------- */
      
      void setDielectricOffset( RealOpenMM dielectricOffset );

   public:

      /**---------------------------------------------------------------------------------------
      
         ObcParameters constructor (Simbios) 
      
         @param numberOfAtoms       number of atoms
      
         --------------------------------------------------------------------------------------- */
      
       ObcParameters( int numberOfAtoms, ObcParameters::ObcType obcType = ObcTypeII );

      /**---------------------------------------------------------------------------------------
      
         ObcParameters destructor (Simbios) 
      
         --------------------------------------------------------------------------------------- */
      
       ~ObcParameters( );

      // override of new/delete

      //static void* operator new( size_t size );
      //static void  operator delete( void *p );

      //static void* operator new[]( size_t size );
      //static void  operator delete[]( void *p );

      /**---------------------------------------------------------------------------------------
      
         Get OBC type
      
         @return OBC type
      
         --------------------------------------------------------------------------------------- */
      
      ObcParameters::ObcType getObcType( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set OBC type specific parameters
      
         @param obcType OBC type (ObcTypeI or ObcTypeII -- Eq. 7 or 8)
      
         --------------------------------------------------------------------------------------- */
      
      void setObcTypeParameters( ObcParameters::ObcType obcType );
      
      /**---------------------------------------------------------------------------------------
      
         Get alpha OBC (Eqs. 6 & 7) in Proteins paper
      
         @return alphaObc
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getAlphaObc( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get beta OBC (Eqs. 6 & 7) in Proteins paper
      
         @return betaObc
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getBetaObc( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get gamma OBC (Eqs. 6 & 7) in Proteins paper
      
         @return gammaObc
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getGammaObc( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get solvent dielectric (Simbios) 
      
         @return dielectricOffset dielectric offset
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getDielectricOffset( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Return OBC scale factors
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      const RealOpenMM* getScaledRadiusFactors( void ) const;
        
      /**---------------------------------------------------------------------------------------
      
         Return OBC scale factors
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      void setScaledRadiusFactors( RealOpenMM* scaledRadiusFactors );
      void setScaledRadiusFactors( const RealOpenMMVector& scaledRadiusFactors );
        
      /**---------------------------------------------------------------------------------------
      
         Set flag indicating whether scale factors arra should be deleted
      
         @param ownScaledRadiusFactors flag indicating whether scale factors 
                                       array should be deleted
      
         --------------------------------------------------------------------------------------- */

      void setOwnScaleFactors( int ownScaledRadiusFactors );
      
      /**--------------------------------------------------------------------------------------- 
      
         Assign standard radii for GB/SA methods other than ACE;
         taken from Macromodel and OPLS-AA, except for hydrogens (Simbios)
      
         Logic based on logic in Tinker's ksolv.f
      
         Currently only works for standard amino acid atoms
         If invalid atom name is encountered, a message is printed to log file and the
         radius for that atom is set to 1.0f
      
         @param numberOfAtoms       number of atoms
         @param atomNames           array of atom names from GMX top data struct
         @param radii               array to store Macromodel radii for each atom
         @param log                 if set, then print error messages to log file
      
         --------------------------------------------------------------------------------------- */
      
      void getMacroModelAtomicRadii( int numberOfAtoms,
                                    char*** atomNames, RealOpenMM* radii, FILE* log );

      /**---------------------------------------------------------------------------------------
      
         Get AtomicRadii array w/ dielectric offset applied
      
         @return array of atom volumes
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM* getAtomicRadii( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set AtomicRadii array
      
         @param atomicRadii array of atomic radii
      
         --------------------------------------------------------------------------------------- */

      void setAtomicRadii( RealOpenMM* atomicRadii );

      /**---------------------------------------------------------------------------------------
      
         Set AtomicRadii array
      
         @param atomicRadii vector of atomic radii
      
         --------------------------------------------------------------------------------------- */

      void setAtomicRadii( const RealOpenMMVector& atomicRadii );

      /**---------------------------------------------------------------------------------------
      
         Map Gmx atom name to Tinker atom number (Simbios)
      
         @param atomName            atom name (CA, HA, ...); upper and lower case should both work
         @param log                 if set, then print error messages to log file
      
         return Tinker atom number if atom name is valid; else return -1
      
         --------------------------------------------------------------------------------------- */
            
      int mapGmxAtomNameToTinkerAtomNumber( const char* atomName, FILE* log ) const;

      /**---------------------------------------------------------------------------------------
            
         Get string w/ state
         
         @param title               title (optional)
            
         @return string
            
         --------------------------------------------------------------------------------------- */
      
      std::string getStateString( const char* title ) const;

      /**---------------------------------------------------------------------------------------

         Set the force to use a cutoff.

         @param distance            the cutoff distance

         --------------------------------------------------------------------------------------- */

      void setUseCutoff( RealOpenMM distance );

      /**---------------------------------------------------------------------------------------

         Get whether to use a cutoff.

         --------------------------------------------------------------------------------------- */

      bool getUseCutoff();

      /**---------------------------------------------------------------------------------------

         Get the cutoff distance.

         --------------------------------------------------------------------------------------- */

      RealOpenMM getCutoffDistance();

      /**---------------------------------------------------------------------------------------

         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.

         @param boxSize             the X, Y, and Z widths of the periodic box

         --------------------------------------------------------------------------------------- */

      void setPeriodic( OpenMM::RealVec& boxSize );

      /**---------------------------------------------------------------------------------------

         Get whether to use periodic boundary conditions.

         --------------------------------------------------------------------------------------- */

      bool getPeriodic();

      /**---------------------------------------------------------------------------------------

         Get the periodic box dimension

         --------------------------------------------------------------------------------------- */

      const RealOpenMM* getPeriodicBox();

};
   
/**---------------------------------------------------------------------------------------
      
   Qsort/heapsort integer comparison (Simbios) 
      
   @parma a first value to compare
   @param b second value to compare

   @return -1, 0, 1
      
--------------------------------------------------------------------------------------- */

int integerComparison( const void *a, const void *b);

// ---------------------------------------------------------------------------------------

#endif // __ObcParameters_H__

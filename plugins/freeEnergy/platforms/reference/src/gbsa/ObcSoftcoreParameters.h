
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

#ifndef __ObcSoftcoreParameters_H__
#define __ObcSoftcoreParameters_H__

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "gbsa/ImplicitSolventParameters.h"

// ---------------------------------------------------------------------------------------

class ObcSoftcoreParameters : public ImplicitSolventParameters {

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

      RealOpenMM _nonPolarPreFactor;

     // scaling factors for nonpolar term

      int _ownNonPolarScaleFactors;
      RealOpenMM* _nonPolarScaleFactors;

      // scaled radius factors (S_kk in HCT paper)

      int _ownScaledRadiusFactors;
      RealOpenMM* _scaledRadiusFactors;

      // cutoff and periodic boundary conditions
      
      bool cutoff;
      bool periodic;
      RealOpenMM periodicBoxSize[3];
      RealOpenMM cutoffDistance;

   public:

      /**---------------------------------------------------------------------------------------
      
         ObcSoftcoreParameters constructor (Simbios) 
      
         @param numberOfAtoms       number of atoms
      
         --------------------------------------------------------------------------------------- */
      
       ObcSoftcoreParameters( int numberOfAtoms, ObcSoftcoreParameters::ObcType obcType = ObcTypeII );

      /**---------------------------------------------------------------------------------------
      
         ObcSoftcoreParameters destructor (Simbios) 
      
         --------------------------------------------------------------------------------------- */
      
       ~ObcSoftcoreParameters( );

      /**---------------------------------------------------------------------------------------
      
         Get OBC type
      
         @return OBC type
      
         --------------------------------------------------------------------------------------- */
      
      ObcSoftcoreParameters::ObcType getObcType( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set OBC type specific parameters
      
         @param obcType OBC type (ObcTypeI or ObcTypeII -- Eq. 7 or 8)
      
         --------------------------------------------------------------------------------------- */
      
      void setObcTypeParameters( ObcSoftcoreParameters::ObcType obcType );
      
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
#if RealOpenMMType == 0
      void setScaledRadiusFactors( float* scaledRadiusFactors );
#endif
      void setScaledRadiusFactors( const RealOpenMMVector& scaledRadiusFactors );
        
      /**---------------------------------------------------------------------------------------
      
         Set flag indicating whether scale factors arra should be deleted
      
         @param ownScaledRadiusFactors flag indicating whether scale factors 
                                       array should be deleted
      
         --------------------------------------------------------------------------------------- */

      void setOwnScaleFactors( int ownScaledRadiusFactors );

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

      void setPeriodic( RealOpenMM* boxSize );

      /**---------------------------------------------------------------------------------------

         Get whether to use periodic boundary conditions.

         --------------------------------------------------------------------------------------- */

      bool getPeriodic();

      /**---------------------------------------------------------------------------------------

         Get the periodic box dimension

         --------------------------------------------------------------------------------------- */

      const RealOpenMM* getPeriodicBox();

      /**---------------------------------------------------------------------------------------
      
         Set flag indicating whether scale factors array should be deleted
      
         @param ownNonPolarScaleFactors flag indicating whether scale factors 
                                       array should be deleted
      
         --------------------------------------------------------------------------------------- */
      
      void setOwnNonPolarScaleFactors( int ownNonPolarScaleFactors );
      
      /**---------------------------------------------------------------------------------------
      
         Return non-polar scale factors
         If not previously set, allocate space
      
         @return array 
      
         --------------------------------------------------------------------------------------- */
      
      const RealOpenMM* getNonPolarScaleFactors( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set non-polar scale factors

         @param nonPolarScaleFactors nonpolar scale factors
      
         --------------------------------------------------------------------------------------- */
      
      void setNonPolarScaleFactors( const RealOpenMMVector& nonPolarScaleFactors );

      /**---------------------------------------------------------------------------------------
      
         Set nonPolarPrefactor
      
         @param nonPolarPrefactor solute dielectric
      
         --------------------------------------------------------------------------------------- */

      void setNonPolarPrefactor( RealOpenMM nonPolarPrefactor );

};
   
/**---------------------------------------------------------------------------------------
      
   Qsort/heapsort integer comparison (Simbios) 
      
   @parma a first value to compare
   @param b second value to compare

   @return -1, 0, 1
      
--------------------------------------------------------------------------------------- */

int integerComparison( const void *a, const void *b);

// ---------------------------------------------------------------------------------------

#endif // __ObcSoftcoreParameters_H__

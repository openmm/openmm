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

// ---------------------------------------------------------------------------------------

class ObcSoftcoreParameters {

   public:

       // OBC types

       enum ObcType { ObcTypeI, ObcTypeII };

   private:

      int _numberOfAtoms;

      // OBC constants & parameters
   
      RealOpenMM _dielectricOffset;
      RealOpenMM _alphaObc;
      RealOpenMM _betaObc;
      RealOpenMM _gammaObc;
      RealOpenMM _probeRadius;
      RealOpenMM _pi4Asolv;

      ObcType _obcType;

      RealOpenMM _nonPolarPreFactor;

      RealOpenMM _solventDielectric;
      RealOpenMM _soluteDielectric;
      RealOpenMM _electricConstant;

     // scaling factors for nonpolar term

      RealOpenMMVector _atomicRadii;
      RealOpenMMVector _nonPolarScaleFactors;

      // scaled radius factors (S_kk in HCT paper)

      RealOpenMMVector _scaledRadiusFactors;

      // cutoff and periodic boundary conditions
      
      bool _cutoff;
      bool _periodic;
      RealOpenMM _periodicBoxSize[3];
      RealOpenMM _cutoffDistance;

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
      
         Get number of atoms
      
         @return number of atoms
      
         --------------------------------------------------------------------------------------- */

      int getNumberOfAtoms( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Get electric constant
      
         @return electric constant
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM getElectricConstant( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Get probe radius (Simbios) 
      
         @return probeRadius
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM getProbeRadius( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set probe radius (Simbios) 
      
         @param probeRadius   probe radius
      
         --------------------------------------------------------------------------------------- */

      void setProbeRadius( RealOpenMM probeRadius );

      /**---------------------------------------------------------------------------------------
      
         Get pi4Asolv:  used in ACE approximation for nonpolar term  
            ((RealOpenMM) M_PI)*4.0f*0.0049f*1000.0f; (Simbios) 
      
         @return pi4Asolv
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM getPi4Asolv( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Get solvent dielectric
      
         @return solvent dielectric
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM getSolventDielectric( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set solvent dielectric
      
         @param solventDielectric solvent dielectric
      
         --------------------------------------------------------------------------------------- */

      void setSolventDielectric( RealOpenMM solventDielectric );

      /**---------------------------------------------------------------------------------------
      
         Get solute dielectric
      
         @return soluteDielectric
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM getSoluteDielectric( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set solute dielectric
      
         @param soluteDielectric solute dielectric
      
         --------------------------------------------------------------------------------------- */

      void setSoluteDielectric( RealOpenMM soluteDielectric );

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
      
      const RealOpenMMVector& getScaledRadiusFactors( void ) const;
        
      /**---------------------------------------------------------------------------------------
      
         Set OBC scale factors
      
         @param input vector of radius factors
      
         --------------------------------------------------------------------------------------- */
      
      void setScaledRadiusFactors( const RealOpenMMVector& scaledRadiusFactors );
        
      /**---------------------------------------------------------------------------------------
      
         Get AtomicRadii array w/ dielectric offset applied
      
         @return array of atom volumes
      
         --------------------------------------------------------------------------------------- */

      const RealOpenMMVector& getAtomicRadii( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set AtomicRadii array
      
         @param atomicRadii vector of atomic radii
      
         --------------------------------------------------------------------------------------- */

      void setAtomicRadii( const RealOpenMMVector& atomicRadii );

      /**---------------------------------------------------------------------------------------

         Set the force to use a cutoff.

         @param distance            the cutoff distance

         --------------------------------------------------------------------------------------- */

      void setUseCutoff( RealOpenMM distance );

      /**---------------------------------------------------------------------------------------

         Get whether to use a cutoff.

         --------------------------------------------------------------------------------------- */

      bool getUseCutoff( void );

      /**---------------------------------------------------------------------------------------

         Get the cutoff distance.

         --------------------------------------------------------------------------------------- */

      RealOpenMM getCutoffDistance( void );

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

      bool getPeriodic( void );

      /**---------------------------------------------------------------------------------------

         Get the periodic box dimension

         --------------------------------------------------------------------------------------- */

      const RealOpenMM* getPeriodicBox( void );

      /**---------------------------------------------------------------------------------------
      
         Return non-polar scale factors
      
         @return array 
      
         --------------------------------------------------------------------------------------- */
      
      const RealOpenMMVector& getNonPolarScaleFactors( void ) const;

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
   
// ---------------------------------------------------------------------------------------

#endif // __ObcSoftcoreParameters_H__

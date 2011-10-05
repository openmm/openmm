
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

#ifndef __GBVISoftcoreParameters_H__
#define __GBVISoftcoreParameters_H__

#include "SimTKUtilities/SimTKOpenMMCommon.h"

// ---------------------------------------------------------------------------------------

class GBVISoftcoreParameters {

   public:

    /** 
     * This is an enumeration of the different methods that may be used for scaling of the Born radii.
     */
    enum BornRadiusScalingSoftcoreMethod {
        /**
         * No scaling method is applied.
         */
        NoScaling          = 0,
        /**
         * Use quintic spline scaling function
         */
        QuinticSpline       = 1
    };  

   private:

      int _numberOfAtoms;

      // parameters:
      //    scaled radii
      //    gamma parameters
      //    BornRadiusScaleFactors parameters

      RealOpenMMVector _scaledRadii;
      RealOpenMMVector _atomicRadii;
      RealOpenMMVector _gammaParameters;
      RealOpenMMVector _bornRadiusScaleFactors;

      RealOpenMM _solventDielectric;
      RealOpenMM _soluteDielectric;
      RealOpenMM _electricConstant;

      // cutoff and periodic boundary conditions
      
      bool _cutoff;
      bool _periodic;
      RealOpenMM _periodicBoxSize[3];
      RealOpenMM _cutoffDistance;

      // Born radii switching function params

      BornRadiusScalingSoftcoreMethod _bornRadiusScalingSoftcoreMethod; 

      RealOpenMM _quinticLowerLimitFactor;
      RealOpenMM _quinticUpperBornRadiusLimit;
      RealOpenMM _quinticUpperSplineLimit;

   public:

      /**---------------------------------------------------------------------------------------
      
         GBVISoftcoreParameters constructor (Simbios) 
      
         @param numberOfAtoms       number of atoms
      
         --------------------------------------------------------------------------------------- */
      
       GBVISoftcoreParameters( int numberOfAtoms );

      /**---------------------------------------------------------------------------------------
      
         GBVISoftcoreParameters destructor (Simbios) 
      
         --------------------------------------------------------------------------------------- */
      
       ~GBVISoftcoreParameters( );

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
      
         Return scaled radii
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      const RealOpenMMVector& getScaledRadii( void ) const;
        
      /**---------------------------------------------------------------------------------------
      
         Set scaled radii
      
         @param vector of scaled radii
      
         --------------------------------------------------------------------------------------- */
      
      void setScaledRadii( const RealOpenMMVector& radii );
        
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
      
         Get GammaParameters array
      
         @return array of gamma values
      
         --------------------------------------------------------------------------------------- */

      const RealOpenMMVector& getGammaParameters( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set GammaParameters array
      
         @param gammaParameters   array of gamma parameters
      
         --------------------------------------------------------------------------------------- */

      void setGammaParameters( const RealOpenMMVector& gammaParameters );

      /**---------------------------------------------------------------------------------------
      
         Get BornRadiusScaleFactors array
      
         @return array of bornRadiusScaleFactor values
      
         --------------------------------------------------------------------------------------- */

      const RealOpenMMVector& getBornRadiusScaleFactors( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set BornRadiusScaleFactors array
      
         @param bornRadiusScaleFactors   array of bornRadiusScaleFactor parameters
      
         --------------------------------------------------------------------------------------- */

      void setBornRadiusScaleFactors( const RealOpenMMVector& bornRadiusScaleFactors );

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

         Get the quintic spline lower limit factor 

         --------------------------------------------------------------------------------------- */

      RealOpenMM getQuinticLowerLimitFactor( void ) const;

      /**---------------------------------------------------------------------------------------

         Set the quintic spline lower limit factor 

         --------------------------------------------------------------------------------------- */

      void setQuinticLowerLimitFactor( RealOpenMM quinticLowerLimitFactor );

      /**---------------------------------------------------------------------------------------

         Get the quintic spline upper limit

         --------------------------------------------------------------------------------------- */

      RealOpenMM getQuinticUpperBornRadiusLimit( void ) const;
      RealOpenMM getQuinticUpperSplineLimit( void ) const;

      /**---------------------------------------------------------------------------------------

         Set the quintic spline upper limit

         --------------------------------------------------------------------------------------- */

      void setQuinticUpperBornRadiusLimit( RealOpenMM quinticUpperBornRadiusLimit );

      /**---------------------------------------------------------------------------------------
      
         Get tau prefactor
      
         @return (1/e1 - 1/e0), where e1 = solute dielectric, e0 = solvent dielectric
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getTau( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Get BornRadiusScalingSoftcoreMethod
      
         @return scaling method
      
         --------------------------------------------------------------------------------------- */
      
      BornRadiusScalingSoftcoreMethod getBornRadiusScalingSoftcoreMethod( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set BornRadiusScalingSoftcoreMethod
      
         @param scaling method
      
         --------------------------------------------------------------------------------------- */
      
      void setBornRadiusScalingSoftcoreMethod( BornRadiusScalingSoftcoreMethod method );

};
   
#endif // __GBVISoftcoreParameters_H__

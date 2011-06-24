
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

#ifndef __GBVIParameters_H__
#define __GBVIParameters_H__

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "ImplicitSolventParameters.h"

// ---------------------------------------------------------------------------------------

class GBVIParameters : public ImplicitSolventParameters {

   public:

       static const std::string ParameterFileName;

    /** 
     * This is an enumeration of the different methods that may be used for scaling of the Born radii.
     */
    enum BornRadiusScalingMethod {
        /**
         * No scaling method is applied.
         */
        NoScaling          = 0,
        /**
         * Use the method outlined in Proteins 55, 383-394 (2004), Eq. 6
         */
        Tanh               = 2,
        /**
         * Use quintic spline scaling function
         */
        QuinticSpline       = 1
    };  

   private:

      // scaled radii

      int _ownScaledRadii;
      RealOpenMM* _scaledRadii;

      // gamma parameters
      int _ownGammaParameters;
      RealOpenMM* _gammaParameters;

      // cutoff and periodic boundary conditions
      
      bool cutoff;
      bool periodic;
      RealOpenMM periodicBoxSize[3];
      RealOpenMM cutoffDistance;

      int _bornRadiusScalingMethod;
      RealOpenMM _quinticLowerLimitFactor;
      RealOpenMM _quinticUpperBornRadiusLimit;

   public:

      /**---------------------------------------------------------------------------------------
      
         GBVIParameters constructor (Simbios) 
      
         @param numberOfAtoms       number of atoms
      
         --------------------------------------------------------------------------------------- */
      
       GBVIParameters( int numberOfAtoms );

      /**---------------------------------------------------------------------------------------
      
         GBVIParameters destructor (Simbios) 
      
         --------------------------------------------------------------------------------------- */
      
       ~GBVIParameters( );

      /**---------------------------------------------------------------------------------------
      
         Return scaled radii
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      const RealOpenMM* getScaledRadii( void ) const;
        
      /**---------------------------------------------------------------------------------------
      
         Return scaled radii
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      void setScaledRadii( RealOpenMM* scaledRadii );
      void setScaledRadii( const RealOpenMMVector& scaledRadii );
        
      /**---------------------------------------------------------------------------------------
      
         Set flag indicating whether scaled radii array should be deleted
      
         @param ownScaledRadiusFactors flag indicating whether scaled radii
                                       array should be deleted
      
         --------------------------------------------------------------------------------------- */

      void setOwnScaledRadii( int ownScaledRadii );
      
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
      
         Set flag indicating whether gamma parameter array should be deleted
      
         @param ownGammaParameters   flag indicating whether gamma parameter
                                     array should be deleted
      
         --------------------------------------------------------------------------------------- */

      void setOwnGammaParameters( int ownGammaParameters );
      
      /**---------------------------------------------------------------------------------------
      
         Get GammaParameters array
      
         @return array of gamma values
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM* getGammaParameters( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set GammaParameters array
      
         @param gammaParameters    array of gamma parameters
      
         --------------------------------------------------------------------------------------- */

      void setGammaParameters( RealOpenMM* gammaParameters );

      /**---------------------------------------------------------------------------------------
      
         Set GammaParameters array
      
         @param gammaParameters   array of gamma parameters
      
         --------------------------------------------------------------------------------------- */

      void setGammaParameters( const RealOpenMMVector& gammaParameters );

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

      /**---------------------------------------------------------------------------------------
      
         Get tau prefactor
      
         @return (1/e1 - 1/e0), where e1 = solute dielectric, e0 = solvent dielectric
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getTau( void ) const;
      
    /**---------------------------------------------------------------------------------------
    
       Get bornRadiusScalingMethod
    
       @return bornRadiusScalingMethod
    
       --------------------------------------------------------------------------------------- */
    
    int getBornRadiusScalingMethod( void ) const;
    
    /**---------------------------------------------------------------------------------------
    
       Set bornRadiusScalingMethod 
    
       @param bornRadiusScalingMethod bornRadiusScalingMethod
    
       --------------------------------------------------------------------------------------- */
    
    void setBornRadiusScalingMethod( int bornRadiusScalingMethod );
    
    /**---------------------------------------------------------------------------------------
    
       Get quinticLowerLimitFactor
    
       @return quinticLowerLimitFactor
    
       --------------------------------------------------------------------------------------- */
    
    RealOpenMM getQuinticLowerLimitFactor( void ) const;
    
    /**---------------------------------------------------------------------------------------
    
       Set quinticLowerLimitFactor 
    
       @param quinticLowerLimitFactor quinticLowerLimitFactor
    
       --------------------------------------------------------------------------------------- */
    
    void setQuinticLowerLimitFactor( RealOpenMM quinticLowerLimitFactor );
    
    /**---------------------------------------------------------------------------------------
    
       Get quinticUpperBornRadiusLimit
    
       @return quinticUpperBornRadiusLimit
    
       --------------------------------------------------------------------------------------- */
    
    RealOpenMM getQuinticUpperBornRadiusLimit( void ) const;
    
    /**---------------------------------------------------------------------------------------
    
       Set quinticUpperBornRadiusLimit 
    
       @param quinticUpperBornRadiusLimit quinticUpperBornRadiusLimit
    
       --------------------------------------------------------------------------------------- */
    
    void setQuinticUpperBornRadiusLimit( RealOpenMM quinticUpperSplineLimit );

};
   
#endif // __GBVIParameters_H__

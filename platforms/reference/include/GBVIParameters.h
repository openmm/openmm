
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

#include "RealVec.h"
#include <vector>

namespace OpenMM {

class GBVIParameters {

   public:

    /** 
     * This is an enumeration of the different methods that may be used for scaling of the Born radii.
     */
    enum BornRadiusScalingMethod {
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

      RealOpenMM _solventDielectric;
      RealOpenMM _soluteDielectric;
      RealOpenMM _electricConstant;

      // parameter vectors

      std::vector<RealOpenMM> _atomicRadii;
      std::vector<RealOpenMM> _scaledRadii;
      std::vector<RealOpenMM> _gammaParameters;

      // cutoff and periodic boundary conditions
      
      bool _cutoff;
      bool _periodic;
      OpenMM::RealVec _periodicBoxVectors[3];
      RealOpenMM _cutoffDistance;

      int _bornRadiusScalingMethod;
      RealOpenMM _quinticLowerLimitFactor;
      RealOpenMM _quinticUpperBornRadiusLimit;

   public:

      /**---------------------------------------------------------------------------------------
      
         GBVIParameters constructor (Simbios) 
      
         @param numberOfAtoms       number of atoms
      
         --------------------------------------------------------------------------------------- */
      
       GBVIParameters(int numberOfAtoms);

      /**---------------------------------------------------------------------------------------
      
         GBVIParameters destructor (Simbios) 
      
         --------------------------------------------------------------------------------------- */
      
       ~GBVIParameters();

      /**---------------------------------------------------------------------------------------
      
         Get number of atoms
      
         @return number of atoms
      
         --------------------------------------------------------------------------------------- */

      int getNumberOfAtoms() const;

      /**---------------------------------------------------------------------------------------
      
         Get electric constant
      
         @return electric constant
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM getElectricConstant() const;

      /**---------------------------------------------------------------------------------------
      
         Get solvent dielectric
      
         @return solvent dielectric
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM getSolventDielectric() const;

      /**---------------------------------------------------------------------------------------
      
         Set solvent dielectric
      
         @param solventDielectric solvent dielectric
      
         --------------------------------------------------------------------------------------- */

      void setSolventDielectric(RealOpenMM solventDielectric);

      /**---------------------------------------------------------------------------------------
      
         Get solute dielectric
      
         @return soluteDielectric
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM getSoluteDielectric() const;

      /**---------------------------------------------------------------------------------------
      
         Set solute dielectric
      
         @param soluteDielectric solute dielectric
      
         --------------------------------------------------------------------------------------- */

      void setSoluteDielectric(RealOpenMM soluteDielectric);

      /**---------------------------------------------------------------------------------------
      
         Return scaled radii
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      const std::vector<RealOpenMM>& getScaledRadii() const;
        
      /**---------------------------------------------------------------------------------------
      
         Return scaled radii
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      void setScaledRadii(const std::vector<RealOpenMM>& scaledRadii);
        
      /**---------------------------------------------------------------------------------------
      
         Get AtomicRadii array w/ dielectric offset applied
      
         @return array of atom volumes
      
         --------------------------------------------------------------------------------------- */

      const std::vector<RealOpenMM>& getAtomicRadii() const;

      /**---------------------------------------------------------------------------------------
      
         Set AtomicRadii array
      
         @param atomicRadii vector of atomic radii
      
         --------------------------------------------------------------------------------------- */

      void setAtomicRadii(const std::vector<RealOpenMM>& atomicRadii);

      /**---------------------------------------------------------------------------------------
      
         Get GammaParameters array
      
         @return array of gamma values
      
         --------------------------------------------------------------------------------------- */

      const std::vector<RealOpenMM>& getGammaParameters() const;

      /**---------------------------------------------------------------------------------------
      
         Set GammaParameters array
      
         @param gammaParameters   array of gamma parameters
      
         --------------------------------------------------------------------------------------- */

      void setGammaParameters(const std::vector<RealOpenMM>& gammaParameters);

      /**---------------------------------------------------------------------------------------

         Set the force to use a cutoff.

         @param distance            the cutoff distance

         --------------------------------------------------------------------------------------- */

      void setUseCutoff(RealOpenMM distance);

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

         @param vectors    the vectors defining the periodic box

         --------------------------------------------------------------------------------------- */

      void setPeriodic(OpenMM::RealVec* vectors);

      /**---------------------------------------------------------------------------------------

         Get whether to use periodic boundary conditions.

         --------------------------------------------------------------------------------------- */

      bool getPeriodic();

      /**---------------------------------------------------------------------------------------

         Get the periodic box dimension

         --------------------------------------------------------------------------------------- */

      const OpenMM::RealVec* getPeriodicBox();

      /**---------------------------------------------------------------------------------------
      
         Get tau prefactor
      
         @return (1/e1 - 1/e0), where e1 = solute dielectric, e0 = solvent dielectric
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getTau() const;
      
    /**---------------------------------------------------------------------------------------
    
       Get bornRadiusScalingMethod
    
       @return bornRadiusScalingMethod
    
       --------------------------------------------------------------------------------------- */
    
    int getBornRadiusScalingMethod() const;
    
    /**---------------------------------------------------------------------------------------
    
       Set bornRadiusScalingMethod 
    
       @param bornRadiusScalingMethod bornRadiusScalingMethod
    
       --------------------------------------------------------------------------------------- */
    
    void setBornRadiusScalingMethod(int bornRadiusScalingMethod);
    
    /**---------------------------------------------------------------------------------------
    
       Get quinticLowerLimitFactor
    
       @return quinticLowerLimitFactor
    
       --------------------------------------------------------------------------------------- */
    
    RealOpenMM getQuinticLowerLimitFactor() const;
    
    /**---------------------------------------------------------------------------------------
    
       Set quinticLowerLimitFactor 
    
       @param quinticLowerLimitFactor quinticLowerLimitFactor
    
       --------------------------------------------------------------------------------------- */
    
    void setQuinticLowerLimitFactor(RealOpenMM quinticLowerLimitFactor);
    
    /**---------------------------------------------------------------------------------------
    
       Get quinticUpperBornRadiusLimit
    
       @return quinticUpperBornRadiusLimit
    
       --------------------------------------------------------------------------------------- */
    
    RealOpenMM getQuinticUpperBornRadiusLimit() const;
    
    /**---------------------------------------------------------------------------------------
    
       Set quinticUpperBornRadiusLimit 
    
       @param quinticUpperBornRadiusLimit quinticUpperBornRadiusLimit
    
       --------------------------------------------------------------------------------------- */
    
    void setQuinticUpperBornRadiusLimit(RealOpenMM quinticUpperSplineLimit);

};

} // namespace OpenMM

#endif // __GBVIParameters_H__

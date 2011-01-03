
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

#ifndef __CpuGBVISoftcore_H__
#define __CpuGBVISoftcore_H__

#include "GBVISoftcoreParameters.h"
#include "gbsa/CpuImplicitSolvent.h"

// ---------------------------------------------------------------------------------------

class CpuGBVISoftcore : public CpuImplicitSolvent {

   private:

      // GB/VI parameters

      GBVISoftcoreParameters* _gbviParameters;

      // arrays containing switching function derivative 

      RealOpenMM* _switchDeriviative;

      // initialize data members (more than
      // one constructor, so centralize intialization here)

      void _initializeGBVISoftcoreDataMembers( void );

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         @param implicitSolventParameters    ImplicitSolventParameters reference
      
         @return CpuImplicitSolvent object
      
         --------------------------------------------------------------------------------------- */

       CpuGBVISoftcore( ImplicitSolventParameters* gbviParameters );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~CpuGBVISoftcore( );

      /**---------------------------------------------------------------------------------------
      
         Return GBVISoftcoreParameters
      
         @return GBVISoftcoreParameters
      
         --------------------------------------------------------------------------------------- */

      GBVISoftcoreParameters* getGBVISoftcoreParameters( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set ImplicitSolventParameters
      
         @param ImplicitSolventParameters
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */

      int setGBVISoftcoreParameters( GBVISoftcoreParameters* gbviParameters );
 
      /**---------------------------------------------------------------------------------------
      
         Return OBC chain derivative: size = _implicitSolventParameters->getNumberOfAtoms()
         On first call, memory for array is allocated if not set
      
         @return array
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM* getSwitchDeriviative( void );
      RealOpenMM* getSwitchDeriviativeConst( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Get Born radii based on Eq. 3 of Labute paper [JCC 29 p. 1693-1698 2008])

         @param atomCoordinates   atomic coordinates
         @param bornRadii         output array of Born radii
         @param switchDeriviative derivative of switch function
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int computeBornRadii( std::vector<OpenMM::RealVec>& atomCoordinates, RealOpenMM* bornRadii,
                            RealOpenMM* switchDeriviative = NULL );
      
      /**---------------------------------------------------------------------------------------
      
         Get state 
      
         title             title (optional)
      
         @return state string
      
         --------------------------------------------------------------------------------------- */
      
      std::string getStateString( const char* title ) const;

      /**---------------------------------------------------------------------------------------
      
         Get volume Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])
      
         @param r                   distance between atoms i & j
         @param R                   atomic radius
         @param S                   scaled atomic radius
      
         @return volume
      
         --------------------------------------------------------------------------------------- */
      
      static RealOpenMM getVolume( RealOpenMM r, RealOpenMM R, RealOpenMM S );
      
      /**---------------------------------------------------------------------------------------
      
         Get L (analytical solution for volume integrals) 
      
         @param r                   distance between atoms i & j
         @param R                   atomic radius
         @param S                   scaled atomic radius
      
         @return L value (Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])
      
         --------------------------------------------------------------------------------------- */
      
      static RealOpenMM getL( RealOpenMM r, RealOpenMM x, RealOpenMM S );

      /**---------------------------------------------------------------------------------------
      
         Get partial derivative of L wrt r
      
         @param r                   distance between atoms i & j
         @param R                   atomic radius
         @param S                   scaled atomic radius
      
         @return partial derivative based on Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])
      
         --------------------------------------------------------------------------------------- */
      
      static RealOpenMM dL_dr( RealOpenMM r, RealOpenMM x, RealOpenMM S );

      /**---------------------------------------------------------------------------------------
      
         Get partial derivative of L wrt x
      
         @param r                   distance between atoms i & j
         @param R                   atomic radius
         @param S                   scaled atomic radius
      
         @return partial derivative based on Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])
      
         --------------------------------------------------------------------------------------- */
      
      static RealOpenMM dL_dx( RealOpenMM r, RealOpenMM x, RealOpenMM S );

      /**---------------------------------------------------------------------------------------
      
         Sgb function
      
         @param t                   r*r*G_i*G_j
      
         @return Sgb (p. 1694 of Labute paper [JCC 29 p. 1693-1698 2008])
      
         --------------------------------------------------------------------------------------- */
      
      static RealOpenMM Sgb( RealOpenMM t );
      
      /**---------------------------------------------------------------------------------------
      
         Get GB/VI energy
      
         @param bornRadii           Born radii
         @param atomCoordinates     atomic coordinates
         @param partialCharges      partial charges
      
         @return energy
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM computeBornEnergy( const RealOpenMM* bornRadii, std::vector<OpenMM::RealVec>& atomCoordinates,
                                    const RealOpenMM* partialCharges );
      
      /**---------------------------------------------------------------------------------------
      
         Get GB/VI forces
      
         @param bornRadii           Born radii
         @param atomCoordinates     atomic coordinates
         @param partialCharges      partial charges
         @param forces              output forces
      
         @return SimTKOpenMMCommon::DefaultReturn;
      
         --------------------------------------------------------------------------------------- */
      
      int computeBornForces( const RealOpenMM* bornRadii, std::vector<OpenMM::RealVec>& atomCoordinates,
                             const RealOpenMM* partialCharges, std::vector<OpenMM::RealVec>& inputForces );
      
      /**---------------------------------------------------------------------------------------
      
         Get volume Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])
      
         @param r                   distance between atoms i & j
         @param R                   atomic radius
         @param S                   scaled atomic radius
      
         @return volume
      
         --------------------------------------------------------------------------------------- */
      
      static double getVolumeD( double r, double R, double S );
      
      /**---------------------------------------------------------------------------------------
      
         Get L (analytical solution for volume integrals) 
      
         @param r                   distance between atoms i & j
         @param R                   atomic radius
         @param S                   scaled atomic radius
      
         @return L value (Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])
      
         --------------------------------------------------------------------------------------- */
      
      static double getLD( double r, double x, double S );

      /**---------------------------------------------------------------------------------------
      
         Get partial derivative of L wrt r
      
         @param r                   distance between atoms i & j
         @param R                   atomic radius
         @param S                   scaled atomic radius
      
         @return partial derivative based on Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])
      
         --------------------------------------------------------------------------------------- */
      
      static double dL_drD( double r, double x, double S );

      /**---------------------------------------------------------------------------------------
      
         Get partial derivative of L wrt x
      
         @param r                   distance between atoms i & j
         @param R                   atomic radius
         @param S                   scaled atomic radius
      
         @return partial derivative based on Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])
      
         --------------------------------------------------------------------------------------- */
      
      static double dL_dxD( double r, double x, double S );

      /**---------------------------------------------------------------------------------------
      
         Compute quintic spline value and associated derviative
      
         @param x                   value to compute spline at
         @param rl                  lower cutoff value
         @param ru                  upper cutoff value
         @param outValue            value of spline at x
         @param outDerivative       value of derivative of spline at x
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int quinticSpline( RealOpenMM x, RealOpenMM rl, RealOpenMM ru, RealOpenMM* outValue, RealOpenMM* outDerivative );
      
      /**---------------------------------------------------------------------------------------
      
         Compute Born radii based on Eq. 3 of Labute paper [JCC 29 p. 1693-1698 2008])
         and quintic splice switching function
      
         @param atomicRadius3       atomic radius cubed
         @param bornSum             Born sum (volume integral)
         @param gbviParameters      Gbvi parameters (parameters used in spline
                                    QuinticLowerLimitFactor & QuinticUpperLimit)
         @param bornRadius          output Born radius
         @param switchDeriviative   output switching function deriviative
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int computeBornRadiiUsingQuinticSpline( RealOpenMM atomicRadius3, RealOpenMM bornSum,
                                              GBVISoftcoreParameters* gbviParameters, 
                                              RealOpenMM* bornRadius, RealOpenMM* switchDeriviative );
      
};

// ---------------------------------------------------------------------------------------

#endif // __CpuGBVISoftcore_H__

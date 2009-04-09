
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

#ifndef __CpuGBVI_H__
#define __CpuGBVI_H__

#include "GBVIParameters.h"
#include "CpuImplicitSolvent.h"

// ---------------------------------------------------------------------------------------

class CpuGBVI : public CpuImplicitSolvent {

   private:

      // GB/VI parameters

      GBVIParameters* _gbviParameters;

      // initialize data members (more than
      // one constructor, so centralize intialization here)

      void _initializeGBVIDataMembers( void );

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         @param implicitSolventParameters    ImplicitSolventParameters reference
      
         @return CpuImplicitSolvent object
      
         --------------------------------------------------------------------------------------- */

       CpuGBVI( ImplicitSolventParameters* gbviParameters );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~CpuGBVI( );

      /**---------------------------------------------------------------------------------------
      
         Return GBVIParameters
      
         @return GBVIParameters
      
         --------------------------------------------------------------------------------------- */

      GBVIParameters* getGBVIParameters( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set ImplicitSolventParameters
      
         @param ImplicitSolventParameters
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */

      int setGBVIParameters( GBVIParameters* gbviParameters );
 
      /**---------------------------------------------------------------------------------------
      
         Get Born radii based on Eq. 3 of Labute paper [JCC 29 p. 1693-1698 2008])

         @param atomCoordinates   atomic coordinates
         @param bornRadii         output array of Born radii
         @param gbviChain         not used
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int computeBornRadii( RealOpenMM** atomCoordinates, RealOpenMM* bornRadii,
                            RealOpenMM* gbviChain = NULL );
      
      /**---------------------------------------------------------------------------------------
      
         Get Born energy and forces (not used)
      
         @param bornRadii         Born radii
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial charges
         @param forces            forces
      
         @return force array
      
         --------------------------------------------------------------------------------------- */
      
      int computeBornEnergyForces( RealOpenMM* bornRadii, RealOpenMM** atomCoordinates,
                                   const RealOpenMM* partialCharges, RealOpenMM** forces );
      
      int computeBornEnergyForcesPrint( RealOpenMM* bornRadii, RealOpenMM** atomCoordinates,
                                        const RealOpenMM* partialCharges, RealOpenMM** forces );
      
      /**---------------------------------------------------------------------------------------
      
         Get state 
      
         title             title (optional)
      
         @return state string
      
         --------------------------------------------------------------------------------------- */
      
      std::string getStateString( const char* title ) const;

      /**---------------------------------------------------------------------------------------
      
         Write Born energy and forces (Simbios)
      
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial atom charges
         @param forces            force array
         @param resultsFileName   output file name
      
         @return SimTKOpenMMCommon::DefaultReturn if file opened; else return SimTKOpenMMCommon::ErrorReturn
      
         --------------------------------------------------------------------------------------- */
          
      int writeBornEnergyForces( RealOpenMM** atomCoordinates,
                                 const RealOpenMM* partialCharges, RealOpenMM** forces,
                                 const std::string& resultsFileName ) const;

      /**---------------------------------------------------------------------------------------
      
         Write  results from first loop
      
         @param atomCoordinates     atomic coordinates
         @param RealOpenMM forces         forces
         @param outputFileName      output file name
      
         @return SimTKOpenMMCommon::DefaultReturn unless
                 file cannot be opened
                 in which case return SimTKOpenMMCommon::ErrorReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int writeForceLoop1( int numberOfAtoms, RealOpenMM** forces, const RealOpenMM* bornForce,
                                  const std::string& outputFileName );
      
      /**---------------------------------------------------------------------------------------
      
         Write results
      
         @param numberOfAtoms       number of atoms
         @param chunkSizes          vector of chunk sizes for realRealOpenMMVector
         @param realRealOpenMMVector      vector of RealOpenMM**
         @param realVector          vector of RealOpenMM*
         @param outputFileName      output file name
      
         @return SimTKOpenMMCommon::DefaultReturn unless
                 file cannot be opened
                 in which case return SimTKOpenMMCommon::ErrorReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int writeForceLoop( int numberOfAtoms, const IntVector& chunkSizes,
                                 const RealOpenMMPtrPtrVector& realRealOpenMMVector, 
                                 const RealOpenMMPtrVector& realVector,
                                 const std::string& outputFileName );
      
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
      
      RealOpenMM computeBornEnergy( const RealOpenMM* bornRadii, RealOpenMM** atomCoordinates,
                                    const RealOpenMM* partialCharges );
      
      /**---------------------------------------------------------------------------------------
      
         Get GB/VI forces
      
         @param bornRadii           Born radii
         @param atomCoordinates     atomic coordinates
         @param partialCharges      partial charges
         @param forces              output forces
      
         @return SimTKOpenMMCommon::DefaultReturn;
      
         --------------------------------------------------------------------------------------- */
      
      int computeBornForces( const RealOpenMM* bornRadii, RealOpenMM** atomCoordinates,
                             const RealOpenMM* partialCharges, RealOpenMM** inputForces );
      
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

};

// ---------------------------------------------------------------------------------------

#endif // __CpuGBVI_H__

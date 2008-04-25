
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

#include <string.h>
#include <sstream>

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "CpuGrycuk.h"
#include <math.h>

/**---------------------------------------------------------------------------------------

   CpuGrycuk constructor

   grycukParameters      grycukParameters object
   
   --------------------------------------------------------------------------------------- */

CpuGrycuk::CpuGrycuk( ImplicitSolventParameters* grycukParameters ) : CpuImplicitSolvent( grycukParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGrycuk::CpuGrycuk";

   // ---------------------------------------------------------------------------------------

   _initializeGrycukDataMembers( );

   _grycukParameters = static_cast<GrycukParameters*> (grycukParameters);

}

/**---------------------------------------------------------------------------------------

   CpuGrycuk destructor

   --------------------------------------------------------------------------------------- */

CpuGrycuk::~CpuGrycuk( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGrycuk::~CpuGrycuk";

   // ---------------------------------------------------------------------------------------

   //if( _grycukParameters != NULL ){
     // delete _grycukParameters;
   //}

}

/**---------------------------------------------------------------------------------------

   Initialize data members

   --------------------------------------------------------------------------------------- */

void CpuGrycuk::_initializeGrycukDataMembers( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGrycuk::initializeDataMembers";

   // ---------------------------------------------------------------------------------------

   _grycukParameters = NULL;
}

/**---------------------------------------------------------------------------------------

   Get GrycukParameters reference

   @return GrycukParameters reference

   --------------------------------------------------------------------------------------- */

GrycukParameters* CpuGrycuk::getGrycukParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGrycuk::getGrycukParameters";

   // ---------------------------------------------------------------------------------------

   return _grycukParameters;
}

/**---------------------------------------------------------------------------------------

   Set GrycukParameters reference

   @param GrycukParameters reference

   @return SimTKOpenMMCommon::DefaultReturn;

   --------------------------------------------------------------------------------------- */

int CpuGrycuk::setGrycukParameters(  GrycukParameters* grycukParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGrycuk::setGrycukParameters";

   // ---------------------------------------------------------------------------------------

   _grycukParameters = grycukParameters;
   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get Born radii based on paper:

      J. Chem Physics 119,9 p. 4817 (2003) (Grycuk paper)

   @param atomCoordinates     atomic coordinates
   @param bornRadii           output array of Born radii

   @return array of Born radii

   --------------------------------------------------------------------------------------- */

int CpuGrycuk::computeBornRadii( RealOpenMM** atomCoordinates, RealOpenMM* bornRadii, RealOpenMM* grycukChain ){

   // ---------------------------------------------------------------------------------------

   static const RealOpenMM zero    = (RealOpenMM) 0.0;
   static const RealOpenMM one     = (RealOpenMM) 1.0;
   static const RealOpenMM two     = (RealOpenMM) 2.0;
   static const RealOpenMM three   = (RealOpenMM) 3.0;
   static const RealOpenMM four    = (RealOpenMM) 4.0;
   static const RealOpenMM six     = (RealOpenMM) 6.0;
   static const RealOpenMM eight   = (RealOpenMM) 8.0;
   static const RealOpenMM half    = (RealOpenMM) 0.5;
   static const RealOpenMM third   = (RealOpenMM) 1.0/ (RealOpenMM) 3.0;
   static const RealOpenMM fourth  = (RealOpenMM) 0.25;

   static const RealOpenMM pi4x3   = (RealOpenMM) (PI_M*4.0/3.0);
   static const RealOpenMM pi4     = (RealOpenMM) (PI_M/4.0);
   static const RealOpenMM pi12    = (RealOpenMM) (PI_M/12.0);

   // static const char* methodName = "\nCpuGrycuk::computeBornRadii";

   // ---------------------------------------------------------------------------------------

   GrycukParameters* grycukParameters    = getGrycukParameters();

   int numberOfAtoms                     = grycukParameters->getNumberOfAtoms();
   RealOpenMM* atomicRadii               = grycukParameters->getAtomicRadii();
   const RealOpenMM* scaledRadiusFactor  = grycukParameters->getScaledRadiusFactors();

   // ---------------------------------------------------------------------------------------

   // calculate Born radii

FILE* logFile = SimTKOpenMMLog::getSimTKOpenMMLogFile( );

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
     
      RealOpenMM radiusI         = atomicRadii[atomI];
      RealOpenMM sum             = POW( radiusI, -three );
                 sum            *= pi4x3;

      RealOpenMM rI3             = POW( radiusI, -three );

      // HCT code

      for( int atomJ = 0; atomJ < numberOfAtoms; atomJ++ ){

         if( atomJ != atomI ){

            RealOpenMM deltaX          = atomCoordinates[atomJ][0] - atomCoordinates[atomI][0];
            RealOpenMM deltaY          = atomCoordinates[atomJ][1] - atomCoordinates[atomI][1];
            RealOpenMM deltaZ          = atomCoordinates[atomJ][2] - atomCoordinates[atomI][2];
 
            RealOpenMM r2              = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
            RealOpenMM r               = SQRT( r2 );
            RealOpenMM scaledRadiusJ   = atomicRadii[atomJ]*scaledRadiusFactor[atomJ];
            RealOpenMM scaledRadiusJ2  = scaledRadiusJ*scaledRadiusJ;
//            RealOpenMM rScaledRadiusJ  = r + scaledRadiusJ;

//          r2 = xr**2 + yr**2 + zr**2
//          r = sqrt(r2)    
//          sk = rk * shct(k)
//          sk2 = sk * sk
//          if (ri+r .lt. sk) then
//             lik = ri
//             uik = sk - r 
//             sum = sum + pi43*(1.0d0/uik**3-1.0d0/lik**3)
//          end if
//          uik = r + sk
//          if (ri+r .lt. sk) then
//             lik = sk - r 
//          else if (r .lt. ri+sk) then
//             lik = ri
//          else
//             lik = r - sk
//          end if
//          l2 = lik * lik 
//          l4 = l2 * l2
//          lr = lik * r 
//          l4r = l4 * r 
//          u2 = uik * uik 
//          u4 = u2 * u2
//          ur = uik * r 
//          u4r = u4 * r 
//          term = (3.0d0*(r2-sk2)+6.0d0*u2-8.0d0*ur)/u4r - (3.0d0*(r2-sk2)+6.0d0*l2-8.0d0*lr)/l4r
//          sum = sum - pi*term/12.0d0


            RealOpenMM rRadiusI    = radiusI + r;

            if( rRadiusI < scaledRadiusJ ){
               RealOpenMM u_ij3    = POW( (scaledRadiusJ - r), -three );
               sum                += pi4x3*(u_ij3 - rI3);
            }

            RealOpenMM l_ij;
            RealOpenMM u_ij = r + scaledRadiusJ;
            if( rRadiusI < scaledRadiusJ ){
               l_ij = scaledRadiusJ - r;
            } else if( r < (radiusI+scaledRadiusJ) ){
               l_ij = radiusI;
            } else {
               l_ij = r - scaledRadiusJ;
            }

            RealOpenMM l_ij2    = l_ij*l_ij;
            RealOpenMM l_ij4    = l_ij2*l_ij2;
            RealOpenMM l_r      = l_ij*r;
            RealOpenMM l4r      = l_ij4*r;

            RealOpenMM u_ij2    = u_ij*u_ij;
            RealOpenMM u_ij4    = u_ij2*u_ij2;
            RealOpenMM u_r      = u_ij*r;
            RealOpenMM u4r      = u_ij4*r;

            RealOpenMM r2Diff3  = three*(r2 - scaledRadiusJ2);
            RealOpenMM term     = ( r2Diff3 + six*u_ij2 - eight*u_r)/u4r - ( r2Diff3 + six*l_ij2 - eight*l_r)/l4r;
            sum                 = sum - pi12*term;


if( atomI >= -217 && atomI <= -219 ){
   (void) fprintf( logFile, "\nRR %d %d r=%.4f rads[%.6f %.6f] scaledJ[%.3f %.3f] sum=%12.6e l/u=[%12.6e %12.6e] term=%12.6e",
                   atomI, atomJ, r, radiusI, atomicRadii[atomJ], scaledRadiusJ, scaledRadiusFactor[atomJ], sum, l_ij, u_ij, term );
}

         }
      }
 
      sum                  /= pi4x3;
      if( sum <= 0.0 ){
         bornRadii[atomI]   = 500.0;
      } else {
         bornRadii[atomI]   = POW( sum, -third );
      }

if( atomI >= 0 ){
   (void) fprintf( logFile, "\nRRQ %d sum=%12.6e radI=%.5f born=%12.6e",
                   atomI, sum, radiusI, bornRadii[atomI] );
}

   }

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Get Grycuk Born energy and forces

   @param bornRadii           Born radii -- optional; if NULL, then GrycukParameters 
                              entry is used
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces

   @return SimTKOpenMMCommon::DefaultReturn;

   The array bornRadii is also updated and the grycukEnergy

   --------------------------------------------------------------------------------------- */

int CpuGrycuk::computeBornEnergyForces( RealOpenMM* bornRadii, RealOpenMM** atomCoordinates,
                                        const RealOpenMM* partialCharges, RealOpenMM** forces ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGrycuk::computeBornEnergyForces";

   static const RealOpenMM zero       = (RealOpenMM) 0.0;
   static const RealOpenMM one        = (RealOpenMM) 1.0;
   static const RealOpenMM two        = (RealOpenMM) 2.0;
   static const RealOpenMM three      = (RealOpenMM) 3.0;
   static const RealOpenMM four       = (RealOpenMM) 4.0;
   static const RealOpenMM five       = (RealOpenMM) 5.0;
   static const RealOpenMM six        = (RealOpenMM) 6.0;
   static const RealOpenMM eight      = (RealOpenMM) 8.0;

   static const RealOpenMM half       = (RealOpenMM) 0.5;
   static const RealOpenMM third      = (RealOpenMM) 1.0/ (RealOpenMM) 3.0;
   static const RealOpenMM fourth     = (RealOpenMM) 0.25;
   static const RealOpenMM eighth     = (RealOpenMM) 0.125;
   static const RealOpenMM nine       = (RealOpenMM) 9.0;
   static const RealOpenMM seventeen  = (RealOpenMM) 9.0;
   static const RealOpenMM pi12       = (RealOpenMM) (  (RealOpenMM) PI_M/ (RealOpenMM) 12.0);

   // ---------------------------------------------------------------------------------------

   const GrycukParameters* grycukParameters = getGrycukParameters();
   const int numberOfAtoms            = grycukParameters->getNumberOfAtoms();

   if( bornRadii == NULL ){
      bornRadii   = getBornRadii();
   }

   // ---------------------------------------------------------------------------------------

   // constants

   const RealOpenMM preFactor           = grycukParameters->getPreFactor();
   const RealOpenMM electricConstant    = grycukParameters->getElectricConstant();

FILE* logFile = SimTKOpenMMLog::getSimTKOpenMMLogFile( );
(void) fprintf( logFile, "electricConstant=%.3f preFact=%.3f", electricConstant, preFactor );

   // ---------------------------------------------------------------------------------------

   // set energy/forces to zero

   RealOpenMM grycukEnergy              = zero;
   const unsigned int arraySzInBytes    = sizeof( RealOpenMM )*numberOfAtoms;

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      memset( forces[ii], 0, 3*sizeof( RealOpenMM ) );
   }

   RealOpenMM* bornForces = getBornForce();
   memset( bornForces, 0, arraySzInBytes );

   // ---------------------------------------------------------------------------------------

   // N*( 8 + pow) ACE
   // compute the nonpolar solvation via ACE approximation
    
   if( includeAceApproximation() ){
      computeAceNonPolarForce( grycukParameters, bornRadii, &grycukEnergy, bornForces );
   }

   // ---------------------------------------------------------------------------------------

   // first main loop

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      RealOpenMM partialChargeI = preFactor*partialCharges[atomI];
      for( int atomJ = atomI; atomJ < numberOfAtoms; atomJ++ ){

         // 3 FLOP

         RealOpenMM deltaX             = atomCoordinates[atomJ][0] - atomCoordinates[atomI][0];
         RealOpenMM deltaY             = atomCoordinates[atomJ][1] - atomCoordinates[atomI][1];
         RealOpenMM deltaZ             = atomCoordinates[atomJ][2] - atomCoordinates[atomI][2];
 
         // 5 FLOP

         RealOpenMM r2                 = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;

         // 3 FLOP

         RealOpenMM alpha2_ij          = bornRadii[atomI]*bornRadii[atomJ];
         RealOpenMM D_ij               = r2/(four*alpha2_ij);

         // exp + 2 + sqrt FLOP 

         RealOpenMM expTerm            = EXP( -D_ij );
         RealOpenMM denominator2       = r2 + alpha2_ij*expTerm; 
         RealOpenMM denominator        = SQRT( denominator2 ); 
         
         // 6 FLOP

         RealOpenMM Gpol               = (partialChargeI*partialCharges[atomJ])/denominator; 
         RealOpenMM dGpol_dr           = -Gpol*( one - fourth*expTerm )/denominator2;  

         // 5 FLOP

         RealOpenMM dGpol_dalpha2_ij   = -half*Gpol*expTerm*( one + D_ij )/denominator2;

         // 11 FLOP

         if( atomI != atomJ ){

             bornForces[atomJ] += dGpol_dalpha2_ij*bornRadii[atomI];

             deltaX            *= dGpol_dr;
             deltaY            *= dGpol_dr;
             deltaZ            *= dGpol_dr;

             forces[atomI][0]  -= deltaX;
             forces[atomI][1]  -= deltaY;
             forces[atomI][2]  -= deltaZ;

             forces[atomJ][0]  += deltaX;
             forces[atomJ][1]  += deltaY;
             forces[atomJ][2]  += deltaZ;

         } else {
            Gpol *= half;
         }

         // 3 FLOP

         grycukEnergy         += Gpol;
         bornForces[atomI]    += dGpol_dalpha2_ij*bornRadii[atomJ];

if( atomI == 0 ){
   (void) fprintf( logFile, "\nVV %d %d r2=%.4f dGpol_dalpha2_ij=%.4e bR=%12.5e bF=%.3e",
                   atomI, atomJ, r2, dGpol_dalpha2_ij, bornRadii[atomJ], bornForces[atomI] );
}
      }
   }

   //grycukEnergy *= getEnergyConversionFactor();

   // ---------------------------------------------------------------------------------------

   // second main loop

   // initialize Born radii & GrycukChain temp arrays -- contain values
   // used in next iteration

   RealOpenMM* bornRadiiTemp             = getBornRadiiTemp();
   memset( bornRadiiTemp, 0, arraySzInBytes );

   const RealOpenMM* atomicRadii         = grycukParameters->getAtomicRadii();
   const RealOpenMM* scaledRadiusFactor  = grycukParameters->getScaledRadiusFactors();

   const RealOpenMM fourThirds           = four*third;
   const RealOpenMM pi43                 = fourThirds*PI_M;
   const RealOpenMM pi4x3                = (RealOpenMM) (PI_M*4.0/3.0);
   const RealOpenMM factor               = -( POW( PI_M, third)* POW( six, two*third) )/nine;
   const RealOpenMM fourthPi             = fourth*PI_M;

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      // radius

      RealOpenMM radiusI        = atomicRadii[atomI];
      RealOpenMM radiusI2       = radiusI*radiusI;
      RealOpenMM radiusI3       = radiusI2*radiusI;
      RealOpenMM radiusI4       = radiusI2*radiusI2;
 
      RealOpenMM term           = POW( bornRadii[atomI], -three );
                 term          *= pi43;
                 term           = POW( term, -fourThirds );
                 term          *= factor;

      // used to compute Born radius for next iteration

      RealOpenMM bornSum        = POW( radiusI, -three );
                 bornSum       *= pi4x3;

      for( int atomJ = 0; atomJ < numberOfAtoms; atomJ++ ){

         if( atomJ != atomI ){

            RealOpenMM deltaX             = atomCoordinates[atomJ][0] - atomCoordinates[atomI][0];
            RealOpenMM deltaY             = atomCoordinates[atomJ][1] - atomCoordinates[atomI][1];
            RealOpenMM deltaZ             = atomCoordinates[atomJ][2] - atomCoordinates[atomI][2];
    
            RealOpenMM r2                 = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
            RealOpenMM r                  = SQRT( r2 );
 
            RealOpenMM scaledRadiusJ      = atomicRadii[atomJ]*scaledRadiusFactor[atomJ];
            RealOpenMM scaledRadiusJ2     = scaledRadiusJ*scaledRadiusJ;
            RealOpenMM rScaledRadiusJ     = r + scaledRadiusJ;

            RealOpenMM rRadiusI           = radiusI + r;

            RealOpenMM rRadiusJ4          = r - scaledRadiusJ;
                       rRadiusJ4          = rRadiusJ4*rRadiusJ4;
                       rRadiusJ4          = rRadiusJ4*rRadiusJ4;
            RealOpenMM rRadiusJ4_2        = one/(rRadiusJ4*r2);

            RealOpenMM de;
            RealOpenMM l_ij;
            RealOpenMM u_ij               =  r + scaledRadiusJ;
            if( rRadiusI < scaledRadiusJ ){

               RealOpenMM rRadiusJ        = scaledRadiusJ - r;

               de                         = -four*PI_M/rRadiusJ4;
               de                        += fourthPi*(scaledRadiusJ2 - four*scaledRadiusJ*r + seventeen*r2)*rRadiusJ4_2;

               RealOpenMM u_ij3           = POW( rRadiusJ, -three );
               bornSum                   += pi4x3*(u_ij3 - radiusI3);
               l_ij                       = rRadiusJ;
            } else if( r < (radiusI + scaledRadiusJ) ){

               de                         = fourthPi*( two*radiusI*radiusI - scaledRadiusJ2 - r2)/(r2*radiusI4);
               l_ij                       = radiusI;

            } else {
               de                         = fourthPi*( scaledRadiusJ2 - four*scaledRadiusJ*r + r2)*rRadiusJ4_2;
               l_ij                       = r - scaledRadiusJ;
            }

            // born R

            RealOpenMM l_ij2              = l_ij*l_ij;
            RealOpenMM l_ij4              = l_ij2*l_ij2;
            RealOpenMM l_r                = l_ij*r;
            RealOpenMM l4r                = l_ij4*r;

            RealOpenMM u_ij2              = u_ij*u_ij;
            RealOpenMM u_ij4              = u_ij2*u_ij2;
            RealOpenMM u_r                = u_ij*r;
            RealOpenMM u4r                = u_ij4*r;

            RealOpenMM r2Diff3            = three*(r2 - scaledRadiusJ2);
            RealOpenMM bornTerm           = ( r2Diff3 + six*u_ij2 - eight*u_r)/u4r - ( r2Diff3 + six*l_ij2 - eight*l_r)/l4r;
            bornSum                       = bornSum - pi12*bornTerm;

            // force

            u_ij4                         = POW( rScaledRadiusJ, four );
            de                           -= fourthPi*(scaledRadiusJ2 + four*scaledRadiusJ*r + r2)/(r2*u_ij4);

            de                           *= (term*bornForces[atomI])/r;

            deltaX                       *= de;
            deltaY                       *= de;
            deltaZ                       *= de;

            forces[atomI][0]             += deltaX;
            forces[atomI][1]             += deltaY;
            forces[atomI][2]             += deltaZ;
  
            forces[atomJ][0]             -= deltaX;
            forces[atomJ][1]             -= deltaY;
            forces[atomJ][2]             -= deltaZ;
 
if( atomI == 0 ){
   (void) fprintf( logFile, "\nWW %d %d r=%.4f rads[%.6f %.6f] scaledJ[%.6f %.3f] de=%12.6e term=%12.6e bF=%12.5e",
                   atomI, atomJ, r, radiusI, atomicRadii[atomJ], scaledRadiusJ, scaledRadiusFactor[atomJ], de, term, bornForces[atomI] );
}

         }
      }

      bornSum                            /= pi4x3;
      if( bornSum <= 0.0 ){
         bornRadiiTemp[atomI]             = 500.0;
      } else {
         bornRadiiTemp[atomI]             = POW( bornSum, -third );
      }

   }

   setEnergy( grycukEnergy );

   // copy new Born radii and grycukChain values into permanent array

   memcpy( bornRadii, bornRadiiTemp, arraySzInBytes );

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------
      
   Get string w/ state 
   
   @param title               title (optional)
      
   @return string containing state
      
   --------------------------------------------------------------------------------------- */

std::string CpuGrycuk::getStateString( const char* title ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getStateString";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   message << CpuImplicitSolvent::getStateString( title );

   return message.str();
}

/**---------------------------------------------------------------------------------------

   Write Born energy and forces (Simbios)

   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces
   @param resultsFileName     output file name

   @return SimTKOpenMMCommon::DefaultReturn unless
           file cannot be opened
           in which case return SimTKOpenMMCommon::ErrorReturn

   --------------------------------------------------------------------------------------- */

int CpuGrycuk::writeBornEnergyForces( RealOpenMM** atomCoordinates,
                                      const RealOpenMM* partialCharges, RealOpenMM** forces,
                                      const std::string& resultsFileName ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nCpuGrycuk::writeBornEnergyForces";

   // ---------------------------------------------------------------------------------------

   ImplicitSolventParameters* implicitSolventParameters = getImplicitSolventParameters();
   const GrycukParameters* grycukParameters             = static_cast<const GrycukParameters*>(implicitSolventParameters);
   

   int numberOfAtoms                                    = grycukParameters->getNumberOfAtoms();
   const RealOpenMM* atomicRadii                        = grycukParameters->getAtomicRadii();
   const RealOpenMM* bornRadii                          = getBornRadiiConst();
   const RealOpenMM* scaledRadiusFactor                 = grycukParameters->getScaledRadiusFactors();

   // ---------------------------------------------------------------------------------------

   // open file -- return if unsuccessful

   FILE* implicitSolventResultsFile = NULL;
#ifdef WIN32
   fopen_s( &implicitSolventResultsFile, resultsFileName.c_str(), "w" );
#else
   implicitSolventResultsFile = fopen( resultsFileName.c_str(), "w" );
#endif

   // diganostics

   std::stringstream message;
   message << methodName;
   if( implicitSolventResultsFile != NULL ){
      std::stringstream message;
      message << methodName;
      message << " Opened file=<" << resultsFileName << ">.";
      SimTKOpenMMLog::printMessage( message );
   } else {
      std::stringstream message;
      message << methodName;
      message << "  could not open file=<" << resultsFileName << "> -- abort output.";
      SimTKOpenMMLog::printMessage( message );
      return SimTKOpenMMCommon::ErrorReturn;
   }

   // header

   (void) fprintf( implicitSolventResultsFile, "# %d atoms format: coords(3) bornRadii(input) q atomicRadii forces\n", numberOfAtoms );

   RealOpenMM forceConversion  = (RealOpenMM) 1.0;
   RealOpenMM lengthConversion = (RealOpenMM) 1.0;

   // output

   if( forces != NULL && atomCoordinates != NULL && partialCharges != NULL && atomicRadii != NULL ){
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
            (void) fprintf( implicitSolventResultsFile, "%.7e %.7e %.7e %.7e %.5f %.5f %.3f %.7e %.7e %.7e\n",
                            lengthConversion*atomCoordinates[ii][0],
                            lengthConversion*atomCoordinates[ii][1], 
                            lengthConversion*atomCoordinates[ii][2],
                           (bornRadii != NULL ? lengthConversion*bornRadii[ii] : 0.0),
                            partialCharges[ii], lengthConversion*atomicRadii[ii],
                            scaledRadiusFactor[ii],
                            forceConversion*forces[ii][0],
                            forceConversion*forces[ii][1],
                            forceConversion*forces[ii][2]
                          );
      }
   }
   (void) fclose( implicitSolventResultsFile );

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Write  results from first loop

   @param numberOfAtoms       number of atoms
   @param forces              forces
   @param bornForce           Born force prefactor
   @param outputFileName      output file name

   @return SimTKOpenMMCommon::DefaultReturn unless
           file cannot be opened
           in which case return SimTKOpenMMCommon::ErrorReturn

   --------------------------------------------------------------------------------------- */

int CpuGrycuk::writeForceLoop1( int numberOfAtoms, RealOpenMM** forces, const RealOpenMM* bornForce,
                                const std::string& outputFileName ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nCpuGrycuk::writeForceLoop1";

   // ---------------------------------------------------------------------------------------

   int chunkSize;
   if( bornForce ){
      chunkSize = 3;
   } else {
      chunkSize = 4;
   }

   StringVector lineVector;
   std::stringstream header;
   lineVector.push_back( "# bornF F" );
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      std::stringstream line;
      line << (atomI+1) << " ";
      SimTKOpenMMUtilities::formatRealStringStream( line, forces[atomI], chunkSize );
      if( bornForce ){
         line << " " << bornForce[atomI];
      }
      lineVector.push_back( line.str() );
   }
   return SimTKOpenMMUtilities::writeFile( lineVector, outputFileName );

}

/**---------------------------------------------------------------------------------------

   Write results

   @param numberOfAtoms        number of atoms
   @param chunkSizes           vector of chunk sizes for realRealOpenMMVector
   @param realRealOpenMMVector vector of RealOpenMM**
   @param realVector           vector of RealOpenMM*
   @param outputFileName       output file name

   @return SimTKOpenMMCommon::DefaultReturn unless
           file cannot be opened
           in which case return SimTKOpenMMCommon::ErrorReturn

   --------------------------------------------------------------------------------------- */

int CpuGrycuk::writeForceLoop( int numberOfAtoms, const IntVector& chunkSizes,
                               const RealOpenMMPtrPtrVector& realRealOpenMMVector, 
                               const RealOpenMMPtrVector& realVector,
                               const std::string& outputFileName ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nCpuGrycuk::writeForceLoop";

   static const int maxChunks = 10;
   int chunks[maxChunks];

   // ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < (int) chunkSizes.size(); ii++ ){
      chunks[ii] = chunkSizes[ii];
   }
   for( int ii = (int) chunkSizes.size(); ii < maxChunks; ii++ ){
      chunks[ii] = 3;
   }

   StringVector lineVector;
   std::stringstream header;
   // lineVector.push_back( "# " );

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      std::stringstream line;
      line << (atomI+1) << " ";

      int index = 0;
      for( RealOpenMMPtrPtrVectorCI ii = realRealOpenMMVector.begin(); ii != realRealOpenMMVector.end(); ii++ ){
         RealOpenMM** forces = *ii;
         SimTKOpenMMUtilities::formatRealStringStream( line, forces[atomI], chunks[index++] );
         line << " ";
      }

      for( RealOpenMMPtrVectorCI ii = realVector.begin(); ii != realVector.end(); ii++ ){
         RealOpenMM* array = *ii;
         line << array[atomI] << " ";
      }

      lineVector.push_back( line.str() );
   }
   return SimTKOpenMMUtilities::writeFile( lineVector, outputFileName );

}

/**---------------------------------------------------------------------------------------

   Get Grycuk Born energy and forces -- used debugging

   @param bornRadii           Born radii -- optional; if NULL, then GrycukParameters 
                              entry is used
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces

   @return SimTKOpenMMCommon::DefaultReturn;

   The array bornRadii is also updated and the grycukEnergy

   --------------------------------------------------------------------------------------- */

int CpuGrycuk::computeBornEnergyForcesPrint( RealOpenMM* bornRadii, RealOpenMM** atomCoordinates,
                                             const RealOpenMM* partialCharges, RealOpenMM** forces ){
 
   // ---------------------------------------------------------------------------------------


   return SimTKOpenMMCommon::DefaultReturn;

}

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Chris Sweet                                                       *
 * Contributors: Christopher Bruns, Pande Group                               *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include <string.h>
#include <sstream>

#include "SimTKUtilities/SimTKOpenMMCommon.h"
#include "SimTKUtilities/SimTKOpenMMLog.h"
#include "SimTKUtilities/SimTKOpenMMUtilities.h"
#include "ReferenceNMLDynamics.h"

/**---------------------------------------------------------------------------------------

   ReferenceNMLDynamics constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         delta t for dynamics
   @param tau            viscosity(?)
   @param temperature    temperature

   --------------------------------------------------------------------------------------- */

ReferenceNMLDynamics::ReferenceNMLDynamics( int numberOfAtoms,
                                              RealOpenMM deltaT, RealOpenMM tau,
                                              RealOpenMM temperature,
                                              RealOpenMM* projectionVectors, 
                                              unsigned int numProjectionVectors, 
                                              RealOpenMM minimumLimit, RealOpenMM maxEig  ) : 
            ReferenceDynamics( numberOfAtoms, deltaT, temperature ), _tau( tau ),
            _projectionVectors(projectionVectors), _numProjectionVectors(numProjectionVectors), 
            _minimumLimit(minimumLimit), _maxEig(maxEig)  {

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceNMLDynamics::ReferenceNMLDynamics";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;

   // ---------------------------------------------------------------------------------------

   // insure tau is not zero -- if it is print warning message

   if( _tau == zero ){

      std::stringstream message;
      message << methodName;
      message << " input tau value=" << tau << " is invalid -- setting to 1.";
      SimTKOpenMMLog::printError( message );

      _tau = one;
     
   }
   _setFixedParameters( );

   allocate2DArrays( numberOfAtoms, 3, Max2DArrays );
   allocate1DArrays( numberOfAtoms, Max1DArrays );
   
}

/**---------------------------------------------------------------------------------------

   ReferenceNMLDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceNMLDynamics::~ReferenceNMLDynamics( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceNMLDynamics::~ReferenceNMLDynamics";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Set fixed parameters

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceNMLDynamics::_setFixedParameters( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceNMLDynamics::_setFixedParameters";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;
   static const RealOpenMM two        =  2.0;
   static const RealOpenMM half       =  0.5;
    static const RealOpenMM eighty       =  80.0;

   // ---------------------------------------------------------------------------------------

    //dt, myGamma, fdt, vdt, ndt, sqrtFCoverM from Langevin Leapfrog
    
    _fixedParameters[DT]         = getTimeStep();
    _fixedParameters[GAMMA]      = (_tau == zero ? eighty : one/_tau);
    _fixedParameters[FDT]        = ( one - EXP( -half * _fixedParameters[GAMMA] * _fixedParameters[DT] ) ) / _fixedParameters[GAMMA];
    _fixedParameters[VDT]        = EXP(-half * _fixedParameters[GAMMA] * _fixedParameters[DT]);
    _fixedParameters[NDT]        = SQRT( ( one - EXP( -_fixedParameters[GAMMA] * _fixedParameters[DT] ) ) / (two * _fixedParameters[GAMMA]) ); 
    _fixedParameters[SQRTFCOVERM]        = SQRT( two * (RealOpenMM) BOLTZ * getTemperature() * _fixedParameters[GAMMA] );
    
   return 0; // ReferenceDynamics::DefaultReturn;

};

/**---------------------------------------------------------------------------------------

   Get tau

   @return tau

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceNMLDynamics::getTau( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceNMLDynamics::getTau";

   // ---------------------------------------------------------------------------------------

   return _tau;
}

/**---------------------------------------------------------------------------------------

   Get array of fixed parameters indexed by 'FixedParameters' enums

   @return array

   --------------------------------------------------------------------------------------- */
   
const RealOpenMM* ReferenceNMLDynamics::getFixedParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceNMLDynamics::getFixedParameters";

   // ---------------------------------------------------------------------------------------

   return _fixedParameters;
}

/**---------------------------------------------------------------------------------------

   Print parameters

   @param message             message

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceNMLDynamics::printParameters( std::stringstream& message ) const {

   // ---------------------------------------------------------------------------------------
    //dt, myGamma, fdt, vdt, ndt, sqrtFCoverM from Langevin Leapfrog

   static const char* methodName  = "\nReferenceNMLDynamics::printParameters";
   static const char* parameterNames[MaxFixedParameters] = { "dt", "myGamma", "fdt", "ndt", "sqrtFCoverM" };

   // ---------------------------------------------------------------------------------------

   // print parameters

   ReferenceDynamics::printParameters( message );
   message << " tau=" << getTau();
   message << " T=" << getTemperature();
   int cut = 3;
   for( int ii = 0; ii < MaxFixedParameters; ii++ ){
      message << " " << parameterNames[ii] << "=" << _fixedParameters[ii];
      if( cut++ > 5 ){
         cut = 0;
         message << std::endl;
      }
   }

   return 0; // ReferenceDynamics::DefaultReturn;

}

/**---------------------------------------------------------------------------------------
 
 Drift update
 
 @param numberOfAtoms       number of atoms
 @param atomCoordinates     atom coordinates
 @param velocities          velocities
 @param forces              forces
 @param inverseMasses       inverse atom masses
 @param xPrime              xPrime
 @param oldVelocities       previous velocities
 @param xVector             xVector
 @param vVector             vVector
 
 @return ReferenceDynamics::DefaultReturn
 
 --------------------------------------------------------------------------------------- */

void ReferenceNMLDynamics::halfKick( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                      RealOpenMM** velocities,
                                      RealOpenMM** forces, RealOpenMM* masses, 
                                      RealOpenMM* inverseMasses, RealOpenMM** xVector ){
    
    //dt, myGamma, fdt, vdt, ndt, sqrtFCoverM from Langevin Leapfrog
    static const RealOpenMM dt = getTimeStep(); // in ps
    static const RealOpenMM myGamma = (_tau == 0.0 ? 80.0 : 1.0/_tau);
    static const RealOpenMM fdt = ( 1.0 - EXP( -0.5 * myGamma * dt ) ) / myGamma;
    static const RealOpenMM vdt = EXP(-0.5 * myGamma * dt);
    static const RealOpenMM ndt = SQRT( ( 1.0 - EXP( -myGamma * dt ) ) / (2.0 * myGamma) ); 
    static const RealOpenMM sqrtFCoverM = SQRT( 2.0 * (RealOpenMM) BOLTZ * getTemperature() * myGamma );
    
    //project forces into sub space. TODO put this somewhere callable from NMLIntegrator, else projecting twice.
    subspaceProjection(forces, forces, numberOfAtoms, masses, inverseMasses, 0);

    for( int ii = 0; ii < numberOfAtoms; ii++ ){
        //RealOpenMM sqrtInverseMass = SQRT( inverseMasses[ii] )*fixedParameters[X];
        for( int jj = 0; jj < 3; jj++ ){
            xVector[ii][jj] = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber() * SQRT( masses[ii] );
        }
    }
    
    subspaceProjection(xVector, xVector, numberOfAtoms, masses, inverseMasses, 0);
    
    for( int ii = 0; ii < numberOfAtoms; ii++ ){
        //RealOpenMM sqrtInverseMass = SQRT( inverseMasses[ii] )*fixedParameters[X];
        for( int jj = 0; jj < 3; jj++ ){
            xVector[ii][jj] *= inverseMasses[ii] ;
        }
    }
    
    
    
    for( int i = 0; i < numberOfAtoms; i++ ) {
        // semi-update velocities
        for( int j = 0; j < 3; j++ ) {
            velocities[i][j] = velocities[i][j] * vdt
                                    + forces[i][j] * fdt * inverseMasses[i]
                                        + xVector[i][j] * sqrtFCoverM * ndt;
                                        //+ SimTKOpenMMUtilities::getNormallyDistributedRandomNumber() * sqrtFCoverM * ndt;
        }
    }
    
    //project forces into complement space, put in xPrime
    subspaceProjection(velocities, velocities, numberOfAtoms, masses, inverseMasses, 1);    //in subspace but position/velocity projection

}

void ReferenceNMLDynamics::drift( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                   RealOpenMM** velocities,
                                   RealOpenMM** forces, RealOpenMM* masses, RealOpenMM* inverseMasses,
                                   RealOpenMM** xPrime ){
    
    static const RealOpenMM deltaT = getTimeStep();
    
    for( int i = 0; i < numberOfAtoms; i++ ){
        for( int j = 0; j < 3; j++ ) {
        
            xPrime[i][j] = atomCoordinates[i][j] + deltaT * velocities[i][j];
            
        }
    }
    

}
    
/**---------------------------------------------------------------------------------------

   Update -- driver routine for performing stochastic dynamics update of coordinates
   and velocities

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceNMLDynamics::update( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                          RealOpenMM** velocities,
                                          RealOpenMM** forces, RealOpenMM* masses, 
                                          const RealOpenMM currentPE, const int stepType ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceNMLDynamics::update";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;

   static int debug                   =  0;

   // ---------------------------------------------------------------------------------------

   // get work arrays

   RealOpenMM** xPrime          = get2DArrayAtIndex( xPrime2D );
   RealOpenMM** oldVelocities   = get2DArrayAtIndex( OldV );
   RealOpenMM** xVector         = get2DArrayAtIndex( X2D );
   RealOpenMM** vVector         = get2DArrayAtIndex( V2D );

   RealOpenMM* inverseMasses    = get1DArrayAtIndex( InverseMasses );

   // first-time-through initialization

   if( getTimeStep() == 0 ){

      std::stringstream message;
      message << methodName;
      int errors = 0;

      // invert masses

      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         if( masses[ii] <= zero ){
            message << "mass at atom index=" << ii << " (" << masses[ii] << ") is <= 0" << std::endl;
            errors++;
         } else {
            inverseMasses[ii] = one/masses[ii];
         }
      }

      // set xVector 

      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         //RealOpenMM sqrtInverseMass = SQRT( inverseMasses[ii] )*fixedParameters[X];
         for( int jj = 0; jj < 3; jj++ ){
            xVector[ii][jj] = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber() * SQRT( masses[ii] );
         }
      }

      subspaceProjection(xVector, xVector, numberOfAtoms, masses, inverseMasses, 0);
       
       for( int ii = 0; ii < numberOfAtoms; ii++ ){
           //RealOpenMM sqrtInverseMass = SQRT( inverseMasses[ii] )*fixedParameters[X];
           for( int jj = 0; jj < 3; jj++ ){
               xVector[ii][jj] *= inverseMasses[ii] ;
           }
       }
       
      
      // exit if errors

      if( errors ){
         SimTKOpenMMLog::printError( message );
      }
   }

    std::stringstream messagea;
    messagea << "Tau " << _tau << "\n";
    //messagea << "Eigs " << _projectionVectors[0] << "," << _projectionVectors[1] << "," << _projectionVectors[2] << "," << _projectionVectors[3] << ","  <<"\n";
    //messagea << "Num vect " << _numProjectionVectors << "\n";
    SimTKOpenMMLog::printMessage( messagea );
    
    switch(stepType){
        case 1: 
            halfKick( numberOfAtoms, atomCoordinates,
                           velocities,
                           forces, masses, inverseMasses, xVector );
            
            drift( numberOfAtoms, atomCoordinates,
                     velocities,
                     forces, masses, inverseMasses,
                     xPrime );
            
            break;
            
        case 2: 
            halfKick( numberOfAtoms, atomCoordinates,
                                  velocities,
                                  forces, masses, inverseMasses, xVector );
            break;
            
        case 3:{
            //save current PE in case quadratic required
            lastPE = currentPE;
            //project forces into complement space, put in xPrime
            subspaceProjection(forces, xPrime, numberOfAtoms, masses, inverseMasses, 2);

            //just equal to minus dot product of 'proposed position move' and forces (=-\nabla PE)
            lastSlope = 0.0;
            for( int ii = 0; ii < numberOfAtoms; ii++ ){
                for( int jj = 0; jj < 3; jj++ ){
                    lastSlope -= xPrime[ii][jj] * forces[ii][jj] * inverseMasses[ii]; //posTemp[k].dot((*myForcesP)[k]);
                }
            }
            
            std::stringstream message1;
            message1 << "Last slope " << lastSlope  <<"\n";
            
            SimTKOpenMMLog::printMessage( message1 );
            
            // copy xPrime -> atomCoordinates
            for( int ii = 0; ii < numberOfAtoms; ii++ ){
                const RealOpenMM factor = inverseMasses[ii] / _maxEig;
    
                xPrime[ii][0] = atomCoordinates[ii][0] + factor * xPrime[ii][0];
                xPrime[ii][1] = atomCoordinates[ii][1] + factor * xPrime[ii][1];
                xPrime[ii][2] = atomCoordinates[ii][2] + factor * xPrime[ii][2];
            }
            break;
        }
        case 4:
            //project forces into complement space, put in xPrime
            subspaceProjection(forces, xPrime, numberOfAtoms, masses, inverseMasses, 2);
            
            //find slope of PE with /lambda here
            RealOpenMM lambda = 1.0 / _maxEig;
            RealOpenMM oldLambda = lambda;
            std::stringstream message;
            message << "Lambda "<< lambda << " oldLambda " << oldLambda  <<"\n";
            
            SimTKOpenMMLog::printMessage( message );
            
            //just equal to minus dot product of 'proposed position move' and forces (=-\nabla PE)
            //RealOpenMM lambdaSlp = 0.0;
            //for( int ii = 0; ii < numberOfAtoms; ii++ ){
            //    for( int jj = 0; jj < 3; jj++ ){
            //        lambdaSlp -= xPrime[ii][jj] * forces[ii][jj] * inverseMasses[ii]; //posTemp[k].dot((*myForcesP)[k]);
            //    }
            //}
            
            //solve for minimum for quadratic fit using two PE vales and the slope with /lambda=0
            RealOpenMM a, b;
            
            //a = -(((currentPE - lastPE) / lambda - lambdaSlp) / lambda);
            //b = lambdaSlp - 2.0 * a * lambda;
            a = (((currentPE - lastPE) / lambda - lastSlope) / lambda);
            b = lastSlope;
            
            if( a != 0.0){
                lambda = -b / (2 * a);
            }else{
                lambda = oldLambda / 2.0;
            }
            
            //test if lambda negative, if so just use smaller lambda
            if(lambda <= 0.0) lambda = oldLambda / 2.0;
            
            //std::stringstream message;
            //message << "Lambda "<< lambda << " oldLambda " << oldLambda << " slope " << lambdaSlp << " lastSlope " << lastSlope<<"\n";
            message << "Lambda "<< lambda << " oldLambda " << oldLambda << " lastSlope " << lastSlope<<"\n";
            
            SimTKOpenMMLog::printMessage( message );
            
            // copy xPrime -> atomCoordinates
            for( int ii = 0; ii < numberOfAtoms; ii++ ){
                const RealOpenMM factor = inverseMasses[ii] * (lambda-oldLambda);
                
                xPrime[ii][0] = atomCoordinates[ii][0] + factor * xPrime[ii][0];
                xPrime[ii][1] = atomCoordinates[ii][1] + factor * xPrime[ii][1];
                xPrime[ii][2] = atomCoordinates[ii][2] + factor * xPrime[ii][2];
            }
            break;
    }
    
  //New NML minimizer
  
  std::stringstream message;
    message << "Array     "<<forces[0][0]<<" " <<forces[0][1]<<" " <<forces[0][2]<<" "<<xPrime[1][0]<<" " <<xPrime[1][1]<<" " <<xPrime[1][2]<<" " <<"\n";
    RealOpenMM* mydata = &forces[0][0];
    message << "Lin Array "<<mydata[0]<<" " <<mydata[1]<<" " <<mydata[2]<<" "<<mydata[3]<<" " <<mydata[4]<<" " <<mydata[5]<<" " <<"\n";
    forces[0][0] = 1.0; forces[0][1] = 1.1; xPrime[1][0] = 2.0;
    message << "Lin Array "<<mydata[0]<<" " <<mydata[1]<<" " <<mydata[2]<<" "<<mydata[3]<<" " <<mydata[4]<<" " <<mydata[5]<<" " <<"\n";
    
    SimTKOpenMMLog::printMessage( message );
    
  //End NML Minimizer
  

   if( debug ){
      int maxPrint = 5;
      std::stringstream message;
      message << methodName << " Post SD2 atoms=" << numberOfAtoms << "\n";
      RealOpenMM** oldVelocities   = get2DArrayAtIndex( OldV );
      for( int ii = 0; ii < maxPrint; ii++ ){
         message << " x[";
         SimTKOpenMMUtilities::formatRealStringStream( message, atomCoordinates[ii], 3, one );
         message << "] xp[";
         SimTKOpenMMUtilities::formatRealStringStream( message, xPrime[ii], 3, one );
         message << "] v[";
         SimTKOpenMMUtilities::formatRealStringStream( message, velocities[ii], 3, one );
         message << "] ov[";
         SimTKOpenMMUtilities::formatRealStringStream( message, oldVelocities[ii], 3, one );
         message << "]\n";
      }
      SimTKOpenMMLog::printMessage( message );
   }
#if 0
  ReferenceConstraintAlgorithm* referenceConstraintAlgorithm = getReferenceConstraintAlgorithm();
  if( referenceConstraintAlgorithm ){
    
    /*
     std::stringstream message;
     message << methodName;
     message << " calling constrain1\n";
     SimTKOpenMMLog::printMessage( message );
     */
    
    referenceConstraintAlgorithm->apply( numberOfAtoms, atomCoordinates, xPrime,
                                        inverseMasses );
    
  }
#endif  
   // copy xPrime -> atomCoordinates

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      atomCoordinates[ii][0] = xPrime[ii][0];
      atomCoordinates[ii][1] = xPrime[ii][1];
      atomCoordinates[ii][2] = xPrime[ii][2];
   }

   incrementTimeStep();

   return 0; // ReferenceDynamics::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Write state

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses
   @param state               0 if initial state; otherwise nonzero
   @param baseFileName        base file name

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceNMLDynamics::writeState( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                             RealOpenMM** velocities,
                                             RealOpenMM** forces, RealOpenMM* masses,
                                             int state, const std::string& baseFileName ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceNMLDynamics::writeState";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;
   static const int threeI            =  3;

   // ---------------------------------------------------------------------------------------

   std::stringstream stateFileName;

   stateFileName << baseFileName;
   stateFileName << "_Step" << getTimeStep();
   // stateFileName << "_State" << state;
   stateFileName << ".txt";

   // ---------------------------------------------------------------------------------------

   // open file -- return if unsuccessful

   FILE* stateFile = NULL;
#ifdef WIN32
   fopen_s( &stateFile, stateFileName.str().c_str(), "w" );
#else
   stateFile = fopen( stateFileName.str().c_str(), "w" );
#endif

   // ---------------------------------------------------------------------------------------

   // diagnostics

   if( stateFile != NULL ){
      std::stringstream message;
      message << methodName;
      message << " Opened file=<" << stateFileName.str() << ">.\n";
      SimTKOpenMMLog::printMessage( message );
   } else {
      std::stringstream message;
      message << methodName;
      message << " could not open file=<" << stateFileName.str() << "> -- abort output.\n";
      SimTKOpenMMLog::printMessage( message );
      return -1; // ReferenceDynamics::ErrorReturn;
   }   

   // ---------------------------------------------------------------------------------------

   StringVector scalarNameI;
   IntVector scalarI;

   StringVector scalarNameR;
   RealOpenMMVector scalarR;

   StringVector scalarNameR1;
   RealOpenMMPtrVector scalarR1;

   StringVector scalarNameR2;
   RealOpenMMPtrPtrVector scalarR2;

   scalarI.push_back( getNumberOfAtoms() );
   scalarNameI.push_back( "Atoms" );

   scalarI.push_back( getTimeStep() );
   scalarNameI.push_back( "Timestep" );

   if( state == 0 || state == -1 ){

      scalarR.push_back( getDeltaT() );
      scalarNameR.push_back( "delta_t" );

      scalarR.push_back( getTemperature() );
      scalarNameR.push_back( "T" );

      scalarR.push_back( getTau() );
      scalarNameR.push_back( "Tau" );

      scalarR1.push_back( masses );
      scalarNameR1.push_back( "mass" );

      scalarR2.push_back( atomCoordinates );
      scalarNameR2.push_back( "coord" );

      scalarR2.push_back( velocities );
      scalarNameR2.push_back( "velocities" );

      scalarR2.push_back( forces );
      scalarNameR2.push_back( "forces" );

      if( state == -1 ){

         RealOpenMM** xPrime          = get2DArrayAtIndex( xPrime2D );
         RealOpenMM** oldVelocities   = get2DArrayAtIndex( OldV );
         RealOpenMM** xVector         = get2DArrayAtIndex( X2D );
         RealOpenMM** vVector         = get2DArrayAtIndex( V2D );

         scalarR2.push_back( xPrime );
         scalarNameR2.push_back( "xPrime" );

         scalarR2.push_back( oldVelocities);
         scalarNameR2.push_back( "vold" );

         scalarR2.push_back( xVector );
         scalarNameR2.push_back( "xVector" );

         scalarR2.push_back( vVector );
         scalarNameR2.push_back( "vVector" );
      }
      
   } else {

      scalarR2.push_back( atomCoordinates );
      scalarNameR2.push_back( "coord" );

      scalarR2.push_back( velocities );
      scalarNameR2.push_back( "velocities" );

   }

   writeStateToFile( stateFile, scalarNameI, scalarI, scalarNameR, scalarR, getNumberOfAtoms(), scalarNameR1, scalarR1, threeI, scalarNameR2, scalarR2 ); 

   (void) fclose( stateFile );

   return 0; // ReferenceDynamics::DefaultReturn;

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Find forces OR positions inside subspace (defined as the span of the 'eigenvectors' Q)
// Take 'array' as input, 'outArray' as output (may be the same vector).
// 'projectionMode' determines:
// a) if the projection is for force or positions, bit 1 set=force.
// b) if the projection is for the sub-space or-complement space, bit 2 set=sub-space.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void ReferenceNMLDynamics::subspaceProjection(  RealOpenMM** arrayParam, 
                                                RealOpenMM** outArrayParam, 
                                                int numberOfAtoms,
                                                RealOpenMM* masses,
                                                RealOpenMM* inverseMasses,
                                                int projectionMode){
  
  //move to linear arrays for Blas etc.
  RealOpenMM* array = &arrayParam[0][0];
  RealOpenMM* outArray = &outArrayParam[0][0];
  
  //If 'array' and 'outArray are not the same array
  //copy 'array' into outArray
  const unsigned int _3N = numberOfAtoms * 3;
  if (array != outArray) {
    for (unsigned int i=0; i<_3N; i++)
      outArray[i] = array[i];
  }
  
  //Temporary array for the 'mode' values, denoted 'c'
  //the size is equal to the number of eigenvectors
  RealOpenMM* tmpC = new RealOpenMM[_numProjectionVectors];
  
  //~~~~We need to calculate M^{1/2}QQ^TM^{-1/2}force or M^{-1/2}QQ^TM^{1/2}positions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  //First we weight the array by multiplying by 
  //the square-root of the atomic masses for positions
  //or the the inverse of the square-root of the atomic masses for forces.
  //
  //a'=M^{-1/2}*f for forces, OR a'= M^{1/2}*x for positions
  //
  if(!(projectionMode & 0x01)) { //Forces if zero so weight is sqrt(1/m)
    
    //array is 3N long, for N atoms
    for( int i=0; i < numberOfAtoms; i++) {           //N loops
      RealOpenMM  massWt = SQRT( inverseMasses[i] );
      for (unsigned int j=0; j<3; j++)                //times 3 loops
        outArray[i*3 + j] *= massWt;
    }
    
  } else { //positions if first bit set, so weight is sqrt(m)
    
    //array is 3N long, for N atoms
    for( int i=0; i < numberOfAtoms; i++) {           //N loops
      RealOpenMM  massWt = SQRT( masses[i] );
      for (unsigned int j=0; j<3; j++)                //times 3 loops
        outArray[i*3 + j] *= massWt;
    }
    
  }
  
  //Project onto mode space by taking the matrix product of 
  //the transpose of the eigenvectors Q with the array.
  //
  //c=Q^T*a', a' from last algorithm step
  //
#if defined(HAVE_LAPACK)
  //~~~~Blas here~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  char *transA = "T";							// Transpose, LAPACK checks only first character N/V
  int m = numberOfAtoms * 3; //eigenvector vector/array length
  int n = _numProjectionVectors;  //number of eigenvectors
  int incxy = 1;	//sizes
  RealOpenMM alpha = 1.0;	RealOpenMM beta = 0.0;
  
  //if single precision use Blas 'sgemv'
  #if RealOpenMMType == 1 
    //
    sgemv_ (transA, &m, &n, &alpha, _projectionVectors, &m, outArray, &incxy, &beta, tmpC, &incxy);
    //
    
    //Now find projected force/positions a'' by matrix product with Eigenvectors Q
    //or the complementary Eigenvectors
    //a''=Qc
    char *transB = "N"; /* LAPACK checks only first character N/V */
    
    //if projectionMode has bit 2 set then complement space (I-Q^TQ) not (Q^TQ)
    if(!(projectionMode & 0x02)) {
      alpha = 1.0;	beta = 0.0;
    } else {
      alpha = -1.0;	beta = 1.0;
    }
    //
    sgemv_ (transB, &m, &n, &alpha, _projectionVectors, &m, tmpC, &incxy, &beta, outArray, &incxy);
    
    //if double precision use Blas 'dgemv'
  #else
    //
    dgemv_ (transA, &m, &n, &alpha, _projectionVectors), &m, outArray, &incxy, &beta, tmpC, &incxy);
    //
    
    //Now find projected force/positions a'' by matrix product with Eigenvectors Q
    //a''=Qc
    char *transB = "N"; /* LAPACK checks only first character N/V */
    
    //if projectionMode has bit 2 set then complement space (I-Q^TQ) not (Q^TQ)
    if(!(projectionMode & 0x02)) {
      alpha = 1.0;	beta = 0.0;
    } else {    
      alpha = -1.0;	beta = 1.0;
    }
    //
    dgemv_ (transB, &m, &n, &alpha, _projectionVectors, &m, tmpC, &incxy, &beta, outArray, &incxy);
  #endif
#else

  //~~~~NO Blas here~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //If no Blas is available we need to manually find the product c=A*b
  //c_i=\sum_{j=1}^n A_{i,j} b_j
  
  //c=Q^T*a', a' from last algorithm step
  //Q is a linear array in column major format
  //so tmpC_i = \sum_{j=1}^n Q_{j,i} outArray_j
  //Q_{j,i}=_projectionVectors[j * numberOfAtoms * 3 + i]
  
  //const int _3N = numberOfAtoms * 3; //eigenvector vector/array length
  
  for(int i=0; i < _numProjectionVectors; i++){   //over all eigenvectors
    
    tmpC[i] = 0.0;  //clear
    
    for(int j=0; j< _3N; j++){  //over each element in the vector
      //tmpC[i] += _projectionVectors[j * _3N + i] * outArray[j];
        tmpC[i] += _projectionVectors[j  + i * _3N] * outArray[j];
    }
  }
  
  //Now find projected force/positions a'' by matrix product with Eigenvectors Q
  //a''=Qc
  //so outArray_i  = \sum_{j=1}^n Q_{i,j} tmpC_i
  
  //if projectionMode has bit 2 set then complement space (I-Q^TQ) not (Q^TQ)
  const bool pMode = !(projectionMode & 0x02);
  
  //find product
  for(int i=0; i< _3N; i++){  //over each element in the vector
    
    //if sub-space do Q*c
    //else do a'-Q(Q^T a') = (I-QQ^T)a'
      if(pMode){
          outArray[i] = 0.0; //if not complement
    
          for(int j=0; j < _numProjectionVectors; j++){   //over all eigenvectors
              //outArray[i] += _projectionVectors[i * _3N + j] * tmpC[j];
              outArray[i] += _projectionVectors[i + j * _3N] * tmpC[j];
          }
      }else{
          for(int j=0; j < _numProjectionVectors; j++){   //over all eigenvectors
              //outArray[i] -= _projectionVectors[i * _3N + j] * tmpC[j];
              outArray[i] -= _projectionVectors[i + j * _3N] * tmpC[j];
          }
          
      }
    
  }
  
  //~~~~End of no-blas~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#endif
    
  //Finally we weight the array by multiplying by 
  //the inverse of the square-root of the atomic masses for positions
  //or the the square-root of the atomic masses for forces.
  //
  //a'''=M^{1/2}*a'' or a'''=M^{-1/2}*a''
  //
  if(!(projectionMode & 0x01)) { //Forces if zero so weight is sqrt(m)
    for( int i=0; i < numberOfAtoms; i++) {
      RealOpenMM  massUnWt = SQRT( masses[i] );
      for (unsigned int j=0; j<3; j++)
        outArray[i*3 + j] *= massUnWt;
    }
  } else {  //positions if first bit set, so weight is sqrt(1/m)
    for( int i=0; i < numberOfAtoms; i++) {
      RealOpenMM  massUnWt = SQRT( inverseMasses[i] );
      for (unsigned int j=0; j<3; j++)
        outArray[i*3 + j] *= massUnWt;
    }
  }
  
  //Clean up
  
  //delete temporary array
  delete [] tmpC;
  
}

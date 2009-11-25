#ifndef __ReferenceNMLDynamics_H__
#define __ReferenceNMLDynamics_H__

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

#include "SimTKUtilities/SimTKOpenMMRealType.h"
#include "SimTKReference/ReferenceDynamics.h"

// ---------------------------------------------------------------------------------------

class ReferenceNMLDynamics : public ReferenceDynamics {

   private:

      //dt, myGamma, fdt, vdt, ndt, sqrtFCoverM from Langevin Leapfrog
      //enum FixedParameters   { GDT, EPH, EMH, EP, EM, B, C, D, V, X, Yv, Yx, MaxFixedParameters };
      enum FixedParameters   { DT, GAMMA, FDT, VDT, NDT, SQRTFCOVERM, MaxFixedParameters };
      enum TwoDArrayIndicies { X2D, V2D, OldV, xPrime2D, vPrime2D, Max2DArrays };
      enum OneDArrayIndicies { InverseMasses, Max1DArrays };

      RealOpenMM _tau;
      RealOpenMM _fixedParameters[MaxFixedParameters];
      RealOpenMM* _projectionVectors;
      unsigned int _numProjectionVectors;
      RealOpenMM _minimumLimit, _maxEig;
      RealOpenMM lastPE, lastSlope;

      /**---------------------------------------------------------------------------------------
      
         Set fixed values
      
         @return ReferenceDynamics::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int _setFixedParameters( void );
      
   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor

         @param numberOfAtoms  number of atoms
         @param deltaT         delta t for dynamics
         @param tau            viscosity
         @param temperature    temperature
      
         --------------------------------------------------------------------------------------- */

       ReferenceNMLDynamics( int numberOfAtoms, RealOpenMM deltaT, RealOpenMM tau, RealOpenMM temperature,
                   RealOpenMM* projectionVectors, unsigned int numProjectionVectors, 
                   RealOpenMM minimumLimit, RealOpenMM maxEig
                    );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceNMLDynamics( );

      /**---------------------------------------------------------------------------------------
      
         Get tau
      
         @return tau
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getTau( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get array of fixed parameters indexed by 'FixedParameters' enums
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      const RealOpenMM* getFixedParameters( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Print parameters
      
         @param message message

         @return ReferenceDynamics::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int printParameters( std::stringstream& message ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Update
      
         @param numberOfAtoms       number of atoms
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses
      
         @return ReferenceDynamics::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
     
      int update( int numberOfAtoms, RealOpenMM** atomCoordinates,
                  RealOpenMM** velocities, RealOpenMM** forces, RealOpenMM* masses, const RealOpenMM currentPE, const int stepType );
     
    /**---------------------------------------------------------------------------------------
     
     Velocity update; Langevin Leapfrog
     
     @param numberOfAtoms       number of atoms
     @param atomCoordinates     atom coordinates
     @param velocities          velocities
     @param forces              forces
     @param masses              atom masses
     @param inverseMasses       inverse atom masses
     @param xVector             xVector
     
     --------------------------------------------------------------------------------------- */
    void halfKick( int numberOfAtoms, RealOpenMM** atomCoordinates, RealOpenMM** velocities,
                    RealOpenMM** forces, RealOpenMM* masses, RealOpenMM* inverseMasses, RealOpenMM** xVector );
    
    /**---------------------------------------------------------------------------------------
     
     Position update; Langevin Leapfrog
     
     @param numberOfAtoms       number of atoms
     @param atomCoordinates     atom coordinates
     @param velocities          velocities
     @param forces              forces
     @param masses              atom masses
     @param inverseMasses       inverse atom masses
     @param xPrime              xPrime
     
     --------------------------------------------------------------------------------------- */
    void drift( int numberOfAtoms, RealOpenMM** atomCoordinates, RealOpenMM** velocities,
                  RealOpenMM** forces, RealOpenMM* masses, RealOpenMM* inverseMasses,
                  RealOpenMM** xPrime );
          
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
      
      int writeState( int numberOfAtoms, RealOpenMM** atomCoordinates,
                      RealOpenMM** velocities, RealOpenMM** forces, RealOpenMM* masses,
                      int state, const std::string& baseFileName ) const;

      /**---------------------------------------------------------------------------------------
       Find forces OR positions inside subspace (defined as the span of the 'eigenvectors' Q)
       Take 'array' as input, 'outArray' as output (may be the same vector).
       'projectionMode' determines:
       a) if the projection is for force or positions, bit 1 set=force.
       b) if the projection is for the sub-space or-complement space, bit 2 set=sub-space.
       ----------------------------------------------------------------------------------------- */
      void subspaceProjection(  RealOpenMM** arrayParam, 
                                RealOpenMM** outArrayParam, 
                                int numberOfAtoms,
                                RealOpenMM* masses,
                                RealOpenMM* inverseMasses,
                                int projectionMode);
  
};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceNMLDynamics_H__

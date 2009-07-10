
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

#ifndef __ReferenceVariableStochasticDynamics_H__
#define __ReferenceVariableStochasticDynamics_H__

#include "ReferenceDynamics.h"

// ---------------------------------------------------------------------------------------

class ReferenceVariableStochasticDynamics : public ReferenceDynamics {

   private:

      enum FixedParameters   { GDT, EPH, EMH, EP, EM, EM_V, B, C, D, V, X, Yv, Yx, MaxFixedParameters };
      enum TwoDArrayIndicies { X2D, V2D, OldV, xPrime2D, vPrime2D, Max2DArrays };
      enum OneDArrayIndicies { InverseMasses, Max1DArrays };

      RealOpenMM _tau, _accuracy;
      RealOpenMM _fixedParameters[MaxFixedParameters];

      /**---------------------------------------------------------------------------------------

         Set fixed values

         @return ReferenceDynamics::DefaultReturn

         --------------------------------------------------------------------------------------- */

      int _setFixedParameters( RealOpenMM timeStep, RealOpenMM prevTimeStep );

   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         @param numberOfAtoms  number of atoms
         @param tau            viscosity
         @param temperature    temperature
         @param accuracy       required accuracy

         --------------------------------------------------------------------------------------- */

       ReferenceVariableStochasticDynamics( int numberOfAtoms, RealOpenMM tau, RealOpenMM temperature, RealOpenMM accuracy );

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceVariableStochasticDynamics( );

      /**---------------------------------------------------------------------------------------

         Get tau

         @return tau

         --------------------------------------------------------------------------------------- */

      RealOpenMM getTau( void ) const;
      
      /**---------------------------------------------------------------------------------------

         Get the required accuracy

         @return accuracy

         --------------------------------------------------------------------------------------- */

      RealOpenMM getAccuracy( void ) const;

      /**---------------------------------------------------------------------------------------

         Set the required accuracy

         --------------------------------------------------------------------------------------- */

      void setAccuracy( RealOpenMM accuracy );

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
         @param maxStepSize         maximum time step

         @return ReferenceDynamics::DefaultReturn

         --------------------------------------------------------------------------------------- */

      int update( int numberOfAtoms, RealOpenMM** atomCoordinates,
                  RealOpenMM** velocities, RealOpenMM** forces, RealOpenMM* masses, RealOpenMM maxStepSize );

      /**---------------------------------------------------------------------------------------

         First update; based on code in update.c do_update_sd() Gromacs 3.1.4

         @param numberOfAtoms       number of atoms
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses
         @param inverseMasses       inverse atom masses
         @param xPrime              xPrime
         @param oldVelocities       previous velocities
         @param xVector             xVector
         @param vVector             vVector
         @param maxStepSize         maximum time step

         @return ReferenceDynamics::DefaultReturn

         --------------------------------------------------------------------------------------- */

      int updatePart1( int numberOfAtoms, RealOpenMM** atomCoordinates, RealOpenMM** velocities,
                       RealOpenMM** forces, RealOpenMM* masses, RealOpenMM* inverseMasses,
                       RealOpenMM** xPrime, RealOpenMM** oldVelocities,
                       RealOpenMM** xVector, RealOpenMM** vVector, RealOpenMM maxStepSize );

      /**---------------------------------------------------------------------------------------

         Second update

         @param numberOfAtoms       number of atoms
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses

         @return ReferenceDynamics::DefaultReturn

         --------------------------------------------------------------------------------------- */

      int updatePart2( int numberOfAtoms, RealOpenMM** atomCoordinates, RealOpenMM** velocities,
                       RealOpenMM** forces, RealOpenMM* inverseMasses,
                       RealOpenMM** xPrime, RealOpenMM** oldVelocities,
                       RealOpenMM** xVector, RealOpenMM** vVector );
      
};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceVariableStochasticDynamics_H__

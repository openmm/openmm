
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

#ifndef __ReferenceShakeAlgorithm_H__
#define __ReferenceShakeAlgorithm_H__

// ---------------------------------------------------------------------------------------

class ReferenceShakeAlgorithm {

   protected:

      int _maximumNumberOfIterations;
      RealOpenMM _tolerance;

      int _numberOfConstraints;
      int** _atomIndices;
      RealOpenMM** _shakeParameters;

      RealOpenMM** _r_ij;
      RealOpenMM* _d_ij2;
      RealOpenMM* _distanceTolerance;
      RealOpenMM* _reducedMasses;
      
   public:

      /**---------------------------------------------------------------------------------------
      
         ReferenceShakeAlgorithm constructor
      
         @param numberOfConstraints      number of constraints 
         @param atomIndices              atom indices for contraints
         @param shakeParameters          Shake parameters
      
         --------------------------------------------------------------------------------------- */
      
      ReferenceShakeAlgorithm( int numberOfConstraints, int** atomIndices, RealOpenMM** shakeParameters );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceShakeAlgorithm( );

      /**---------------------------------------------------------------------------------------
      
         Get number of constraints
      
         @return number of constraints
      
         --------------------------------------------------------------------------------------- */
      
      int getNumberOfConstraints( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get maximum number of iterations
      
         @return maximum number of iterations
      
         --------------------------------------------------------------------------------------- */
      
      int getMaximumNumberOfIterations( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set maximum number of iterations
      
         @param maximumNumberOfIterations   new maximum number of iterations
         @return ReferenceDynamics::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int setMaximumNumberOfIterations( int maximumNumberOfIterations );
      
      /**---------------------------------------------------------------------------------------
      
         Get tolerance
      
         @return tolerance
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getTolerance( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set tolerance
      
         @param tolerance new tolerance

         @return tolerance
      
         --------------------------------------------------------------------------------------- */
      
      int setTolerance( RealOpenMM tolerance );
      
      /**---------------------------------------------------------------------------------------
      
         Print parameters
      
         @param message message

         @return ReferenceShakeAlgorithm::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
          
      int printParameters( std::stringstream& message ) const;

      /**---------------------------------------------------------------------------------------
      
         Apply Shake algorithm
      
         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param atomCoordinatesP atom coordinates prime
         @param inverseMasses    1/mass
      
         @return ReferenceDynamics::DefaultReturn if converge; else 
          return ReferenceDynamics::ErrorReturn
      
         --------------------------------------------------------------------------------------- */
      
      int applyShake( int numberOfAtoms, RealOpenMM** atomCoordinates,
                       RealOpenMM** atomCoordinatesP, RealOpenMM* inverseMasses );
      
      /**---------------------------------------------------------------------------------------
      
         Report any violated constraints
      
         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param message          report
      
         @return number of violated constraints
      
         --------------------------------------------------------------------------------------- */
      
      int reportShake( int numberOfAtoms, RealOpenMM** atomCoordinates, std::stringstream& message );
};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceShakeAlgorithm_H__


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

#ifndef __ReferenceConstraint_H__
#define __ReferenceConstraint_H__

// ---------------------------------------------------------------------------------------

class ReferenceConstraint {

   private:

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceConstraint( );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceConstraint( );

};

class ReferenceShakeConstraint : public ReferenceConstraint {

   private:

      // force heavy atom into index 0

      int _atomIndices[2];
      RealOpenMM _inverseMasses[2];
      RealOpenMM _constraintDistance;

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor

         @param atomIndex1     atom index 1
         @param atomIndex2     atom index 2
         @param distance       distance constraint
         @param inverseMass1   inverseMass of atom 1
         @param inverseMass2   inverseMass of atom 2
      
         --------------------------------------------------------------------------------------- */

       ReferenceShakeConstraint( int atomIndex1, int atomIndex2, RealOpenMM distance,
                                 RealOpenMM inverseMass1, RealOpenMM inverseMass2 );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceShakeConstraint( );

      /**---------------------------------------------------------------------------------------
      
         Get constraint distance
      
         @return distance
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getConstraintDistance( void );

      /**---------------------------------------------------------------------------------------
      
         Get inverse mass of heavy atom
      
         @return inverse mass
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getHeavyAtomInverseMass( void );
      
      /**---------------------------------------------------------------------------------------
      
         Get inverse mass of light atom
      
         @return inverse mass
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getLightAtomInverseMass( void );

      /**---------------------------------------------------------------------------------------
      
         Get index of heavy atom
      
         @return index
      
         --------------------------------------------------------------------------------------- */
      
      int getHeavyAtomIndex( void );

      /**---------------------------------------------------------------------------------------
      
         Get index of light atom
      
         @return index
      
         --------------------------------------------------------------------------------------- */
      
      int getLightAtomIndex( void );
      
      /**---------------------------------------------------------------------------------------
      
         Print state
      
         @param message             message
      
         @return ReferenceDynamics::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int printState( std::stringstream& message );
      
};

typedef std::vector<ReferenceShakeConstraint*> ShakeVector;
typedef ShakeVector::iterator ShakeVectorI;
typedef std::map<int, ShakeVector*> IntShakeMap;
typedef IntShakeMap::iterator IntShakeMapI;

// ---------------------------------------------------------------------------------------

#endif // __ReferenceShakeConstraint_H__

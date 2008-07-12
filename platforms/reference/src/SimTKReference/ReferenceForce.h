
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

#ifndef __ReferenceForce_H__
#define __ReferenceForce_H__

// ---------------------------------------------------------------------------------------

class ReferenceForce {

   private:
       
       static RealOpenMM periodicDifference(RealOpenMM val1, RealOpenMM val2, RealOpenMM period);

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceForce( );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceForce( );

      /**---------------------------------------------------------------------------------------
      
         Static variables
      
         --------------------------------------------------------------------------------------- */

       static const int XIndex             = 0;
       static const int YIndex             = 1;
       static const int ZIndex             = 2;
       static const int R2Index            = 3;
       static const int RIndex             = 4;
       static const int LastDeltaRIndex    = 5;

       static const int DeltaRMaxIndex     = 5;

       static const int DefaultReturn      = 0;
       static const int ErrorReturn        = -1;
   
      /**---------------------------------------------------------------------------------------
      
         Get deltaR and distance and distance**2 between atomI and atomJ (static method)
         deltaR: j - i
      
         @param atomCoordinatesI    atom i coordinates
         @param atomCoordinatesI    atom j coordinates
         @param deltaR              deltaX, deltaY, deltaZ, R2, R upon return
      
         @return ReferenceForce::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int getDeltaR( const RealOpenMM* atomCoordinatesI, const RealOpenMM* atomCoordinatesJ,
                            RealOpenMM* deltaR );
      
      /**---------------------------------------------------------------------------------------

         Get deltaR and distance and distance**2 between atomI and atomJ, assuming periodic
         boundary conditions (static method); deltaR: j - i

         @param atomCoordinatesI    atom i coordinates
         @param atomCoordinatesI    atom j coordinates
         @param boxSize             X, Y, and Z sizes of the periodic box
         @param deltaR              deltaX, deltaY, deltaZ, R2, R upon return

         @return ReferenceForce::DefaultReturn

         --------------------------------------------------------------------------------------- */

      static int ReferenceForce::getDeltaRPeriodic( const RealOpenMM* atomCoordinatesI, const RealOpenMM* atomCoordinatesJ,
                                             const RealOpenMM* boxSize, RealOpenMM* deltaR );

    /**---------------------------------------------------------------------------------------
      
         Get deltaR between atomI and atomJ (static method): deltaR: j - i
      
         @param atomCoordinatesI    atom i coordinates
         @param atomCoordinatesI    atom j coordinates
         @param deltaR              deltaX, deltaY, deltaZ upon return
      
         @return ReferenceForce::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int getDeltaROnly( const RealOpenMM* atomCoordinatesI, const RealOpenMM* atomCoordinatesJ,
                                RealOpenMM* deltaR );
      
      /**---------------------------------------------------------------------------------------
      
         Write coordinates, energies and forces (Simbios)
      
         @param numberOfAtoms       number of atoms
         @param atomsPerBond        atoms/bond (used to normalize total energy)
         @param atomCoordinates     atomic coordinates
         @param forces              forces
         @param energies            energies (optional)
         @param resultsFileName     output file name
         @param lengthConversion    length conversion (optional defaults to 1)
         @param forceConversion     force  conversion (optional defaults to 1)
         @param energyConversion    energy conversion (optional defaults to 1)
      
         @return ReferenceForce::DefaultReturn unless
                 file cannot be opened
                 in which case return ReferenceForce::ErrorReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int writeForces( int numberOfAtoms, int atomsPerBond, RealOpenMM** atomCoordinates,
                              RealOpenMM** forces, RealOpenMM* energies,
                              const std::string& resultsFileName,
                              RealOpenMM lengthConversion, RealOpenMM forceConversion,
                              RealOpenMM energyConversion );
      
};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceForce_H__

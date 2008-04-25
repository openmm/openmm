
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

#ifndef __GrycukParameters_H__
#define __GrycukParameters_H__

#include "SimTKOpenMMCommon.h"
#include "ImplicitSolventParameters.h"

// ---------------------------------------------------------------------------------------

class GrycukParameters : public ImplicitSolventParameters {

   public:

       static const std::string ParameterFileName;

   private:

      // scaled radius factors (S_kk in HCT paper)

      int _ownScaledRadiusFactors;
      RealOpenMM* _scaledRadiusFactors;

   public:

      /**---------------------------------------------------------------------------------------
      
         GrycukParameters constructor (Simbios) 
      
         @param numberOfAtoms       number of atoms
      
         --------------------------------------------------------------------------------------- */
      
       GrycukParameters( int numberOfAtoms );

      /**---------------------------------------------------------------------------------------
      
         GrycukParameters destructor (Simbios) 
      
         --------------------------------------------------------------------------------------- */
      
       ~GrycukParameters( );

      // override of new/delete

      //static void* operator new( size_t size );
      //static void  operator delete( void *p );

      //static void* operator new[]( size_t size );
      //static void  operator delete[]( void *p );

      /**---------------------------------------------------------------------------------------
      
         Return Grycuk scale factors
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      const RealOpenMM* getScaledRadiusFactors( void ) const;
        
      /**---------------------------------------------------------------------------------------
      
         Return Grycuk scale factors
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      int setScaledRadiusFactors( RealOpenMM* scaledRadiusFactors );
#if RealOpenMMType == 2
      int setScaledRadiusFactors( float* scaledRadiusFactors );
#endif
      int setScaledRadiusFactors( const RealOpenMMVector& scaledRadiusFactors );
        
      /**---------------------------------------------------------------------------------------
      
         Set flag indicating whether scale factors arra should be deleted
      
         @param ownScaledRadiusFactors flag indicating whether scale factors 
                                       array should be deleted

         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */

      int setOwnScaleFactors( int ownScaledRadiusFactors );
      
      /**--------------------------------------------------------------------------------------- 
      
         Assign standard radii for GB/SA methods other than ACE;
         taken from Macromodel and OPLS-AA, except for hydrogens (Simbios)
      
         Logic based on logic in Tinker's ksolv.f
      
         Currently only works for standard amino acid atoms
         If invalid atom name is encountered, a message is printed to log file and the
         radius for that atom is set to 1.0f
      
         @param numberOfAtoms       number of atoms
         @param atomNames           array of atom names from GMX top data struct
         @param radii               array to store Macromodel radii for each atom
         @param log                 if set, then print error messages to log file
      
         @return SimTKOpenMMCommon::DefaultReturn always
      
         --------------------------------------------------------------------------------------- */
      
      int getMacroModelAtomicRadii( int numberOfAtoms,
                                    char*** atomNames, RealOpenMM* radii, FILE* log );

      /**---------------------------------------------------------------------------------------
      
         Get AtomicRadii array w/ dielectric offset applied
      
         @return array of atom volumes
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM* getAtomicRadii( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set AtomicRadii array
      
         @param atomicRadii array of atomic radii
         @param units       units flag: SimTKOpenMMCommon::KcalAngUnits or
                                        SimTKOpenMMCommon::MdUnits 
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */

      int setAtomicRadii( RealOpenMM* atomicRadii, int units = SimTKOpenMMCommon::MdUnits  );

      /**---------------------------------------------------------------------------------------
      
         Set AtomicRadii array
      
         @param atomicRadii vector of atomic radii
         @param units       units flag: SimTKOpenMMCommon::KcalAngUnits or
                                        SimTKOpenMMCommon::MdUnits 
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */

      int setAtomicRadii( const RealOpenMMVector& atomicRadii, int units = SimTKOpenMMCommon::MdUnits );

      /**---------------------------------------------------------------------------------------
      
         Map Gmx atom name to Tinker atom number (Simbios)
      
         @param atomName            atom name (CA, HA, ...); upper and lower case should both work
         @param log                 if set, then print error messages to log file
      
         return Tinker atom number if atom name is valid; else return -1
      
         --------------------------------------------------------------------------------------- */
            
      int mapGmxAtomNameToTinkerAtomNumber( const char* atomName, FILE* log ) const;

      /**---------------------------------------------------------------------------------------
            
         Get string w/ state
         
         @param title               title (optional)
            
         @return string
            
         --------------------------------------------------------------------------------------- */
      
      std::string getStateString( const char* title ) const;

      /**---------------------------------------------------------------------------------------
                  
         Return zero value if all parameters set; else return nonzero
               
         @return ready status
                  
         --------------------------------------------------------------------------------------- */
       
      int isNotReady( void ) const;
      
};
   
/**---------------------------------------------------------------------------------------
      
   Qsort/heapsort integer comparison (Simbios) 
      
   @parma a first value to compare
   @param b second value to compare

   @return -1, 0, 1
      
--------------------------------------------------------------------------------------- */

int integerComparison( const void *a, const void *b);

// ---------------------------------------------------------------------------------------

#endif // __GrycukParameters_H__

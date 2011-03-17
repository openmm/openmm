
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

#ifndef __ImplicitSolventParameters_H__
#define __ImplicitSolventParameters_H__

#include "../SimTKUtilities/SimTKOpenMMRealType.h"
#include "../../../../openmmapi/include/openmm/internal/windowsExport.h"
#include <string>

// ---------------------------------------------------------------------------------------

class OPENMM_EXPORT ImplicitSolventParameters {

   protected:

      // parameters common to implicit solvent models
   
      RealOpenMM _solventDielectric;
      RealOpenMM _soluteDielectric;
      RealOpenMM _kcalA_To_kJNm;
      RealOpenMM _electricConstant;
      RealOpenMM _probeRadius;
      RealOpenMM _pi4Asolv;

      RealOpenMM _preFactor;
   
      // ---------------------------------------------------------------------------------------

      // atom count

      int _numberOfAtoms;

      // flag signalling whether arrays are to be freed

      int _freeArrays;

      // atomic radii

      int _ownAtomicRadii;
      RealOpenMM* _atomicRadii;

      /**---------------------------------------------------------------------------------------
      
         Initialize ImplicitSolvent Parameters (Simbios) 
      
         --------------------------------------------------------------------------------------- */
      
      void _initializeImplicitSolventConstants( void );

      /**---------------------------------------------------------------------------------------
      
         Set KcalA_To_kJNm
      
         @param kcalA_To_kJNm probe radius
      
         --------------------------------------------------------------------------------------- */
      
      void setKcalA_To_kJNm( RealOpenMM kcalA_To_kJNm );

      /**---------------------------------------------------------------------------------------
      
         Set free array flag -- if set then work arrays are freed when destructor is called
      
         @param freeArrays flag
      
         --------------------------------------------------------------------------------------- */

      void setFreeArrays( int freeArrays );

      /**---------------------------------------------------------------------------------------
      
         Reset prefactor (Simbios) 
      
         called when _electricConstant, _soluteDielectric, or _solventDielectric are modified
      
         --------------------------------------------------------------------------------------- */
      
      void _resetPreFactor( void );
      
   public:

      /**---------------------------------------------------------------------------------------
      
         ImplicitSolventParameters constructor (Simbios) 
      
         @param numberOfAtoms       number of atoms
      
         --------------------------------------------------------------------------------------- */
    
       ImplicitSolventParameters( int numberOfAtoms );

      /**---------------------------------------------------------------------------------------
      
         ImplicitSolventParameters destructor (Simbios) 
      
         --------------------------------------------------------------------------------------- */
    
       virtual ~ImplicitSolventParameters( );

      // override of new/delete

      //static void* operator new( size_t size );
      //static void  operator delete( void *p );

      //static void* operator new[]( size_t size );
      //static void  operator delete[]( void *p );

      /**---------------------------------------------------------------------------------------
      
         Get number of atoms
      
         @return number of atoms
      
         --------------------------------------------------------------------------------------- */
    
      int getNumberOfAtoms( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Get free array flag -- if set then work arrays are freed when destructor is called
      
         @return freeArrays flag
      
         --------------------------------------------------------------------------------------- */
      
      int getFreeArrays( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Get solvent dielectric
      
         @return solvent dielectric
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getSolventDielectric( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set solvent dielectric
      
         @param solventDielectric solvent dielectric
      
         --------------------------------------------------------------------------------------- */
      
      void setSolventDielectric( RealOpenMM solventDielectric );
      
      /**---------------------------------------------------------------------------------------
      
         Get solute dielectric
      
         @return soluteDielectric
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getSoluteDielectric( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set solute dielectric
      
         @param soluteDielectric solute dielectric
      
         --------------------------------------------------------------------------------------- */
      
      void setSoluteDielectric( RealOpenMM soluteDielectric );
      
      /**---------------------------------------------------------------------------------------
      
         Get conversion factor for kcal/A -> kJ/nm (Simbios) 
      
         @return kcalA_To_kJNm factor
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getKcalA_To_kJNm( void ) const;
      
      
      /**---------------------------------------------------------------------------------------
      
         Get electric constant (Simbios) 
      
         @return electricConstant
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getElectricConstant( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set electric constant (Simbios) 
      
         @param electricConstant       electric constant
      
         --------------------------------------------------------------------------------------- */
      
      void setElectricConstant( RealOpenMM electricConstant );

      /**---------------------------------------------------------------------------------------
      
         Get probe radius (Simbios) 
      
         @return probeRadius
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getProbeRadius( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set probe radius (Simbios) 
      
         @param probeRadius   probe radius
      
         --------------------------------------------------------------------------------------- */
      
      void setProbeRadius( RealOpenMM probeRadius );
      
      /**---------------------------------------------------------------------------------------
      
         Get pi4Asolv:  used in ACE approximation for nonpolar term  
            ((RealOpenMM) M_PI)*4.0f*0.0049f*1000.0f; (Simbios) 
      
         @return pi4Asolv
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getPi4Asolv( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get prefactor
      
         @returnpreFactor 
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getPreFactor( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set values used in ACE approximation for nonpolar term
            ((RealOpenMM) M_PI)*4.0f*0.0049f*1000.0f; (Simbios) 
      
         @param pi4Asolv   see above
      
         --------------------------------------------------------------------------------------- */
      
      void setPi4Asolv( RealOpenMM pi4Asolv );
      
      /**---------------------------------------------------------------------------------------
      
         Get AtomicRadii array
      
         @return array of atom volumes
      
         --------------------------------------------------------------------------------------- */

      virtual RealOpenMM* getAtomicRadii( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set AtomicRadii array
      
         @param atomicRadii array of atomic radii
      
         --------------------------------------------------------------------------------------- */

      virtual void setAtomicRadii( RealOpenMM* atomicRadii );

      /**---------------------------------------------------------------------------------------
      
         Set AtomicRadii array
      
         @param atomicRadii vector of atomic radii
      
         --------------------------------------------------------------------------------------- */

      virtual void setAtomicRadii( const RealOpenMMVector& atomicRadii );
      
      /**---------------------------------------------------------------------------------------
      
         Set flag indicating whether AtomicRadii array is owned by
         object; if flag is set, then when the object is deleted,
         the array will be freed
      
         @param ownAtomicRadii flag indicating whether array of atomic radii
                               should be freed
      
         --------------------------------------------------------------------------------------- */

      void setOwnAtomicRadii( int ownAtomicRadii );

      /**---------------------------------------------------------------------------------------
            
         Print state to log file (Simbios)
         
         @param title               title (optional)
            
         --------------------------------------------------------------------------------------- */
      
      virtual std::string getStateString( const char* title ) const;

      /**---------------------------------------------------------------------------------------
            
         Get string tab -- centralized
         
         @return tab string
            
         --------------------------------------------------------------------------------------- */
      
      std::string getStringTab( void ) const;
      
};

// ---------------------------------------------------------------------------------------

#endif // __ImplicitSolventParameters_H__

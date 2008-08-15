#ifndef OPENMM_BROOK_COMMON_H_
#define OPENMM_BROOK_COMMON_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
 * Contributors:                                                              *
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

#include <vector>
#include <set>

#include "BrookFloatStreamInternal.h"
#include "BrookIntStreamInternal.h"
#include "BrookPlatform.h"

namespace OpenMM {

/**
 * This kernel is invoked by StandardMMForceField to calculate the forces acting on the system.
 */
class BrookCommon {

   public:
  
      // return values

      static const int DefaultReturnValue = 0;
      static const int ErrorReturnValue   = -1;

      BrookCommon( );
  
      ~BrookCommon();
  
      /**
       * Return number of atoms
       * 
       * @return number of atoms
       *
       */

      int getNumberOfAtoms( void ) const; 

      /** 
       * Get atom ceiling parameter
       * 
       * @return atom ceiling parameter
       *
       */
      
      int getAtomSizeCeiling( void ) const;
      
      /**
       * Get atom stream width
       *
       * @param platform platform reference
       *
       * @return atom stream width
       */

      int getAtomStreamWidth( const Platform& platform ); 

      /**
       * Get atom stream width
       *
       * @return atom stream width
       */

      int getAtomStreamWidth( void ) const; 

      /**
       * Get atom stream height
       *
       * @param platform platform reference
       *
       * @return atom stream height
       */

      int getAtomStreamHeight( const Platform& platform ); 

      /**
       * Get atom stream height
       *
       * @return atom stream height
       */

      int getAtomStreamHeight( void ) const;

      /**
       * Get atom stream size
       * 
       * @param platform platform reference
       *
       * @return atom stream size
       */

      int getAtomStreamSize( const Platform& platform ); 

      /**
       * Get atom stream size
       * 
       * @return atom stream size
       */

      int getAtomStreamSize( void ) const; 

      /** 
       * Set log file reference
       * 
       * @param  log file reference
       *
       * @return DefaultReturnValue
       *
       */
      
      int setLog( FILE* log );

      /* 
       * Get contents of object
       *
       * @param level of dump
       *
       * @return string containing contents
       *
       * */
      
      std::string getContents( int level ) const;

      /** 
       * Get log file reference
       * 
       * @return  log file reference
       *
       */
      
      FILE* getLog( void ) const;
      
   protected:
   
      // number of atoms

      int _numberOfAtoms;

      // atom stream dimensions

      int _atomStreamWidth;
      int _atomStreamHeight;
      int _atomStreamSize;

      // atom size mod

      int _atomSizeModified;

      // log file reference

      FILE* _log;

      /**
       * Set number of atoms
       * 
       * @param numberOfAtoms number of atoms
       *
       */

      int setNumberOfAtoms( int numberOfAtoms ); 

      /** 
       * Get atom stream dimensions
       * 
       * @param platform                  platform
       *
       */
      
      void _getAtomStreamDimensions( const Platform& platform );
      
      /* 
       * Get line
       *
       * @param tab         tab
       * @param description description
       * @param value       value
       *
       * @return string containing contents
       *
       * */
      
      std::string _getLine( const std::string& tab, const std::string& description, const std::string& value ) const;
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_COMMON_H_ */

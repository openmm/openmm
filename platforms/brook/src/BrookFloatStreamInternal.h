#ifndef OPENMM_BROOK_FLOAT_STREAM_INTERNAL_H_
#define OPENMM_BROOK_FLOAT_STREAM_INTERNAL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright ( c ) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files ( the "Software" ), *
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

#include "BrookPlatform.h"
#include "BrookStreamInternal.h"
#include "../../reference/src/SimTKUtilities/SimTKOpenMMRealType.h"

namespace OpenMM {

/**
 * This is the implementation of Float and Double streams in the Brook Platform.
 */

class BrookFloatStreamInternal : public BrookStreamInternal {

   public:

      /** 
       * BrookFloatStreamInternal constructor
       * 
       * @param name                      stream name
       * @param size                      size of array
       * @param streamWidth               stream width
       * @param type                      stream type (float, float2, ...)
       * @param inputDefaultDangleValue   default dangle value
       *
       */
      
      BrookFloatStreamInternal( const std::string& name, int size, int inputStreamWidth, BrookStreamInternal::DataType type, double defaultDangleValue = 0.0 );

      /** 
       * BrookFloatStreamInternal destructor
       * 
       */
      
      ~BrookFloatStreamInternal(  );

      /** 
       * Load data from input array
       * 
       * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
       * and to contain elements of the correct _data type for this stream.  If the stream has a compound _data type, all
       * the values should be packed into a single array: all the values for the first element, followed by all the values
       * for the next element, etc.
       *
       * @throw exception if baseType not float or double
       *
       */

      void loadFromArray( const void* array );

      /** 
       * Copy the contents of an array into this stream.
       * 
       * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
       * and to contain elements of the correct data type for this stream.  If the stream has a compound data type, all
       * the values should be packed into a single array: all the values for the first element, followed by all the values
       * for the next element, etc.
       */
      void loadFromArray( const void* array, BrookStreamInternal::DataType baseType );

      /** 
       * Save data to input array
       * 
       * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
       * and to contain elements of the correct _data type for this stream.  If the stream has a compound _data type, all
       * the values should be packed into a single array: all the values for the first element, followed by all the values
       * for the next element, etc.
       *
       * @throw exception if baseType not float or double
       *
       */

      void saveToArray( void* array );

      /** 
       * Fill data w/ input value
       * 
       * @param  value to set array to
       *
       *
       */

      void fillWithValue( void* value );

      /** 
       * Get data array -- no read from board!
       * 
       * @return  data array
       */
      
      void* getData( void );

      /** 
       * Get data array
       * 
       * @param readFromBoard if set, read values on board 
       *
       * @return  data array
       */
      
      void* getData( int readFromBoard );

      /** 
       * Get array of appropritate size for loading data
       *
       * @return data array -- user's responsibility to free
       */
      void* getDataArray( void );
  
      /** 
       * Get dangle value
       * 
       * @return  dangle value
       */
      
      double getDangleValue( void ) const;

      /*  
       * Get contents of object
       *
       *
       * @param level   level of dump
       *
       * @return string containing contents
       *
       * */

      const std::string getContentsString( int level = 0 ) const;

   private:

      BrookOpenMMFloat _defaultDangleValue;

      float* _data;

      RealOpenMM* _realOpenMMData;
  
      void _loadDanglingValues( void  );
      void _loadDanglingValues( float );

      /*  
       * Print array to file
       *
       * @param log  log file
       *
       * @return  DefaultReturnValue
       *
       * */

      int _bodyPrintToFile( FILE* log );
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_FLOAT_STREAM_INTERNAL_H_ */

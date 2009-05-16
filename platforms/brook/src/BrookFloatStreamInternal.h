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
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs, Mike Houston                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
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

      /** 
       * Sum over stream dimensions
       * 
       * @param stopIndex                 index to stop sum
       * @param sum                       array of size=getWidth()
       *
       * @return DefaultReturnValue
       *
       * @throw exception if stopIndex is too large
       */
      
      int sumByDimension( int stopIndex, double* sum );
      
      /*  
       * Get stats
       *
       * @return statistics vector
       *
       * */
     
      int getStatistics( std::vector<std::vector<double>>&, int maxScan );

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

      int _bodyPrintToFile( FILE* log, int maxPrint );
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_FLOAT_STREAM_INTERNAL_H_ */

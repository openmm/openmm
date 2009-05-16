#ifndef OPENMM_BROOK_STREAM_IMPL_H_
#define OPENMM_BROOK_STREAM_IMPL_H_

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

#include "openmm/StreamImpl.h"
#include "BrookFloatStreamInternal.h"
#include "BrookIntStreamInternal.h"

namespace OpenMM {

/**
 * This is the base class of Float and Double streams in the Brook Platform.
 */

class BrookStreamImpl : public StreamImpl {

   public:

      /**
       * BrookStreamImpl constructor
       * 
       * @param name        name of the stream to create
       * @param size        number of elements in the stream
       * @param streamWidth stream width
       * @param type        data type of each element in the stream
       * @param platform    platform
       *
       */
  
      BrookStreamImpl( const std::string& name, int size, int streamWidth, Stream::DataType type, const Platform& platform );
  
      /**
       * BrookStreamImpl destructor
       * 
       */

      ~BrookStreamImpl( );
  
      /** 
       * BrookStreamImpl constructor
       *
       * @param type        StreamImpl data type (float, float2, ...)
       * @param isFloat     on output = 1 if float
       *                                2 if double
       *                                0 if integer
       *
       * @return BrookStreamInternal::DataType mapping to Stream::DataType type
       *         if no match, return BrookStreamInternal::Unknown
       * 
       */
      
      BrookStreamInternal::DataType getTypeMap( Stream::DataType type, int* isFloat ) const;

      /**
       * Copy the contents of an array into this stream.
       * 
       * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
       * and to contain elements of the correct data type for this stream.  If the stream has a compound data type, all
       * the values should be packed into a single array: all the values for the first element, followed by all the values
       * for the next element, etc.
       */
      void loadFromArray( const void* array );
  
      /**
       * Copy the contents of this stream into an array.
       * 
       * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
       * and to contain elements of the correct data type for this stream.  If the stream has a compound data type, all
       * the values should be packed into a single array: all the values for the first element, followed by all the values
       * for the next element, etc.
       */
      void saveToArray( void* array );
  
      /**
       * Set every element of this stream to the same value.
       * 
       * @param a pointer to the value.  It is assumed to be of the correct data type for this stream.
       */
      void fillWithValue( void* value );
  
      /**
       * Get data array
       * 
       * @return data array
       */
      void* getData( void );
      void* getData( int readFromBoard );

      /**
       * Get Brook stream
       * 
       * @return Brook stream reference 
       */
      brook::stream& getBrookStream( void );
  
      /** 
       * Get width
       * 
       * @return width
       */

      int getWidth( void ) const;
    
      /** 
       * Get stream width
       * 
       * @return stream width
       */

      int getStreamWidth( void ) const;

      /** 
       * Get stream height
       * 
       * @return stream height
       */

      int getStreamHeight( void ) const;

      /** 
       * Get stream size
       * 
       * @return stream size
       */

      int getStreamSize( void ) const;

      /* 
       * Get contents of object
       *
       *
       * @param level   level of dump
       *
       * @return string containing contents
       *
       * */

      //const std::string getContentsString( int level = 0 ) const;

      /** 
       * Get Brook stream impl 
       * 
       * @return Brook stream impl
       */

      BrookStreamInternal* getBrookStreamInternal( void ) const;

   protected:

      BrookStreamInternal* _brookStreamInternal; 
  
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_STREAM_IMPL_H_ */

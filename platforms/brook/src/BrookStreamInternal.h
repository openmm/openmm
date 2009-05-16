#ifndef OPENMM_BROOK_STREAM_INTERNAL_H_
#define OPENMM_BROOK_STREAM_INTERNAL_H_

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

#include <brook/brook.hpp>

namespace OpenMM {

/**
 * This is the implementation of Float and Double streams in the Brook Platform.
 */

class BrookStreamInternal {

   public:

      /** 
        * This is an enumeration of the allowed data types for a Stream.
        */
      enum DataType { Float, Float2, Float3, Float4, Double, Double2, Double3, Double4, Integer, Integer2, Integer3, Integer4, Unknown };

      // return values

      static const int DefaultReturnValue = 0;
      static const int ErrorReturnValue   = -1; 

      // ---------------------------------------------------------------------------------------

      /**
       * BrookStreamInternal constructor
       * 
       * @param name        name of the stream to create
       * @param size        number of elements in the stream
       * @param streamWidth stream width
       * @param type        data type of each element in the stream
       *
       */
  
      BrookStreamInternal( const std::string& name, int size, int streamWidth, BrookStreamInternal::DataType type );
  
      /**
       * BrookStreamInternal destructor
       * 
       */

      ~BrookStreamInternal( );
  
      /**
       * Get the name of this stream.
       */
  
      const std::string& getName( void ) const;
  
      /**
       * Get the number of elements in this stream.
       */
  
      int getSize( void ) const;
  
      /**
       * Get the data type of each element in the stream.
       */
  
      BrookStreamInternal::DataType getDataType( void ) const;
  
      /**
       * Get base data type of each element in the stream ( float, double, int )
       */
  
      BrookStreamInternal::DataType getBaseDataType( void ) const;
  
      /**
       * Copy the contents of an array into this stream.
       * 
       * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
       * and to contain elements of the correct data type for this stream.  If the stream has a compound data type, all
       * the values should be packed into a single array: all the values for the first element, followed by all the values
       * for the next element, etc.
       */
      virtual void loadFromArray( const void* array ) = 0;
  
      /**
       * Copy the contents of an array into this stream.
       * 
       * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
       * and to contain elements of the correct data type for this stream.  If the stream has a compound data type, all
       * the values should be packed into a single array: all the values for the first element, followed by all the values
       * for the next element, etc.
       */
      virtual void loadFromArray( const void* array, BrookStreamInternal::DataType baseType ) = 0;
  
      /**
       * Copy the contents of this stream into an array.
       * 
       * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
       * and to contain elements of the correct data type for this stream.  If the stream has a compound data type, all
       * the values should be packed into a single array: all the values for the first element, followed by all the values
       * for the next element, etc.
       */
      virtual void saveToArray( void* array ) = 0;
  
      /**
       * Set every element of this stream to the same value.
       * 
       * @param a pointer to the value.  It is assumed to be of the correct data type for this stream.
       */
      virtual void fillWithValue( void* value ) = 0;
  
      /**
       * Get data
       * 
       * @return data array
       */
      virtual void* getData( void ) = 0;
  
      /**
       * Get data
       *
       * @param readFromBoard  read data from board
       * 
       * @return data array
       */
      virtual void* getData( int readFromBoard ) = 0;
  
      /**
       * Get array of appropritate size for loading data
       *
       * @return data array -- user's responsibility to free
       */
      virtual void* getDataArray( void ) = 0;
  
      /** 
       * Get type string
       *
       * @param type        BrookStreamInternal data type (float, float2, ...)
       *
       * @return string matching type or "Unknown"
       * 
       */
      
      std::string getTypeString( BrookStreamInternal::DataType type ) const;
      
      /** 
       * Get Brook stream reference
       * 
       * @return  Brook stream reference
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

      const std::string getContentsString( int level = 0 ) const;

      /* 
       * Print to file
       *
       * @param log         log file
       * @param maxPrint    max values to print; if < 0, then all values printed; default value is -1
       *
       * @return  DefaultReturnValue
       *
       * */

      int printToFile( FILE* log, int maxPrint = -1 );

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
       * @param statistics  output vector of stats
       * @param maxScan     number of points to use in computing stats
       *
       * @return statistics vector
       *
       * */
     
      virtual int getStatistics( std::vector<std::vector<double>>& statVector, int maxScan );

      /* 
       * Get stat string
       *
       * @param         tag    id tag
       * @param  statistics    stat vector
       *
       * @return stat string
       *
       **/
      
      std::string printStatistics( std::string tag, std::vector<std::vector<double> >& statistics ) const;
      
      /* 
       * Read stream from file
       *
       * @param fileName     file name
       *
       * @return DefaultReturnValue or ErrorReturnValue if problems
       *
       **/
      
      int loadStreamGivenFileName( std::string& filename );
      
      /* 
       * Print streams to file
       *
       * @param fileName     file name
       * @param streams      streams to print
       *
       * @return DefaultReturnValue
       *
       **/
       
      static int printStreamsToFile( std::string fileName, std::vector<BrookStreamInternal*>& streams );

      /* 
       * Check for NANs
       *
       * @return number of Nans found
       *
       **/
      
      int checkForNans( void );
      
      /* 
       * Sum columns
       *
       * @param  sums   output vector of column sums
       *
       * @return DefaultReturnValue
       *
       **/
      
      int sumColumns( std::vector<float>& sums );
      
      
   protected:

      std::string _name;

      BrookStreamInternal::DataType _type;
      BrookStreamInternal::DataType _baseType;

      int _size;
      int _width;
      int _streamWidth;
      int _streamHeight;
      int _streamSize;
  
      brook::stream _aStream;
  
      /* 
       * Get contents of object
       *
       * @param tab         tab
       * @param description description
       * @param value       value
       *
       * @return string containing contents
       *
       * */
      
      std::string _getLine( const std::string& tab, const std::string& description, 
                            const std::string& value ) const;
      
      /* 
       * Print array to file
       *
       * @param log  log file
       *
       * @return  DefaultReturnValue
       *
       * */

      virtual int _bodyPrintToFile( FILE* log,  int maxPrint ) = 0;
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_STREAM_INTERNAL_H_ */

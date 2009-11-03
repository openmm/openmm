
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

#ifndef __TwoDimArraySimTk_H_
#define __TwoDimArraySimTk_H_

// ---------------------------------------------------------------------------------------

//#include "SimTKOpenMMUtilities.h"
#include "SimTKOpenMMUtilities.h"

/**---------------------------------------------------------------------------------------

   Class for 2D arrays

--------------------------------------------------------------------------------------- */

template <typename T> class TwoDimArraySimTk {

   private:
 
      int _rowSize;
      int _columnSize;
      T** _2Darray;
      T* _2DMemoryBlock;
      std::string _name;

   public:

      TwoDimArraySimTk( int rowSize, int columnSize );

      ~TwoDimArraySimTk();

      /**---------------------------------------------------------------------------------------
      
         Get 2D array
         
         @return ptr to array
         
         --------------------------------------------------------------------------------------- */
      
      T** get2DArray( void ) const { return _2Darray; };

      /**---------------------------------------------------------------------------------------
      
         Get print string
         
         @return  std::string
         
         --------------------------------------------------------------------------------------- */
      
      std::string getPrintString( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Get name
         
         @return array name
         
         --------------------------------------------------------------------------------------- */
      
      std::string getName( void ) const { return _name; };

      /**---------------------------------------------------------------------------------------
      
         Set name
         
         @param array name
         
         --------------------------------------------------------------------------------------- */
      
      void setName( std::string name ){ _name = name; };

};
   
/**---------------------------------------------------------------------------------------

   TwoDimArraySimTk constructor

   @param rowSize                row dimension
   @param columnSize             column dimension

   --------------------------------------------------------------------------------------- */

template <typename T>
TwoDimArraySimTk<T>::TwoDimArraySimTk( int rowSize, int columnSize ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nTwoDimArraySimTk<T>::TwoDimArraySimTk";

   // ---------------------------------------------------------------------------------------

   _name              = "NotSet";
   _rowSize           = rowSize;
   _columnSize        = columnSize;

#ifdef useXMalloc 
   _2Darray           = (T**) SimTKOpenMMUtilities::Xmalloc( "TwoDimArraySimTk", __FILE__, __LINE__, rowSize*sizeof( T* ) );
   _2DMemoryBlock     = (T*)  SimTKOpenMMUtilities::Xmalloc( "block",   __FILE__, __LINE__, columnSize*rowSize*sizeof( T ) );
#else
   _2Darray           = (T**) SimTKOpenMMUtilities::Xmalloc( "TwoDimArraySimTk", __FILE__, __LINE__, rowSize*sizeof( T* ) );
   _2DMemoryBlock     = (T*)  SimTKOpenMMUtilities::Xmalloc( "block",   __FILE__, __LINE__, columnSize*rowSize*sizeof( T ) );
#endif

   T* blockPtr        = _2DMemoryBlock;

   for( int ii = 0; ii < rowSize; ii++ ){
      _2Darray[ii]    = blockPtr;
      blockPtr       += columnSize;
   }

   memset( _2DMemoryBlock, 0, columnSize*rowSize*sizeof( T ) );

// (void) fprintf( stdout, "\nTwoDimArraySimTk %s ", getPrintString().c_str() ); 
// (void) fflush( stdout );

}

/**---------------------------------------------------------------------------------------

   TwoDimArraySimTk destructor

   array[i][j]

   @param rowSize                row dimension
   @param columnSize             column dimension

   --------------------------------------------------------------------------------------- */

template <typename T> TwoDimArraySimTk<T>::~TwoDimArraySimTk( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nTwoDimArraySimTk<T>::~TwoDimArraySimTk";

   // ---------------------------------------------------------------------------------------

   if( _2DMemoryBlock ){
      SimTKOpenMMUtilities::Xfree(  "2DMemoryBlock", __FILE__, __LINE__, _2DMemoryBlock );
   }
   if( _2Darray ){
      SimTKOpenMMUtilities::Xfree( "2Darray", __FILE__, __LINE__, _2Darray );
   }
      
}

/**---------------------------------------------------------------------------------------

   TwoDimArraySimTk print string 

   @return id string

   --------------------------------------------------------------------------------------- */

template <typename T> std::string TwoDimArraySimTk<T>::getPrintString( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nTwoDimArraySimTk<T>::getPrintString";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   message << _name << " [ " << _rowSize << ", " << _columnSize << " ] sizeT=" << sizeof(  T );
   message << " array=" <<  (unsigned int) _2Darray << " block=" << (unsigned int) _2DMemoryBlock;
   // message << " array=" <<  (void *) _2Darray << " block=" << (void*) _2DMemoryBlock;
   return message.str();
      
}

/**---------------------------------------------------------------------------------------

   Typdefs to make code more readable

   --------------------------------------------------------------------------------------- */

typedef TwoDimArraySimTk<Real> TwoDimRealArraySimTk;
typedef std::vector<TwoDimRealArraySimTk*> TwoDimRealArraySimTkVector;
typedef std::vector<TwoDimRealArraySimTk*>::iterator TwoDimRealArraySimTkVectorI;

#endif // __TwoDimArraySimTk_H_

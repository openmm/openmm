
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

#ifndef __GbsaAtomParameter_H__
#define __GbsaAtomParameter_H__

// STL includes

#include <vector>
#include <map>
#include <set>
#include <string>
#include "CommonSimTk.h"
#include "RealTypeSimTk.h"

//#ifndef Real
//#define Real float
//#endif


/*
class StringComparisonForMap : public std::binary_function< const std::string&, const std::string&, bool > {

   public:

      bool operator()( const std::string& str1, const std::string& str2 ) const {
         return strcmp( str1.c_str(), str2.c_str() ) < 0;
      }
}; 

typedef std::vector<std::string> StringVector;
typedef StringVector::iterator StringVectorI;
typedef StringVector::const_iterator StringVectorCI;
*/

// ---------------------------------------------------------------------------------------

class GbsaAtomParameter {

   private:

      // type id -- int Gromacs 'amber99_10', for example

      std::string _typeId;

      // vdw radius

      Real _vdwRadius;

   public:

      /**---------------------------------------------------------------------------------------
      
         GbsaParameters constructor (Simbios) 
      
         @param parameterLineTokens tokens from parameter file
      
         --------------------------------------------------------------------------------------- */
      
      GbsaAtomParameter( const StringVector& parameterLineTokens );

      /**---------------------------------------------------------------------------------------
      
         GbsaAtomParameter destructor (Simbios) 
      
         --------------------------------------------------------------------------------------- */
      
      ~GbsaAtomParameter( );
   
      /**---------------------------------------------------------------------------------------
      
         Get type id (Simbios) 
      
         @return type id
      
         --------------------------------------------------------------------------------------- */
      
      std::string getTypeId( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set type id (Simbios) 
      
         @param typeId type id
      
         @return 0
      
         --------------------------------------------------------------------------------------- */
      
      int setTypeId( std::string typeId );
      
      /**---------------------------------------------------------------------------------------
      
         Get vdw radius (Simbios) 
      
         @return vdw radius
      
         --------------------------------------------------------------------------------------- */
      
      Real getVdwRadius( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set vdw radius (Simbios) 
      
         @param vdwRadius new vdw radius
      
         @return 0
      
         --------------------------------------------------------------------------------------- */
      
      int setVdwRadius( Real vdwRadius );

      /**---------------------------------------------------------------------------------------
      
         Set vdw radius (Simbios) 
      
         @param vdwRadius new vdw radius
      
         @return 0
      
         --------------------------------------------------------------------------------------- */
      
      int setVdwRadius( const std::string& vdwRadius );
};

typedef std::map<std::string, GbsaAtomParameter*, StringComparisonForMap> GbsaAtomParameterMap;
typedef GbsaAtomParameterMap::iterator GbsaAtomParameterMapI;
typedef GbsaAtomParameterMap::const_iterator GbsaAtomParameterMapCI;

// ---------------------------------------------------------------------------------------

#endif // __GbsaAtomParameter_H__

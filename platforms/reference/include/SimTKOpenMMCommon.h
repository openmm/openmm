
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

#ifndef __SimTKOpenMMCommon_H__
#define __SimTKOpenMMCommon_H__

// include file containing entries commonly used

// STL includes

#include <vector>
#include <string>

// ---------------------------------------------------------------------------------------

#include "RealVec.h"

// ---------------------------------------------------------------------------------------

typedef std::vector<RealOpenMM> RealOpenMMVector;
typedef RealOpenMMVector::iterator RealOpenMMVectorI;
typedef RealOpenMMVector::const_iterator RealOpenMMVectorCI;

// ---------------------------------------------------------------------------------------

class SimTKOpenMMCommon {

   public:

      static const std::string NotSet;
      static const RealOpenMM BigCutoffValue;
      static const std::string Comment;
      static const std::string Tab;
      static const std::string YesU;
      static const std::string YesUl;
      static const std::string YesL;

      // subroutine returns

      static const int DefaultReturn;
      static const int ErrorReturn;
};

#endif // __SimTKOpenMMCommon_H__

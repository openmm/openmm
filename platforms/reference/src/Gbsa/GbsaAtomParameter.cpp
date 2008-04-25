
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

#include <math.h>

// #include "UtilitiesSimTk.h"
#include "GbsaAtomParameter.h"

/**---------------------------------------------------------------------------------------

   GbsaAtomParameter:

		Calculates for each atom

			(1) the van der Waal radii
         (2) volume
         (3) fixed terms in Gbsa equation gPol
         (4) list of atoms that should be excluded in calculating
				 force -- nonbonded atoms (1-2, and 1-3 atoms)

   --------------------------------------------------------------------------------------- */

/**---------------------------------------------------------------------------------------

   GbsaParameters constructor (Simbios) 

   @param parameterLineTokens tokens from parameter file

   --------------------------------------------------------------------------------------- */

GbsaAtomParameter::GbsaAtomParameter( const StringVector& parameterLineTokens ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaAtomParameter::GbsaAtomParameter";
   
   // ---------------------------------------------------------------------------------------

   if( parameterLineTokens.size() < 2 ){
      // XXX
      return;
   }
   setTypeId( parameterLineTokens[0] );
   setVdwRadius( parameterLineTokens[1] );

}

/**---------------------------------------------------------------------------------------

   GbsaAtomParameter destructor (Simbios) 

   --------------------------------------------------------------------------------------- */

GbsaAtomParameter::~GbsaAtomParameter( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaAtomParameter::~GbsaAtomParameter";
   
   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Get type id (Simbios) 

   @return type id

   --------------------------------------------------------------------------------------- */

std::string GbsaAtomParameter::getTypeId( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaAtomParameter::getTypeId:";

   // ---------------------------------------------------------------------------------------

   return _typeId;

}

/**---------------------------------------------------------------------------------------

   Set type id (Simbios) 

   @param typeId type id

   @return 0

   --------------------------------------------------------------------------------------- */

int GbsaAtomParameter::setTypeId( std::string typeId ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaAtomParameter::setTypeId:";

   // ---------------------------------------------------------------------------------------

   _typeId = typeId; 
   return 0;

}

/**---------------------------------------------------------------------------------------

   Get vdw radius (Simbios) 

   @return vdw radius

   --------------------------------------------------------------------------------------- */

Real GbsaAtomParameter::getVdwRadius( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaAtomParameter::getVdwRadius:";

   // ---------------------------------------------------------------------------------------

   return _vdwRadius;

}

/**---------------------------------------------------------------------------------------

   Set vdw radius (Simbios) 

   @param vdwRadius new vdw radius

   @return 0

   --------------------------------------------------------------------------------------- */

int GbsaAtomParameter::setVdwRadius( Real vdwRadius ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaAtomParameter::setVdwRadius:";

   // ---------------------------------------------------------------------------------------

   _vdwRadius = vdwRadius; 
   return 0;

}

/**---------------------------------------------------------------------------------------

   Set vdw radius (Simbios) 

   @param vdwRadius new vdw radius (string)

   @return 0

   --------------------------------------------------------------------------------------- */

int GbsaAtomParameter::setVdwRadius( const std::string& vdwRadius ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaAtomParameter::setVdwRadius(s):";

   // ---------------------------------------------------------------------------------------

    // bool SimTKOpenMMUtilities::isValidReal( std::string stringToCheck );
    return setVdwRadius( (Real) atof( vdwRadius.c_str() ) );

}


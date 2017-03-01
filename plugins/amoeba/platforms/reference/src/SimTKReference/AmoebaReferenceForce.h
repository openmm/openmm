
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

#ifndef __AmoebaReferenceForce_H__
#define __AmoebaReferenceForce_H__

#include "openmm/Vec3.h"
#include <vector>

namespace OpenMM {

class AmoebaReferenceForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    AmoebaReferenceForce();
 
    /**---------------------------------------------------------------------------------------
       
          Destructor
       
       --------------------------------------------------------------------------------------- */
 
    ~AmoebaReferenceForce();
 
 
    /**---------------------------------------------------------------------------------------
    
       Load delta of two vectors
    
       @param xVector      first vector
       @param yVector      second vector
       @param deltaR       output vector: y - x
    
       --------------------------------------------------------------------------------------- */
    
    static void loadDeltaR(const OpenMM::Vec3& xVector, const OpenMM::Vec3& yVector,
                           std::vector<double>& deltaR);
    /**---------------------------------------------------------------------------------------

       Load delta of two vectors, applying periodic boundary conditions

       @param xVector      first vector
       @param yVector      second vector
       @param deltaR       output vector: y - x
       @param boxVectors   periodic box vectors

       --------------------------------------------------------------------------------------- */

    static void loadDeltaRPeriodic(const Vec3& xVector, const Vec3& yVector, std::vector<double>& deltaR, const Vec3* boxVectors);
    
    /**---------------------------------------------------------------------------------------
    
       Calculate norm squared of 3d vector
    
       @param inputVector    vector whose norm squared is to be computed
    
       @return norm squared
    
       --------------------------------------------------------------------------------------- */
    
    static double getNormSquared3(const std::vector<double>& inputVector);
    static double getNormSquared3(const double* inputVector);
    
    /**---------------------------------------------------------------------------------------
    
       Calculate norm of 3d vector
    
       @param inputVector            vector whose norm is to be computed
    
       @return norm
    
       --------------------------------------------------------------------------------------- */
    
    static double getNorm3(const std::vector<double>& inputVector);
    static double getNorm3(const double* inputVector);
    
    /**---------------------------------------------------------------------------------------
    
       Normalize 3d vector
    
       @param inputVector            vector to normalize

       @return norm
    
       --------------------------------------------------------------------------------------- */
    
    static double normalizeVector3(double* inputVector);
    
    /**---------------------------------------------------------------------------------------
    
       Calculate dot product of 3d vectors
    
       @param xVector   first vector
       @param yVector   second vector
    
       @return dot product
    
       --------------------------------------------------------------------------------------- */
    
    static double getDotProduct3(const std::vector<double>& xVector, const std::vector<double>& yVector);
    static double getDotProduct3(const double* xVector,              const double* yVector);
    static double getDotProduct3(const double* xVector,              const OpenMM::Vec3& yVector);
    static double getDotProduct3(unsigned int vectorOffset, const std::vector<double>& xVector, const double* yVector);
    
    /**---------------------------------------------------------------------------------------
    
       Calculate z = x X y
    
       @param xVector      input vector
       @param yVector      input vector
       @param zVector      output vector: z = x X y
    
       --------------------------------------------------------------------------------------- */
    
    static void getCrossProduct(const std::vector<double>& xVector, const std::vector<double>& yVector,
                                std::vector<double>& zVector);
    
    static void getCrossProduct(const double* xVector, const double* yVector, double* zVector);


};

} // namespace OpenMM

#endif // _AmoebaReferenceForce___

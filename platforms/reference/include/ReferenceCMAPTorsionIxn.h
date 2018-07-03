
/* Portions copyright (c) 2010-2016 Stanford University and Simbios.
 * Contributors: Peter Eastman
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

#ifndef __ReferenceCMAPTorsionIxn_H__
#define __ReferenceCMAPTorsionIxn_H__

#include "SimTKOpenMMUtilities.h"
#include "ReferenceBondIxn.h"
#include <vector>

namespace OpenMM {

class ReferenceCMAPTorsionIxn : public ReferenceBondIxn {

private:

    std::vector<std::vector<std::vector<double> > > coeff;
    std::vector<int> torsionMaps;
    std::vector<std::vector<int> > torsionIndices;
    bool usePeriodic;
    Vec3 boxVectors[3];

    /**---------------------------------------------------------------------------------------

       Calculate the interaction due to a single torsion pair

       @param index            the index of the torsion
       @param atomCoordinates  atom coordinates
       @param forces           force array (forces added)
       @param totalEnergy      total energy

         --------------------------------------------------------------------------------------- */

    void calculateOneIxn(int index, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& forces,
                         double* totalEnergy) const;

public:

    /**---------------------------------------------------------------------------------------

       Constructor

       --------------------------------------------------------------------------------------- */

    ReferenceCMAPTorsionIxn(const std::vector<std::vector<std::vector<double> > >& coeff,
                            const std::vector<int>& torsionMaps,
                            const std::vector<std::vector<int> >& torsionIndices);

       /**---------------------------------------------------------------------------------------
      
         Set the force to use periodic boundary conditions.
      
         @param vectors    the vectors defining the periodic box
      
         --------------------------------------------------------------------------------------- */
      
       void setPeriodic(OpenMM::Vec3* vectors);

    /**---------------------------------------------------------------------------------------

       Calculate torsion interaction

       @param atomCoordinates    atom coordinates
       @param forces             force array (forces added)
       @param totalEnergy        total energy

         --------------------------------------------------------------------------------------- */

    void calculateIxn(std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& forces, double* totalEnergy) const;

    /**---------------------------------------------------------------------------------------

       This is present only because we must define it to subclass ReferenceBondIxn.  It is never called.

       --------------------------------------------------------------------------------------- */

    void calculateBondIxn(std::vector<int>& atomIndices, std::vector<OpenMM::Vec3>& atomCoordinates,
                         std::vector<double>& parameters, std::vector<OpenMM::Vec3>& forces,
                         double* totalEnergy, double* energyParamDerivs);

// ---------------------------------------------------------------------------------------

};

} // namespace OpenMM

#endif // __ReferenceCMAPTorsionIxn_H__


/* Portions copyright (c) 2006-2016 Stanford University and Simbios.
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

#ifndef __AmoebaReferenceTorsionTorsionForce_H__
#define __AmoebaReferenceTorsionTorsionForce_H__

#include "openmm/Vec3.h"
#include <vector>

namespace OpenMM {

class AmoebaReferenceTorsionTorsionForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    AmoebaReferenceTorsionTorsionForce() : usePeriodic(false) {};
 
    /**---------------------------------------------------------------------------------------
       
       Destructor
       
          --------------------------------------------------------------------------------------- */
 
    ~AmoebaReferenceTorsionTorsionForce() {};

    /**---------------------------------------------------------------------------------------

       Set the force to use periodic boundary conditions.
      
       @param vectors    the vectors defining the periodic box
      
       --------------------------------------------------------------------------------------- */
      
    void setPeriodic(OpenMM::Vec3* vectors);

     /**---------------------------------------------------------------------------------------
     
        Calculate Amoeba torsion-torsion ixns (force and energy)
     
        @param numTorsions             number of torsions
        @param posData                 particle positions
        @param particle1               particle 1 indices
        @param particle2               particle 2 indices
        @param particle3               particle 3 indices
        @param particle4               particle 4 indices
        @param particle5               particle 5 indices
        @param chiralCheckAtom         chiralCheckAtom (-1 if no chiral check is to be performed)
        @param gridIndices             grid indices
        @param torsionTorsionGrids     torsion-torsion grids
        @param forces                  output force vector
     
        @return total energy

     
        --------------------------------------------------------------------------------------- */

    double calculateForceAndEnergy(int numTorsionTorsions, std::vector<OpenMM::Vec3>& posData,
                                   const std::vector<int>&  particle1,
                                   const std::vector<int>&  particle2,
                                   const std::vector<int>&  particle3,
                                   const std::vector<int>&  particle4,
                                   const std::vector<int>&  particle5,
                                   const std::vector<int>&  chiralCheckAtom,
                                   const std::vector<int>&  gridIndices,
                                   const std::vector< std::vector< std::vector< std::vector<double> > > >& torsionTorsionGrids,
                                   std::vector<OpenMM::Vec3>& forceData) const;

private:

    bool usePeriodic;
    Vec3 boxVectors[3];

    /**---------------------------------------------------------------------------------------
    
       Load grid values from rectangle enclosing angles
    
       @param angle1    angle in first dimension
       @param angle2    angle in second dimension
    
       @param corners   on return contains the coordinates of the rectangle
    
       @param fValues   on return contains the values of the function at the vertices of the rectangle
       @param fValues1  on return contains the values of the derivative of the function wrt first dimension
       @param fValues2  on return contains the values of the derivative of the function wrt second dimension
       @param fValues12 on return contains the values of the derivative of the function wrt first & second dimension
    
       On first call a check is performed to see if the grid is valid
    
       --------------------------------------------------------------------------------------- */
    
    void loadGridValuesFromEnclosingRectangle(
               const std::vector< std::vector< std::vector<double> > >& grid,
               double angle1, double angle2, double corners[2][2],
               double* fValues, double* fValues1, double* fValues2, double* fValues12) const;
    
    /**---------------------------------------------------------------------------------------
    
       Determines the coefficient matrix needed for bicubic
       interpolation of a function, gradients and cross derivatives
    
       Reference:
    
         W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
         Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
         University Press, 1992, Section 3.6
    
       @param y       y          values
       @param y1      dy/dx1     values
       @param y2      dy/dx2     values
       @param y12     d2y/dx1dx2 values
       @param d1      d1Upper - d1Lower
       @param d2      d2Upper - d2Lower
       @param  c      4x4 return coefficient matrix
    
       --------------------------------------------------------------------------------------- */
    
    void getBicubicCoefficientMatrix(const double* y, const double* y1, const double* y2, const double* y12,
                                     const double d1, const double d2, double c[4][4]) const;
    
     /**---------------------------------------------------------------------------------------
     
        Determines the coefficient matrix needed for bicubic
        interpolation of a function, gradients and cross derivatives
     
        Reference:
     
          W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
          Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
          University Press, 1992, Section 3.6
     
        @param y       y          values
        @param y1      dy/dx1     values
        @param y2      dy/dx2     values
        @param y12     d2y/dx1dx2 values
     
        @param x1Upper upper x1
        @param x1Lower lower x1
     
        @param x2Upper upper x2
        @param x2Lower lower x2
     
        @param  gridValue1 grid value 1
        @param  gridValue2 grid value 2
     
        @param functionValue   function value (energy)
        @param functionValue1  d(energy)/dx1
        @param functionValue2  d(energy)/dx2
     
        --------------------------------------------------------------------------------------- */
    
    void getBicubicValues(
               const double* y, const double* y1, const double* y2, const double* y12,
               const double x1Lower, const double x1Upper,
               const double x2Lower, const double x2Upper,
               const double gridValue1, const double gridValue2,
               double* functionValue, double* functionValue1, double* functionValue2) const;
     
    /**---------------------------------------------------------------------------------------
     
        Tests the attached atoms at a torsion-torsion central
        site and returns -1 if the site is chiral; else return 1
     
        @param atomA   a-atom
        @param atomB   b-atom
        @param atomC   c-atom (central atom)
        @param atomD   d-atom
     
        @return 1.0 or -1.0 depending on whether chirality has an inverted sign
     
        --------------------------------------------------------------------------------------- */
    
    int checkTorsionSign(const OpenMM::Vec3& positionAtomA, const OpenMM::Vec3& positionAtomB,
                         const OpenMM::Vec3& positionAtomC, const OpenMM::Vec3& positionAtomD) const;
    
    /**---------------------------------------------------------------------------------------
    
       Calculate Amoeba torsion-torsion ixn (force and energy)
    
       @param positionAtomA           Cartesian coordinates of atom A
       @param positionAtomB           Cartesian coordinates of atom B
       @param positionAtomC           Cartesian coordinates of atom C
       @param positionAtomD           Cartesian coordinates of atom D
       @param positionAtomE           Cartesian coordinates of atom E
       @param positionChiralCheckAtom Cartesian coordinates of atom to be used in chiral check;
                                      if NULL, then no check is performed 
       @param grid                    grid vector
       @param forces                  force vector
    
       @return energy
    
       --------------------------------------------------------------------------------------- */
    
    double calculateTorsionTorsionIxn(const OpenMM::Vec3& positionAtomA, const OpenMM::Vec3& positionAtomB,
                                      const OpenMM::Vec3& positionAtomC, const OpenMM::Vec3& positionAtomD,
                                      const OpenMM::Vec3& positionAtomE, const OpenMM::Vec3* chiralCheckAtom,
                                      const std::vector< std::vector< std::vector<double> > >& grid,
                                      OpenMM::Vec3* forces) const;
         
};

} // namespace OpenMM

#endif // _AmoebaReferenceTorsionTorsionForce___

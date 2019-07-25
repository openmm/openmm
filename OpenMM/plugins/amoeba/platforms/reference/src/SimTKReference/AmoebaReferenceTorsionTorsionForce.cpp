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

#include "AmoebaReferenceForce.h"
#include "AmoebaReferenceTorsionTorsionForce.h"

using std::vector;
using namespace OpenMM;

void AmoebaReferenceTorsionTorsionForce::setPeriodic(OpenMM::RealVec* vectors) {
    usePeriodic = true;
    boxVectors[0] = vectors[0];
    boxVectors[1] = vectors[1];
    boxVectors[2] = vectors[2];
}

/**---------------------------------------------------------------------------------------

   Load grid values from rectenclosing angles

   @param grid        grid values [ angle1, angle2, f, df1, df2, df12 ] 
   @param angle1      angle in first dimension
   @param angle2      angle in second dimension

   @param corners   on return contains the coordinates of the rectangle

   @param fValues   on return contains the values of the function at the vertices of the rectangle
   @param fValues1  on return contains the values of the derivative of the function wrt first dimension
   @param fValues2  on return contains the values of the derivative of the function wrt second dimension
   @param fValues12 on return contains the values of the derivative of the function wrt first & second dimension

   On first call a check is performed to see if the grid is valid

   --------------------------------------------------------------------------------------- */

void AmoebaReferenceTorsionTorsionForce::loadGridValuesFromEnclosingRectangle(
           const std::vector< std::vector< std::vector<RealOpenMM> > >& grid,
           RealOpenMM angle1, RealOpenMM angle2, RealOpenMM corners[2][2],
           RealOpenMM* fValues, RealOpenMM* fValues1, RealOpenMM* fValues2, RealOpenMM* fValues12) const {
 
    // ---------------------------------------------------------------------------------------
 
    // static const std::string methodName = "loadGridValuesFromEnclosingRectangle";
 
    // ---------------------------------------------------------------------------------------

    // get 2 opposing grid indices for rectangle

    unsigned int gridSize  = grid.size();
    double gridSpacingI    = static_cast<double>(gridSize-1)/360.0;
    const int X_gridIndex  =  (int) ((angle1 - grid[0][0][0])*gridSpacingI + 1.0e-06);
    const int Y_gridIndex  =  (int) ((angle2 - grid[0][0][1])*gridSpacingI + 1.0e-06);

    // get coordinates of corner indices
  
    corners[0][0]         = grid[X_gridIndex][Y_gridIndex ][0]; 
    corners[0][1]         = grid[X_gridIndex+1][Y_gridIndex][0]; 

    corners[1][0]         = grid[X_gridIndex][Y_gridIndex  ][1]; 
    corners[1][1]         = grid[X_gridIndex+1][Y_gridIndex+1][1]; 

    for (int ii = 0; ii < 4; ii++) {
       int gridX = X_gridIndex;
       int gridY = Y_gridIndex;
       if (ii == 1) {
           gridX++;
       } else if (ii == 2) {
           gridX++;
           gridY++;
       } else if (ii == 3) {
           gridY++;
       }
           
       fValues[ii]   = grid[gridX][gridY][2];
       fValues1[ii]  = grid[gridX][gridY][3];
       fValues2[ii]  = grid[gridX][gridY][4];
       fValues12[ii] = grid[gridX][gridY][5];

    }

    return;

}

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

void AmoebaReferenceTorsionTorsionForce::getBicubicCoefficientMatrix(const RealOpenMM* y,
        const RealOpenMM* y1, const RealOpenMM* y2, const RealOpenMM* y12, const RealOpenMM d1, const RealOpenMM d2,
        RealOpenMM c[4][4]) const {

   // ---------------------------------------------------------------------------------------

    //static const std::string methodName   = "AmoebaReferenceTorsionTorsionForce::getBicubicCoefficientMatrix";

    static const RealOpenMM zero          = 0.0;
    static const RealOpenMM one           = 1.0;
    static const RealOpenMM two           = 2.0;
    static const RealOpenMM three         = 3.0;
    static const RealOpenMM four          = 4.0;
    static const RealOpenMM five          = 5.0;
    static const RealOpenMM six           = 6.0;
    static const RealOpenMM nine          = 9.0;
 
    // transpose of matrix in Tinker due to difference in C/Fotran row/column major
    // change indices when multiplying by weightMatrix
  
    static const RealOpenMM weightMatrix[16][16] = {
      {  one,  zero, -three,  two, zero, zero, zero,   zero, -three,   zero,   nine,   -six,  two, zero,   -six,  four },
      { zero,  zero, zero,   zero, zero, zero, zero,   zero,  three,   zero,  -nine,    six, -two, zero,    six, -four },
      { zero,  zero, zero,   zero, zero, zero, zero,   zero,   zero,   zero,   nine,   -six, zero, zero,   -six,  four },
      { zero,  zero, three,  -two, zero, zero, zero,   zero,   zero,   zero,  -nine,    six, zero, zero,    six, -four },
      { zero,  zero, zero,  zero,   one, zero, -three, two,    -two,   zero,    six,  -four,  one, zero, -three,   two },
      { zero,  zero, zero,  zero,  zero, zero, zero,   zero,   -one,   zero,  three,   -two,  one, zero, -three,   two },
      { zero,  zero, zero,  zero,  zero, zero, zero,   zero,   zero,   zero, -three,    two, zero, zero,  three,  -two },
      { zero,  zero, zero,  zero,  zero, zero, three,  -two,   zero,   zero,   -six,   four, zero, zero,  three,  -two },
      { zero,   one, -two,   one,  zero, zero, zero,   zero,   zero, -three,    six, -three, zero,  two,  -four,   two },
      { zero,  zero, zero,  zero,  zero, zero, zero,   zero,   zero,  three,   -six,  three, zero, -two,   four,  -two },
      { zero,  zero, zero,  zero,  zero, zero, zero,   zero,   zero,   zero, -three,  three, zero, zero,    two,  -two },
      { zero,  zero, -one,   one,  zero, zero, zero,   zero,   zero,   zero,  three, -three, zero, zero,   -two,   two },
      { zero,  zero, zero,  zero,  zero, one,  -two,    one,   zero,   -two,   four,   -two, zero,  one,   -two,   one },
      { zero,  zero, zero,  zero,  zero, zero, zero,   zero,   zero,   -one,    two,   -one, zero,  one,   -two,   one },
      { zero,  zero, zero,  zero,  zero, zero, zero,   zero,   zero,   zero,    one,   -one, zero, zero,   -one,   one },
      { zero,  zero, zero,  zero,  zero, zero, -one,    one,   zero,   zero,    two,   -two, zero, zero,   -one,   one } };
 
    // ---------------------------------------------------------------------------------------
 
    // pack y, y1, y2, y12 into single vector of dimension 16
 
    std::vector<RealOpenMM> x(16);
    RealOpenMM d1d2 = d1*d2;
    for (int ii = 0; ii < 4; ii++) {
       x[ii]     =   y[ii];
       x[ii+4]   =  y1[ii]*d1;
       x[ii+8]   =  y2[ii]*d2;
       x[ii+12]  = y12[ii]*d1d2;
    }
 
    // matrix multiply by the transpose of the stored weight table
 
    int rowIndex = 0;
    int colIndex = 0;
    for (int ii = 0; ii < 16; ii++) {
       RealOpenMM sum = weightMatrix[0][ii]*x[0];
       for (int jj = 1; jj < 16; jj++) {
          sum += weightMatrix[jj][ii]*x[jj];
       }
       c[rowIndex][colIndex++] = sum;
       if (!(colIndex % 4)) {
          colIndex = 0;
          rowIndex++;
       }
    }
 }
 
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
 
void AmoebaReferenceTorsionTorsionForce::getBicubicValues(
         const RealOpenMM* y, const RealOpenMM* y1, const RealOpenMM* y2, const RealOpenMM* y12,
         const RealOpenMM x1Lower, const RealOpenMM x1Upper,
         const RealOpenMM x2Lower, const RealOpenMM x2Upper, 
         const RealOpenMM gridValue1, const RealOpenMM gridValue2,
         RealOpenMM* functionValue, RealOpenMM* functionValue1, RealOpenMM* functionValue2) const { 
 
    // ---------------------------------------------------------------------------------------
 
    // static const std::string methodName = "getBicubicValues";

    static const RealOpenMM zero          = 0.0;
    static const RealOpenMM two           = 2.0;
    static const RealOpenMM three         = 3.0;
 
    // ---------------------------------------------------------------------------------------
 
    // get coefficent matrix
 
    RealOpenMM coefficientMatrix[4][4];
    getBicubicCoefficientMatrix(y, y1, y2, y12, x1Upper-x1Lower, x2Upper-x2Lower, coefficientMatrix);
 
    // apply coefficent matrix
 
    RealOpenMM t = (gridValue1 - x1Lower)/(x1Upper - x1Lower);
    RealOpenMM u = (gridValue2 - x2Lower)/(x2Upper - x2Lower);
 
    *functionValue  = zero;
    *functionValue1 = zero;
    *functionValue2 = zero;
 
    for (int ii = 3; ii >= 0; ii--) {
       *functionValue   = t*(*functionValue)   + (     (coefficientMatrix[ii][3]*u +     coefficientMatrix[ii][2])*u + coefficientMatrix[ii][1])*u + coefficientMatrix[ii][0];
       *functionValue1  = u*(*functionValue1)  + (three*coefficientMatrix[3][ii]*t + two*coefficientMatrix[2][ii])*t + coefficientMatrix[1][ii];
       *functionValue2  = t*(*functionValue2)  + (three*coefficientMatrix[ii][3]*u + two*coefficientMatrix[ii][2])*u + coefficientMatrix[ii][1];
    }
 
    *functionValue1 /= (x1Upper - x1Lower);
    *functionValue2 /= (x2Upper - x2Lower);
 
    return;
 }
 
 /**---------------------------------------------------------------------------------------
 
    Tests the attached atoms at a torsion-torsion central
    site and returns -1 if the site is chiral; else return 1
 
    @param atomA   a-atom
    @param atomB   b-atom
    @param atomC   c-atom (central atom)
    @param atomD   d-atom
 
    @return 1.0 or -1.0 depending on whether chirality has an inverted sign
 
    --------------------------------------------------------------------------------------- */
 
int AmoebaReferenceTorsionTorsionForce::checkTorsionSign(
         const RealVec& positionAtomA, const RealVec& positionAtomB,
         const RealVec& positionAtomC, const RealVec& positionAtomD) const {
 
    // ---------------------------------------------------------------------------------------
 
    // static const std::string methodName = "AmoebaReferenceTorsionTorsionForce::checkTorsionSign";
 
    static const RealOpenMM zero          = 0.0;
    static const int one                  = 1;

    // ---------------------------------------------------------------------------------------
 
    // compute parallelpiped volume at atomC and return sign based on sign of volume
  
    enum { CA, CB, CD, LastDeltaIndex };
    std::vector<RealOpenMM> deltaR[LastDeltaIndex];
    for (unsigned int ii = 0; ii < LastDeltaIndex; ii++) {
        deltaR[ii].resize(3);
    }   

    if (usePeriodic) {
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomC, positionAtomA, deltaR[CA], boxVectors);
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomC, positionAtomB, deltaR[CB], boxVectors);
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomC, positionAtomD, deltaR[CD], boxVectors);
    }
    else {
        AmoebaReferenceForce::loadDeltaR(positionAtomC, positionAtomA, deltaR[CA]);
        AmoebaReferenceForce::loadDeltaR(positionAtomC, positionAtomB, deltaR[CB]);
        AmoebaReferenceForce::loadDeltaR(positionAtomC, positionAtomD, deltaR[CD]);
    }

    RealOpenMM volume = deltaR[CA][0]*(deltaR[CB][1]*deltaR[CD][2] - deltaR[CB][2]*deltaR[CD][1]) +
                        deltaR[CB][0]*(deltaR[CD][1]*deltaR[CA][2] - deltaR[CD][2]*deltaR[CA][1]) +
                        deltaR[CD][0]*(deltaR[CA][1]*deltaR[CB][2] - deltaR[CA][2]*deltaR[CB][1]);

    return (volume >= zero ? one : -one);
 
 }
 
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

RealOpenMM AmoebaReferenceTorsionTorsionForce::calculateTorsionTorsionIxn(const RealVec& positionAtomA, const RealVec& positionAtomB,
                                                                          const RealVec& positionAtomC, const RealVec& positionAtomD,
                                                                          const RealVec& positionAtomE, const RealVec* positionChiralCheckAtom,
                                                                          const std::vector< std::vector< std::vector<RealOpenMM> > >& grid,
                                                                          RealVec* forces) const {
 
   // ---------------------------------------------------------------------------------------

    //static const std::string methodName = "AmoebaReferenceTorsionTorsionForce::calculateForceAndEnergy";
 
    static const RealOpenMM zero          = 0.0;
    static const RealOpenMM one           = 1.0;

    // ---------------------------------------------------------------------------------------
 
    enum { A, B, C, D, E, LastAtomIndex };
 
    // get deltaR between various combinations of the 4 atoms
    // and various intermediate terms

    enum { BA, CB, DC, ED, T, U, V, UxV, CA, DB, EC, dT, dU, dU2, dV2, LastDeltaIndex };
 
    std::vector<RealOpenMM> deltaR[LastDeltaIndex];
    for (unsigned int ii = 0; ii < LastDeltaIndex; ii++) {
        deltaR[ii].resize(3);
    }   

    if (usePeriodic) {
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomA, positionAtomB, deltaR[BA], boxVectors);
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomB, positionAtomC, deltaR[CB], boxVectors);
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomC, positionAtomD, deltaR[DC], boxVectors);
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomD, positionAtomE, deltaR[ED], boxVectors);
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomA, positionAtomC, deltaR[CA], boxVectors);
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomB, positionAtomD, deltaR[DB], boxVectors);
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomC, positionAtomE, deltaR[EC], boxVectors);
    }
    else {
        AmoebaReferenceForce::loadDeltaR(positionAtomA, positionAtomB, deltaR[BA]);
        AmoebaReferenceForce::loadDeltaR(positionAtomB, positionAtomC, deltaR[CB]);
        AmoebaReferenceForce::loadDeltaR(positionAtomC, positionAtomD, deltaR[DC]);
        AmoebaReferenceForce::loadDeltaR(positionAtomD, positionAtomE, deltaR[ED]);
        AmoebaReferenceForce::loadDeltaR(positionAtomA, positionAtomC, deltaR[CA]);
        AmoebaReferenceForce::loadDeltaR(positionAtomB, positionAtomD, deltaR[DB]);
        AmoebaReferenceForce::loadDeltaR(positionAtomC, positionAtomE, deltaR[EC]);
    }

    std::vector<RealOpenMM> d[LastAtomIndex];
    for (unsigned int ii = 0; ii < LastAtomIndex; ii++) {
        d[ii].resize(3);
    }   
 
    AmoebaReferenceForce::getCrossProduct(deltaR[BA], deltaR[CB], deltaR[T]);
    AmoebaReferenceForce::getCrossProduct(deltaR[CB], deltaR[DC], deltaR[U]);
    AmoebaReferenceForce::getCrossProduct(deltaR[DC], deltaR[ED], deltaR[V]);
 
    AmoebaReferenceForce::getCrossProduct(deltaR[U],  deltaR[V],  deltaR[UxV]);
 
    RealOpenMM rT2  = AmoebaReferenceForce::getNormSquared3(deltaR[T]);
    RealOpenMM rU2  = AmoebaReferenceForce::getNormSquared3(deltaR[U]);
    RealOpenMM rV2  = AmoebaReferenceForce::getNormSquared3(deltaR[V]);
    RealOpenMM rUrV = SQRT(rU2*rV2);
 
    RealOpenMM rTrU = SQRT(rT2*rU2);
 
    if (rTrU <= zero || rUrV <= zero) {
       return zero;
    }
    RealOpenMM rCB         = AmoebaReferenceForce::getNorm3(deltaR[CB]);
    RealOpenMM cosine1     = AmoebaReferenceForce::getDotProduct3(deltaR[T], deltaR[U]);
               cosine1    /= rTrU;
 
    RealOpenMM angle1;
    if (cosine1 <= -one) {
       angle1 = PI_M*RADIAN;
    } else if (cosine1 >= one) {
       angle1 = zero;
    } else {
       angle1 = RADIAN*ACOS(cosine1);
    }
 
    RealOpenMM sign = AmoebaReferenceForce::getDotProduct3(deltaR[BA], deltaR[U]);
    if (sign < zero) {
       angle1 = -angle1; 
    }
 
    // value1 = angle1;
 
    RealOpenMM rDC         = AmoebaReferenceForce::getNorm3(deltaR[DC]);
    RealOpenMM cosine2     = AmoebaReferenceForce::getDotProduct3(deltaR[U], deltaR[V]);
               cosine2    /= rUrV;
 
    RealOpenMM angle2;
    if (cosine2 <= -one) {
       angle2 = PI_M*RADIAN;
    } else if (cosine1 >= one) {
       angle2 = zero;
    } else {
       angle2 = RADIAN*ACOS(cosine2);
    }
 
    sign = AmoebaReferenceForce::getDotProduct3(deltaR[CB], deltaR[V]);
    if (sign < zero) {
       angle2 = -angle2; 
    }
 
    // swap signs of angles if chirality at central atom
    // is 'negative'
  
    if (positionChiralCheckAtom) {
        sign = checkTorsionSign(*positionChiralCheckAtom, positionAtomB, positionAtomC, positionAtomD);
        if (sign < zero) {
           angle1 = -angle1;
           angle2 = -angle2;
        }
    } else {
        sign = one;
    }

    // bicubic interpolation
 
    RealOpenMM corners[2][2];
    RealOpenMM eValues[4][4];
    enum { E0, E1, E2, E12, LastEIndex };
    loadGridValuesFromEnclosingRectangle(grid, angle1, angle2, corners, eValues[E0], eValues[E1], eValues[E2], eValues[E12]);
    
    // get coordinates of point in grid closest to angle1 & angle2
  
    // get corners of grid encompassing point
    // get width/height of encompassing rectangle
 
    RealOpenMM gridEnergy;
    RealOpenMM dEdAngle1;
    RealOpenMM dEdAngle2;
    AmoebaReferenceTorsionTorsionForce::getBicubicValues(
          eValues[E0], eValues[E1], eValues[E2], eValues[E12],
          corners[0][0], corners[0][1], corners[1][0], corners[1][1],
          angle1, angle2, &gridEnergy, &dEdAngle1, &dEdAngle2); 
 
    dEdAngle1 = sign*RADIAN*dEdAngle1;
    dEdAngle2 = sign*RADIAN*dEdAngle2;
 
    AmoebaReferenceForce::getCrossProduct(deltaR[T], deltaR[CB], deltaR[dT]);
    AmoebaReferenceForce::getCrossProduct(deltaR[U], deltaR[CB], deltaR[dU]);
 
    RealOpenMM factorT =  dEdAngle1/(rCB*rT2);
    RealOpenMM factorU = -dEdAngle1/(rCB*rU2);
 
    deltaR[dT][0] *= factorT;
    deltaR[dT][1] *= factorT;
    deltaR[dT][2] *= factorT;

    deltaR[dU][0] *= factorU;
    deltaR[dU][1] *= factorU;
    deltaR[dU][2] *= factorU;
 
    AmoebaReferenceForce::getCrossProduct(deltaR[dT], deltaR[CB], d[A]);
    AmoebaReferenceForce::getCrossProduct(deltaR[dU], deltaR[CB], d[D]);
 
    std::vector<RealOpenMM> tmp[3];
    for (unsigned int ii = 0; ii < 3; ii++) {
        tmp[ii].resize(3);
    }   
    AmoebaReferenceForce::getCrossProduct(deltaR[CA], deltaR[dT], d[B]);
    AmoebaReferenceForce::getCrossProduct(deltaR[dU], deltaR[DC], tmp[0]);
 
    AmoebaReferenceForce::getCrossProduct(deltaR[dT], deltaR[BA], d[C]);
    AmoebaReferenceForce::getCrossProduct(deltaR[DB], deltaR[dU], tmp[1]);

    d[B][0] += tmp[0][0];
    d[B][1] += tmp[0][1];
    d[B][2] += tmp[0][2];

    d[C][0] += tmp[1][0];
    d[C][1] += tmp[1][1];
    d[C][2] += tmp[1][2];
 
    // angle2
 
    AmoebaReferenceForce::getCrossProduct(deltaR[U], deltaR[DC], deltaR[dU2]);
    AmoebaReferenceForce::getCrossProduct(deltaR[V], deltaR[DC], deltaR[dV2]);
 
    RealOpenMM factorU2 =  dEdAngle2/(rDC*rU2);
    RealOpenMM factorV2 = -dEdAngle2/(rDC*rV2);
 
    deltaR[dU2][0] *= factorU2;
    deltaR[dU2][1] *= factorU2;
    deltaR[dU2][2] *= factorU2;

    deltaR[dV2][0] *= factorV2;
    deltaR[dV2][1] *= factorV2;
    deltaR[dV2][2] *= factorV2;
 
    // dB2
  
    AmoebaReferenceForce::getCrossProduct(deltaR[dU2], deltaR[DC],  tmp[0]);
 
    // dC2
  
    AmoebaReferenceForce::getCrossProduct(deltaR[DB],  deltaR[dU2], tmp[1]);
    AmoebaReferenceForce::getCrossProduct(deltaR[dV2], deltaR[ED],  tmp[2]);

    d[B][0] += tmp[0][0];
    d[B][1] += tmp[0][1];
    d[B][2] += tmp[0][2];

    d[C][0] += tmp[1][0] + tmp[2][0];
    d[C][1] += tmp[1][1] + tmp[2][1];
    d[C][2] += tmp[1][2] + tmp[2][2];
 
    // dD2
  
    AmoebaReferenceForce::getCrossProduct(deltaR[dU2],  deltaR[CB], tmp[0]);
    AmoebaReferenceForce::getCrossProduct(deltaR[EC], deltaR[dV2],  tmp[1]);

    d[D][0] += tmp[0][0] + tmp[1][0];
    d[D][1] += tmp[0][1] + tmp[1][1];
    d[D][2] += tmp[0][2] + tmp[1][2];
 
    // dE
  
    AmoebaReferenceForce::getCrossProduct(deltaR[dV2], deltaR[DC],  d[E]);
 
    // ---------------------------------------------------------------------------------------
 
    // add in forces
 
    for (int jj = 0; jj < LastAtomIndex; jj++) {
        forces[jj][0] = d[jj][0];
        forces[jj][1] = d[jj][1];
        forces[jj][2] = d[jj][2];
    }
 
    // ---------------------------------------------------------------------------------------
 
    return gridEnergy;
}

RealOpenMM AmoebaReferenceTorsionTorsionForce::calculateForceAndEnergy(int numTorsionTorsions, vector<RealVec>& posData,
                                                                       const std::vector<int>&  particle1,
                                                                       const std::vector<int>&  particle2,
                                                                       const std::vector<int>&  particle3,
                                                                       const std::vector<int>&  particle4,
                                                                       const std::vector<int>&  particle5,
                                                                       const std::vector<int>&  chiralCheckAtom,
                                                                       const std::vector<int>&  gridIndices,
                                                                       const std::vector< std::vector< std::vector< std::vector<RealOpenMM> > > >& torsionTorsionGrids,
                                                                       vector<RealVec>& forceData) const {
    RealOpenMM energy  = 0.0; 
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(numTorsionTorsions); ii++) {

        int particle1Index       = particle1[ii];
        int particle2Index       = particle2[ii];
        int particle3Index       = particle3[ii];
        int particle4Index       = particle4[ii];
        int particle5Index       = particle5[ii];

        int chiralCheckAtomIndex = chiralCheckAtom[ii];

        int gridIndex            = gridIndices[ii];

        RealVec forces[5];
        RealVec* chiralCheckAtom;
        if (chiralCheckAtomIndex > -1) {
            chiralCheckAtom = &posData[chiralCheckAtomIndex];
        } else {
            chiralCheckAtom = NULL;
        }
        energy                 += calculateTorsionTorsionIxn(posData[particle1Index], posData[particle2Index],
                                                             posData[particle3Index], posData[particle4Index],
                                                             posData[particle5Index], chiralCheckAtom, torsionTorsionGrids[gridIndex],
                                                             forces);

        // accumulate forces
    
        for (int jj = 0; jj < 3; jj++) {
            forceData[particle1Index][jj] -= forces[0][jj];
            forceData[particle2Index][jj] -= forces[1][jj];
            forceData[particle3Index][jj] -= forces[2][jj];
            forceData[particle4Index][jj] -= forces[3][jj];
            forceData[particle5Index][jj] -= forces[4][jj];
        }

    }   
    return energy;
}

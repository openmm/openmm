/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
 * Authors:                                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/AmoebaTorsionTorsionForce.h"
#include "openmm/internal/AmoebaTorsionTorsionForceImpl.h"
#include "openmm/internal/SplineFitter.h"

using namespace OpenMM;
using namespace std;

AmoebaTorsionTorsionForce::AmoebaTorsionTorsionForce() : usePeriodic(false) {
}

int AmoebaTorsionTorsionForce::addTorsionTorsion(int particle1, int particle2, int particle3,
                                                 int particle4, int particle5, int chiralCheckAtomIndex,
                                                 int gridIndex) {
    torsionTorsions.push_back(TorsionTorsionInfo(particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex));
    return torsionTorsions.size()-1;
}

void AmoebaTorsionTorsionForce::getTorsionTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4,
                                                            int& particle5, int& chiralCheckAtomIndex, int& gridIndex) const {
    particle1                = torsionTorsions[index].particle1;
    particle2                = torsionTorsions[index].particle2;
    particle3                = torsionTorsions[index].particle3;
    particle4                = torsionTorsions[index].particle4;
    particle5                = torsionTorsions[index].particle5;
    chiralCheckAtomIndex     = torsionTorsions[index].chiralCheckAtomIndex;
    gridIndex                = torsionTorsions[index].gridIndex;
}

void AmoebaTorsionTorsionForce::setTorsionTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4,
                                                            int particle5, int chiralCheckAtomIndex, int gridIndex) {
    torsionTorsions[index].particle1              = particle1;
    torsionTorsions[index].particle2              = particle2;
    torsionTorsions[index].particle3              = particle3;
    torsionTorsions[index].particle4              = particle4;
    torsionTorsions[index].particle5              = particle5;
    torsionTorsions[index].chiralCheckAtomIndex   = chiralCheckAtomIndex;
    torsionTorsions[index].gridIndex              = gridIndex;
}

const TorsionTorsionGrid& AmoebaTorsionTorsionForce::getTorsionTorsionGrid(int index) const {
   return torsionTorsionGrids[index].getTorsionTorsionGrid();
}

void AmoebaTorsionTorsionForce::setTorsionTorsionGrid(int index, const TorsionTorsionGrid& grid) {
   if (index >= static_cast<int>(torsionTorsionGrids.size())) {
      torsionTorsionGrids.resize(index + 1);
   }
   torsionTorsionGrids[index] = grid;
}

ForceImpl* AmoebaTorsionTorsionForce::createImpl() const {
    return new AmoebaTorsionTorsionForceImpl(*this);
}

AmoebaTorsionTorsionForce::TorsionTorsionGridInfo::TorsionTorsionGridInfo(const TorsionTorsionGrid& grid) {
    if (grid[0][0][0] != grid[1][0][0])
        _grid = grid;
    else {
        // We need to transpose the grid.
        
        int xsize = grid[0].size();
        int ysize = grid.size();
        _grid.resize(xsize);
        for (int i = 0; i < xsize; i++) {
            _grid[i].resize(ysize);
            for (int j = 0; j < ysize; j++)
                _grid[i][j] = grid[j][i];
        }
    }

    _startValues[0] =  _grid[0][0][0];
    _startValues[1] =  _grid[0][0][1];

    _spacing[0]     = static_cast<double>(_grid.size()-1)/360.0;
    _spacing[1]     = static_cast<double>(_grid.size()-1)/360.0;

    _size[0]        = static_cast<int>(_grid.size());
    _size[1]        = static_cast<int>(_grid[0].size());
    
    if (_grid[0][0].size() == 3) {
        // We need to compute the derivatives ourselves.  First determine if the grid is periodic.
        
        int xsize = _size[0];
        int ysize = _size[1];
        bool periodic = true;
        for (int i = 0; i < xsize; i++)
            if (_grid[i][0][2] != _grid[i][ysize-1][2])
                periodic = false;
        for (int i = 0; i < ysize; i++)
            if (_grid[0][i][2] != _grid[xsize-1][i][2])
                periodic = false;
        
        // Compute derivatives with respect to the first angle.

        vector<double> x(xsize), y(ysize);
        for (int i = 0; i < xsize; i++)
            x[i] = _grid[i][0][0];
        for (int i = 0; i < ysize; i++)
            y[i] = _grid[0][i][1];
        vector<double> d1(xsize*ysize), d2(xsize*ysize), d12(xsize*ysize);
        vector<double> t(xsize), deriv(xsize);
        for (int i = 0; i < ysize; i++) {
            for (int j = 0; j < xsize; j++)
                t[j] = _grid[j][i][2];
            if (periodic)
                SplineFitter::createPeriodicSpline(x, t, deriv);
            else
                SplineFitter::createNaturalSpline(x, t, deriv);
            for (int j = 0; j < xsize; j++)
                d1[j+xsize*i] = SplineFitter::evaluateSplineDerivative(x, t, deriv, x[j]);
        }

        // Compute derivatives with respect to the second angle.

        t.resize(ysize);
        deriv.resize(ysize);
        for (int i = 0; i < xsize; i++) {
            for (int j = 0; j < ysize; j++)
                t[j] = _grid[i][j][2];
            if (periodic)
                SplineFitter::createPeriodicSpline(y, t, deriv);
            else
                SplineFitter::createNaturalSpline(y, t, deriv);
            for (int j = 0; j < ysize; j++)
                d2[i+xsize*j] = SplineFitter::evaluateSplineDerivative(y, t, deriv, y[j]);
        }

        // Compute cross derivatives.

        t.resize(xsize);
        deriv.resize(xsize);
        for (int i = 0; i < ysize; i++) {
            for (int j = 0; j < xsize; j++)
                t[j] = d2[j+xsize*i];
            if (periodic)
                SplineFitter::createPeriodicSpline(x, t, deriv);
            else
                SplineFitter::createNaturalSpline(x, t, deriv);
            for (int j = 0; j < xsize; j++)
                d12[j+xsize*i] = SplineFitter::evaluateSplineDerivative(x, t, deriv, x[j]);
        }
        
        // Add the derivatives to the grid.
        
        for (int i = 0; i < xsize; i++)
            for (int j = 0; j < ysize; j++) {
                _grid[i][j].push_back(d1[i+xsize*j]);
                _grid[i][j].push_back(d2[i+xsize*j]);
                _grid[i][j].push_back(d12[i+xsize*j]);
            }
    }
}

void AmoebaTorsionTorsionForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool AmoebaTorsionTorsionForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}

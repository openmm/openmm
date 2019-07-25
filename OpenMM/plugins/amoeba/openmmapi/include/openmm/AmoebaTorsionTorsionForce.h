#ifndef OPENMM_AMOEBA_TORSION_TORSION_FORCE_H_
#define OPENMM_AMOEBA_TORSION_TORSION_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs, Peter Eastman                                    *
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
#include "internal/windowsExportAmoeba.h"

#include <vector>
#include <cmath>

namespace OpenMM {

typedef std::vector< std::vector< std::vector<double> > > TorsionTorsionGrid;
typedef std::vector< std::vector< std::vector<float> > > TorsionTorsionGridFloat;

/**
 * This class implements the Amoeba torsion-torsion interaction.
 *
 * To use it, create an AmoebaTorsionTorsionForce object then call addTorsionTorsion() once for each torsion-torsion.  After
 * a torsion-torsion has been added, you can modify its force field parameters by calling setTorsionTorsionParameters().
 */

class OPENMM_EXPORT_AMOEBA AmoebaTorsionTorsionForce : public Force {

public:
    /**
     * Create an AmoebaTorsionTorsionForce.
     */
    AmoebaTorsionTorsionForce(void);

    /**
     * Get the number of torsion-torsion terms in the potential function
     */
    int getNumTorsionTorsions(void) const {
        return torsionTorsions.size();
    }

    /**
     * Get the number of torsion-torsion grids
     */
    int getNumTorsionTorsionGrids(void) const {
        return torsionTorsionGrids.size();
    }

    /**
     * Add a torsion-torsion term to the force field.
     *
     * @param particle1                 the index of the first particle connected by the torsion-torsion
     * @param particle2                 the index of the second particle connected by the torsion-torsion
     * @param particle3                 the index of the third particle connected by the torsion-torsion
     * @param particle4                 the index of the fourth particle connected by the torsion-torsion
     * @param particle5                 the index of the fifth particle connected by the torsion-torsion
     * @param chiralCheckAtomIndex      the index of the particle connected to particle3, but not particle2 or particle4 to be used in chirality check
     * @param gridIndex                 the index to the grid to be used
     * @return                          the index of the torsion-torsion that was added
     */
    int addTorsionTorsion(int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex);

    /**
     * Get the force field parameters for a torsion-torsion term.
     *
     * @param index                          the index of the torsion-torsion for which to get parameters
     * @param[out] particle1                 the index of the first particle connected by the torsion-torsion
     * @param[out] particle2                 the index of the second particle connected by the torsion-torsion
     * @param[out] particle3                 the index of the third particle connected by the torsion-torsion
     * @param[out] particle4                 the index of the fourth particle connected by the torsion-torsion
     * @param[out] particle5                 the index of the fifth particle connected by the torsion-torsion
     * @param[out] chiralCheckAtomIndex      the index of the particle connected to particle3, but not particle2 or particle4 to be used in chirality check
     * @param[out] gridIndex                 the grid index
     */
    void getTorsionTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, int& particle5, int& chiralCheckAtomIndex, int& gridIndex) const;

    /**
     * Set the force field parameters for a torsion-torsion term.
     *
     * @param index                     the index of the torsion-torsion for which to set parameters
     * @param particle1                 the index of the first particle connected by the torsion-torsion
     * @param particle2                 the index of the second particle connected by the torsion-torsion
     * @param particle3                 the index of the third particle connected by the torsion-torsion
     * @param particle4                 the index of the fourth particle connected by the torsion-torsion
     * @param particle5                 the index of the fifth particle connected by the torsion-torsion
     * @param chiralCheckAtomIndex      the index of the particle connected to particle3, but not particle2 or particle4 to be used in chirality check
     * @param gridIndex                 the grid index
     */
    void setTorsionTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex);

    /**
     * Get the torsion-torsion grid at the specified index
     *
     * @param  index        the grid index
     * @return grid         return grid reference
     */
    const std::vector<std::vector<std::vector<double> > >& getTorsionTorsionGrid(int index) const;

    /**
     * Set the torsion-torsion grid at the specified index
     *
     * @param index         the index of the torsion-torsion for which to get parameters
     * @param grid          either 3 or 6 values may be specified per grid point.  If the derivatives
     *                      are omitted, they are calculated automatically by fitting a 2D spline to
     *                      the energies.
     *                         grid[x][y][0] = x value
     *                         grid[x][y][1] = y value
     *                         grid[x][y][2] = energy
     *                         grid[x][y][3] = dEdx value
     *                         grid[x][y][4] = dEdy value
     *                         grid[x][y][5] = dEd(xy) value
     */
    void setTorsionTorsionGrid(int index, const std::vector<std::vector<std::vector<double> > >& grid);
    /**
     * Set whether this force should apply periodic boundary conditions when calculating displacements.
     * Usually this is not appropriate for bonded forces, but there are situations when it can be useful.
     */
    void setUsesPeriodicBoundaryConditions(bool periodic);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const;
protected:
    ForceImpl* createImpl() const;
private:
    class TorsionTorsionInfo;
    class TorsionTorsionGridInfo;
    std::vector<TorsionTorsionInfo> torsionTorsions;
    std::vector<TorsionTorsionGridInfo> torsionTorsionGrids;
    bool usePeriodic;
};

/**
 * This is an internal class used to record information about a torsion-torsion term.
 * @private
 */
class AmoebaTorsionTorsionForce::TorsionTorsionInfo {

public:

    int particle1, particle2, particle3, particle4, particle5;
    int chiralCheckAtomIndex;
    int gridIndex;
    TorsionTorsionInfo() {
        particle1 = particle2  = particle3 = particle4 = particle5 = chiralCheckAtomIndex = -1;
        gridIndex = 0;
    }
    TorsionTorsionInfo(int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex) :
                       particle1(particle1), particle2(particle2), particle3(particle3),
                       particle4(particle4), particle5(particle5), chiralCheckAtomIndex(chiralCheckAtomIndex), gridIndex(gridIndex) {

    }
};

/**
 * This is an internal class used to record information about a grid.
 * @private
 */
class AmoebaTorsionTorsionForce::TorsionTorsionGridInfo {

public:

    TorsionTorsionGridInfo() {
        _size[0]        = _size[1]        = 0;
        _startValues[0] = _startValues[1] = 0.0;
        _spacing[0]     = _spacing[1]     = 1.0;
    }

    TorsionTorsionGridInfo(const TorsionTorsionGrid& grid);

    const TorsionTorsionGrid& getTorsionTorsionGrid(void) const {
        return _grid;
    }
    int getDimensionSize(int index) const {
        return _size[index];
    }
    double getStartValue(int index) const {
        return _startValues[index];
    }
    double getSpacing(int index) const {
        return _spacing[index];
    }

private:

    TorsionTorsionGrid _grid;
    int _size[2];
    double _startValues[2];
    double _spacing[2];
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_TORSION_TORSION_FORCE_H_*/

#ifndef OPENMM_AMOEBA_TORSION_TORSION_FORCE_H_
#define OPENMM_AMOEBA_TORSION_TORSION_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                              AmoebaOpenMM                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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
#include <vector>
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

typedef std::vector< std::vector< std::vector<double> > > TorsionTorsionGrid;

/**
 * This class implements the Amoeba Out-of-plane bend interaction
 * To use it, create a TorsionTorsionForce object then call addTorsionTorsion() once for each torsionTorsion.  After
 * a torsionTorsion has been added, you can modify its force field parameters by calling setTorsionTorsionParameters().
 */

class OPENMM_EXPORT AmoebaTorsionTorsionForce : public Force {
public:
    /**
     * Create a Amoeba TorsionTorsionForce.
     */
    AmoebaTorsionTorsionForce();

    /**
     * Get the number of torsionTorsion terms in the potential function
     */
    int getNumTorsionTorsions() const {
        return torsionTorsions.size();
    }

    /**
     * Get the number of torsionTorsion grids
     */
    int getNumTorsionTorsionGrids() const {
        return torsionTorsionGrids.size();
    }

    /**
     * Add a torsionTorsion term to the force field.
     *
     * @param particle1                 the index of the first particle connected by the torsionTorsion
     * @param particle2                 the index of the second particle connected by the torsionTorsion
     * @param particle3                 the index of the third particle connected by the torsionTorsion
     * @param particle4                 the index of the fourth particle connected by the torsionTorsion
     * @param particle5                 the index of the fifth particle connected by the torsionTorsion
     * @param chiralCheckAtomIndex      the index of the particle connected to particle3, but not particle2 or particle4 to be used in chirality check
     * @param gridIndex                 the grid index
     * @return             the index of the torsionTorsion that was added
     */
    int addTorsionTorsion(int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex );

    /**
     * Get             the force field parameters for a torsionTorsion term.
     * 
     * @param index                     the index of the torsionTorsion for which to get parameters
     * @param particle1                 the index of the first particle connected by the torsionTorsion
     * @param particle2                 the index of the second particle connected by the torsionTorsion
     * @param particle3                 the index of the third particle connected by the torsionTorsion
     * @param particle4                 the index of the fourth particle connected by the torsionTorsion
     * @param particle5                 the index of the fifth particle connected by the torsionTorsion
     * @param chiralCheckAtomIndex      the index of the particle connected to particle3, but not particle2 or particle4 to be used in chirality check
     * @param gridIndex                 the grid index
     */
    void getTorsionTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, int& particle5, int& chiralCheckAtomIndex, int& gridIndex ) const;

    /**
     * Set             the force field parameters for a torsionTorsion term.
     * 
     * @param index                     the index of the torsionTorsion for which to set parameters
     * @param particle1                 the index of the first particle connected by the torsionTorsion
     * @param particle2                 the index of the second particle connected by the torsionTorsion
     * @param particle3                 the index of the third particle connected by the torsionTorsion
     * @param particle4                 the index of the fourth particle connected by the torsionTorsion
     * @param particle5                 the index of the fifth particle connected by the torsionTorsion
     * @param chiralCheckAtomIndex      the index of the particle connected to particle3, but not particle2 or particle4 to be used in chirality check
     * @param gridIndex                 the grid index
     */
    void setTorsionTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex );

    /**
     * Get the torsion-torsion grid at the specified index
     * 
     * @param gridIndex     the grid index
     * @param grid          the grid
     */
    void getTorsionTorsionGrid(int index, TorsionTorsionGrid& grid ) const;

    /**
     * Set the torsion-torsion grid at the specified index
     * 
     * @param index         the index of the torsionTorsion for which to get parameters
     * @param grid          grid
     *                         grid[x][y][0] = x value 
     *                         grid[x][y][1] = y value 
     *                         grid[x][y][2] = function value 
     *                         grid[x][y][3] = dfdx value 
     *                         grid[x][y][4] = dfdy value 
     *                         grid[x][y][5] = dfd(xy) value 
     */
    void setTorsionTorsionGrid(int index, TorsionTorsionGrid& grid );

protected:
    ForceImpl* createImpl();
private:

    class TorsionTorsionInfo;

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    std::vector<TorsionTorsionInfo> torsionTorsions;
    std::vector<TorsionTorsionGrid> torsionTorsionGrids;
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

class AmoebaTorsionTorsionForce::TorsionTorsionInfo {
public:
    int particle1, particle2, particle3, particle4, particle5;
    int chiralCheckAtomIndex;
    int gridIndex;
    TorsionTorsionInfo() {
        particle1 = particle2  = particle3 = particle4 = particle5 = chiralCheckAtomIndex = -1;
        gridIndex = 0;
    }
    TorsionTorsionInfo(int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex ) :
                       particle1(particle1), particle2(particle2), particle3(particle3),
                       particle4(particle4), particle5(particle5), gridIndex(gridIndex), chiralCheckAtomIndex(chiralCheckAtomIndex) {
     
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_TORSION_TORSION_FORCE_H_*/

#ifndef OPENMM_CMAPTORSIONFORCE_H_
#define OPENMM_CMAPTORSIONFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2016 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

#include "Force.h"
#include "Vec3.h"
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements an interaction between pairs of dihedral angles.  The interaction energy is
 * defined by an "energy correction map" (CMAP), which is simply a set of tabulated energy values
 * on a regular grid of (phi, psi) angles.  Natural cubic spline interpolation is used to compute
 * forces and energies at arbitrary values of the two angles.
 *
 * To use this class, first create one or more energy correction maps by calling addMap().  For each
 * one, you provide an array of energies at uniformly spaced values of the two angles.  Next,
 * add interactions by calling addTorsion().  For each one, you specify the sequence of particles used
 * to calculate each of the two dihedral angles, and the index of the map used to calculate their
 * interaction energy.
 */

class OPENMM_EXPORT CMAPTorsionForce : public Force {
public:
    /**
     * Create a CMAPTorsionForce.
     */
    CMAPTorsionForce();
    /**
     * Get the number of maps that have been defined.
     */
    int getNumMaps() const {
        return maps.size();
    }
    /**
     * Get the number of CMAP torsion terms in the potential function
     */
    int getNumTorsions() const {
        return torsions.size();
    }
    /**
     * Create a new map that can be used for torsion pairs.
     *
     * @param size    the size of the map along each dimension
     * @param energy  the energy values for the map.  This must be of length size*size.
     *                The element energy[i+size*j] contains the energy when the first
     *                torsion angle equals i*2*PI/size and the second torsion angle
     *                equals j*2*PI/size.
     * @return the index of the map that was added
     */
    int addMap(int size, const std::vector<double>& energy);
    /**
     * Get the energy values of a map.
     *
     * @param index        the index of the map for which to get energy values
     * @param[out] size    the size of the map along each dimension
     * @param[out] energy  the energy values for the map.  This must be of length size*size.
     *                     The element energy[i+size*j] contains the energy when the first
     *                     torsion angle equals i*2*PI/size and the second torsion angle
     *                     equals j*2*PI/size.
     */
    void getMapParameters(int index, int& size, std::vector<double>& energy) const;
    /**
     * Set the energy values of a map.
     *
     * @param index   the index of the map for which to set energy values
     * @param size    the size of the map along each dimension
     * @param energy  the energy values for the map.  This must be of length size*size.
     *                The element energy[i+size*j] contains the energy when the first
     *                torsion angle equals i*2*PI/size and the second torsion angle
     *                equals j*2*PI/size.
     */
    void setMapParameters(int index, int size, const std::vector<double>& energy);
    /**
     * Add a CMAP torsion term to the force field.
     *
     * @param map   the index of the map to use for this term
     * @param a1    the index of the first particle forming the first torsion
     * @param a2    the index of the second particle forming the first torsion
     * @param a3    the index of the third particle forming the first torsion
     * @param a4    the index of the fourth particle forming the first torsion
     * @param b1    the index of the first particle forming the second torsion
     * @param b2    the index of the second particle forming the second torsion
     * @param b3    the index of the third particle forming the second torsion
     * @param b4    the index of the fourth particle forming the second torsion
     * @return the index of the torsion that was added
     */
    int addTorsion(int map, int a1, int a2, int a3, int a4, int b1, int b2, int b3, int b4);
    /**
     * Get the force field parameters for a CMAP torsion term.
     *
     * @param index      the index of the torsion for which to get parameters
     * @param[out] map   the index of the map to use for this term
     * @param[out] a1    the index of the first particle forming the first torsion
     * @param[out] a2    the index of the second particle forming the first torsion
     * @param[out] a3    the index of the third particle forming the first torsion
     * @param[out] a4    the index of the fourth particle forming the first torsion
     * @param[out] b1    the index of the first particle forming the second torsion
     * @param[out] b2    the index of the second particle forming the second torsion
     * @param[out] b3    the index of the third particle forming the second torsion
     * @param[out] b4    the index of the fourth particle forming the second torsion
     */
    void getTorsionParameters(int index, int& map, int& a1, int& a2, int& a3, int& a4, int& b1, int& b2, int& b3, int& b4) const;
    /**
     * Set the force field parameters for a CMAP torsion term.
     *
     * @param index the index of the torsion for which to set parameters
     * @param map   the index of the map to use for this term
     * @param a1    the index of the first particle forming the first torsion
     * @param a2    the index of the second particle forming the first torsion
     * @param a3    the index of the third particle forming the first torsion
     * @param a4    the index of the fourth particle forming the first torsion
     * @param b1    the index of the first particle forming the second torsion
     * @param b2    the index of the second particle forming the second torsion
     * @param b3    the index of the third particle forming the second torsion
     * @param b4    the index of the fourth particle forming the second torsion
     */
    void setTorsionParameters(int index, int map, int a1, int a2, int a3, int a4, int b1, int b2, int b3, int b4);
    /**
     * Update the map and torsion parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setMapParameters() and setTorsionParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     *
     * The only information that can be updated with this method is the energy values for a map, and the map index
     * for a torsion.  The size of a map and the set of particles involved in a torsion cannot be changed.  Also,
     * new bonds and torsions cannot be added.
     */
    void updateParametersInContext(Context& context);
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
    class MapInfo;
    class CMAPTorsionInfo;
    std::vector<MapInfo> maps;
    std::vector<CMAPTorsionInfo> torsions;
    bool usePeriodic;
};

/**
 * This is an internal class used to record information about a map.
 * @private
 */
class CMAPTorsionForce::MapInfo {
public:
    int size;
    std::vector<double> energy;
    MapInfo() {
        size = -1;
    }
    MapInfo(int size, const std::vector<double>& energy) :
        size(size), energy(energy) {
    }
};

/**
 * This is an internal class used to record information about a torsion.
 * @private
 */
class CMAPTorsionForce::CMAPTorsionInfo {
public:
    int map, a1, a2, a3, a4, b1, b2, b3, b4;
    CMAPTorsionInfo() {
        map = a1 = a2 = a3 = a4 = b1 = b2 = b3 = b4 = -1;
    }
    CMAPTorsionInfo(int map, int a1, int a2, int a3, int a4, int b1, int b2, int b3, int b4) :
        map(map), a1(a1), a2(a2), a3(a3), a4(a4), b1(b1), b2(b2), b3(b3), b4(b4) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_CMAPTORSIONFORCE_H_*/

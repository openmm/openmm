#ifndef OPENMM_AMOEBA_MULTIPOLE_FORCE_H_
#define OPENMM_AMOEBA_MULTIPOLE_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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
#include "openmm/OpenMMException.h"
#include "internal/windowsExportAmoeba.h"
#include "openmm/Vec3.h"

#include <sstream>
#include <vector>

namespace OpenMM {

/**
 * This class implements the Amoeba multipole interaction.
 *
 * To use it, create an AmoebaMultipoleForce object then call addMultipole() once for each atom.  After
 * an entry has been added, you can modify its force field parameters by calling setMultipoleParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 */

class OPENMM_EXPORT_AMOEBA AmoebaMultipoleForce : public Force {

public:

    enum NonbondedMethod {

        /**
         * No cutoff is applied to nonbonded interactions.  The full set of N^2 interactions is computed exactly.
         * This necessarily means that periodic boundary conditions cannot be used.  This is the default.
         */
        NoCutoff = 0,

        /**
         * Periodic boundary conditions are used, and Particle-Mesh Ewald (PME) summation is used to compute the interaction of each particle
         * with all periodic copies of every other particle.
         */
        PME = 1
    };

    enum PolarizationType {

        /**
         * Mutual polarization
         */
        Mutual = 0,

        /**
         * Direct polarization
         */
        Direct = 1
    };

    enum MultipoleAxisTypes { ZThenX = 0, Bisector = 1, ZBisect = 2, ThreeFold = 3, ZOnly = 4, NoAxisType = 5, LastAxisTypeIndex = 6 };

    enum CovalentType {
                          Covalent12 = 0, Covalent13 = 1, Covalent14 = 2, Covalent15 = 3,
                          PolarizationCovalent11 = 4, PolarizationCovalent12 = 5, PolarizationCovalent13 = 6, PolarizationCovalent14 = 7, CovalentEnd = 8 };

    /**
     * Create an AmoebaMultipoleForce.
     */
    AmoebaMultipoleForce();

    /**
     * Get the number of particles in the potential function
     */
    int getNumMultipoles() const {
        return multipoles.size();
    }

    /**
     * Get the method used for handling long-range nonbonded interactions.
     */
    NonbondedMethod getNonbondedMethod() const;

    /**
     * Set the method used for handling long-range nonbonded interactions.
     */
    void setNonbondedMethod(NonbondedMethod method);

    /**
     * Get polarization type
     */
    PolarizationType getPolarizationType() const;

    /**
     * Set the polarization type
     */
    void setPolarizationType(PolarizationType type);

    /**
     * Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     *
     * @return the cutoff distance, measured in nm
     */
    double getCutoffDistance() const;

    /**
     * Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     *
     * @param distance    the cutoff distance, measured in nm
     */
    void setCutoffDistance(double distance);

    /**
     * Get the Ewald alpha parameter.  If this is 0 (the default), a value is chosen automatically
     * based on the Ewald error tolerance.
     *
     * @return the Ewald alpha parameter
     */
    double getAEwald() const;

    /**
     * Set the Ewald alpha parameter.  If this is 0 (the default), a value is chosen automatically
     * based on the Ewald error tolerance.
     *
     * @param Ewald alpha parameter
     */
    void setAEwald(double aewald);

    /**
     * Get the B-spline order to use for PME charge spreading
     *
     * @return the B-spline order
     */
    int getPmeBSplineOrder() const;

    /**
     * Get the PME grid dimensions.  If Ewald alpha is 0 (the default), this is ignored and grid dimensions
     * are chosen automatically based on the Ewald error tolerance.
     *
     * @return the PME grid dimensions
     */
    void getPmeGridDimensions(std::vector<int>& gridDimension) const;

   /**
     * Set the PME grid dimensions.  If Ewald alpha is 0 (the default), this is ignored and grid dimensions
     * are chosen automatically based on the Ewald error tolerance.
     *
     * @param the PME grid dimensions
     */
    void setPmeGridDimensions(const std::vector<int>& gridDimension);

    /**
     * Add multipole-related info for a particle
     *
     * @param charge               the particle's charge
     * @param molecularDipole      the particle's molecular dipole (vector of size 3)
     * @param molecularQuadrupole  the particle's molecular quadrupole (vector of size 9)
     * @param axisType             the particle's axis type
     * @param multipoleAtomZ       index of first atom used in constructing lab<->molecular frames
     * @param multipoleAtomX       index of second atom used in constructing lab<->molecular frames
     * @param multipoleAtomY       index of second atom used in constructing lab<->molecular frames
     * @param thole                Thole parameter
     * @param dampingFactor        dampingFactor parameter
     * @param polarity             polarity parameter
     *
     * @return the index of the particle that was added
     */
    int addMultipole(double charge, const std::vector<double>& molecularDipole, const std::vector<double>& molecularQuadrupole, int axisType,
                     int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity);

    /**
     * Get the multipole parameters for a particle.
     *
     * @param index                the index of the atom for which to get parameters
     * @param charge               the particle's charge
     * @param molecularDipole      the particle's molecular dipole (vector of size 3)
     * @param molecularQuadrupole  the particle's molecular quadrupole (vector of size 9)
     * @param axisType             the particle's axis type
     * @param multipoleAtomZ       index of first atom used in constructing lab<->molecular frames
     * @param multipoleAtomX       index of second atom used in constructing lab<->molecular frames
     * @param multipoleAtomY       index of second atom used in constructing lab<->molecular frames
     * @param thole                Thole parameter
     * @param dampingFactor        dampingFactor parameter
     * @param polarity             polarity parameter
     */
    void getMultipoleParameters(int index, double& charge, std::vector<double>& molecularDipole, std::vector<double>& molecularQuadrupole,
                                int& axisType, int& multipoleAtomZ, int& multipoleAtomX, int& multipoleAtomY, double& thole, double& dampingFactor, double& polarity) const;

    /**
     * Set the multipole parameters for a particle.
     *
     * @param index                the index of the atom for which to set parameters
     * @param charge               the particle's charge
     * @param molecularDipole      the particle's molecular dipole (vector of size 3)
     * @param molecularQuadrupole  the particle's molecular quadrupole (vector of size 9)
     * @param axisType             the particle's axis type
     * @param multipoleAtomZ       index of first atom used in constructing lab<->molecular frames
     * @param multipoleAtomX       index of second atom used in constructing lab<->molecular frames
     * @param multipoleAtomY       index of second atom used in constructing lab<->molecular frames
     * @param polarity             polarity parameter
     */
    void setMultipoleParameters(int index, double charge, const std::vector<double>& molecularDipole, const std::vector<double>& molecularQuadrupole,
                                int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity);

    /**
     * Set the CovalentMap for an atom
     *
     * @param index                the index of the atom for which to set parameters
     * @param typeId               CovalentTypes type
     * @param covalentAtoms        vector of covalent atoms associated w/ the specfied CovalentType
     */
    void setCovalentMap(int index, CovalentType typeId, const std::vector<int>& covalentAtoms);

    /**
     * Get the CovalentMap for an atom
     *
     * @param index                the index of the atom for which to set parameters
     * @param typeId               CovalentTypes type
     * @param covalentAtoms        output vector of covalent atoms associated w/ the specfied CovalentType
     */
    void getCovalentMap(int index, CovalentType typeId, std::vector<int>& covalentAtoms) const;

    /**
     * Get the CovalentMap for an atom
     *
     * @param index                the index of the atom for which to set parameters
     * @param covalentLists        output vector of covalent lists of atoms
     */
    void getCovalentMaps(int index, std::vector < std::vector<int> >& covalentLists) const;

    /**
     * Get the max number of iterations to be used in calculating the mutual induced dipoles
     *
     * @return max number of iterations
     */
    int getMutualInducedMaxIterations(void) const;

    /**
     * Set the max number of iterations to be used in calculating the mutual induced dipoles
     *
     * @param max number of iterations
     */
    void setMutualInducedMaxIterations(int inputMutualInducedMaxIterations);

    /**
     * Get the target epsilon to be used to test for convergence of iterative method used in calculating the mutual induced dipoles
     *
     * @return target epsilon
     */
    double getMutualInducedTargetEpsilon(void) const;

    /**
     * Set the target epsilon to be used to test for convergence of iterative method used in calculating the mutual induced dipoles
     *
     * @param target epsilon
     */
    void setMutualInducedTargetEpsilon(double inputMutualInducedTargetEpsilon);

    /**
     * Get the error tolerance for Ewald summation.  This corresponds to the fractional error in the forces
     * which is acceptable.  This value is used to select the grid dimensions and separation (alpha)
     * parameter so that the average error level will be less than the tolerance.  There is not a
     * rigorous guarantee that all forces on all atoms will be less than the tolerance, however.
     * 
     * This can be overridden by explicitly setting an alpha parameter and grid dimensions to use.
     */
    double getEwaldErrorTolerance() const;
    /**
     * Get the error tolerance for Ewald summation.  This corresponds to the fractional error in the forces
     * which is acceptable.  This value is used to select the grid dimensions and separation (alpha)
     * parameter so that the average error level will be less than the tolerance.  There is not a
     * rigorous guarantee that all forces on all atoms will be less than the tolerance, however.
     * 
     * This can be overridden by explicitly setting an alpha parameter and grid dimensions to use.
     */
    void setEwaldErrorTolerance(double tol);

    /**
     * Get the induced dipole moments of all particles.
     * 
     * @param context    the Context for which to get the induced dipoles
     * @param dipoles    the induced dipole moment of particle i is stored into the i'th element
     */
    void getInducedDipoles(Context& context, std::vector<Vec3>& dipoles);
    
    /**
     * Get the electrostatic potential.
     *
     * @param inputGrid    input grid points over which the potential is to be evaluated
     * @param context      context
     * @param outputElectrostaticPotential output potential
     */

    void getElectrostaticPotential(const std::vector< Vec3 >& inputGrid,
                                    Context& context, std::vector< double >& outputElectrostaticPotential);

    /**
     * Get the system multipole moments.
     * 
     * This method is most useful for non-periodic systems.  When called for a periodic system, only the
     * <i>lowest nonvanishing moment</i> has a well defined value.  This means that if the system has a net
     * nonzero charge, the dipole and quadrupole moments are not well defined and should be ignored.  If the
     * net charge is zero, the dipole moment is well defined (and really represents a dipole density), but
     * the quadrupole moment is still undefined and should be ignored.
     *
     * @param context      context
     * @param outputMultipoleMoments (charge,
                                      dipole_x, dipole_y, dipole_z,
                                      quadrupole_xx, quadrupole_xy, quadrupole_xz,
                                      quadrupole_yx, quadrupole_yy, quadrupole_yz,
                                      quadrupole_zx, quadrupole_zy, quadrupole_zz)
     */
    void getSystemMultipoleMoments(Context& context, std::vector< double >& outputMultipoleMoments);
    /**
     * Update the multipole parameters in a Context to match those stored in this Force object.  This method
     * provides an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setMultipoleParameters() to modify this object's parameters, then call updateParametersInState() to
     * copy them over to the Context.
     * 
     * This method has several limitations.  The only information it updates is the parameters of multipoles.
     * All other aspects of the Force (the nonbonded method, the cutoff distance, etc.) are unaffected and can only be
     * changed by reinitializing the Context.  Furthermore, this method cannot be used to add new multipoles,
     * only to change the parameters of existing ones.
     */
    void updateParametersInContext(Context& context);

protected:
    ForceImpl* createImpl() const;
private:
    NonbondedMethod nonbondedMethod;
    PolarizationType polarizationType;
    double cutoffDistance;
    double aewald;
    int pmeBSplineOrder;
    std::vector<int> pmeGridDimension;
    int mutualInducedMaxIterations;
    double mutualInducedTargetEpsilon;
    double scalingDistanceCutoff;
    double electricConstant;
    double ewaldErrorTol;
    class MultipoleInfo;
    std::vector<MultipoleInfo> multipoles;
};

/**
 * This is an internal class used to record information about a multipole.
 * @private
 */
class AmoebaMultipoleForce::MultipoleInfo {
public:

    int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
    double charge, thole, dampingFactor, polarity;

    std::vector<double> molecularDipole;
    std::vector<double> molecularQuadrupole;
    std::vector< std::vector<int> > covalentInfo;

    MultipoleInfo() {
        axisType = multipoleAtomZ = multipoleAtomX = multipoleAtomY = -1;
        charge   = thole          = dampingFactor  = 0.0;

        molecularDipole.resize(3);
        molecularQuadrupole.resize(9);

    }

    MultipoleInfo(double charge, const std::vector<double>& inputMolecularDipole, const std::vector<double>& inputMolecularQuadrupole,
                   int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity) :
        charge(charge), axisType(axisType), multipoleAtomZ(multipoleAtomZ), multipoleAtomX(multipoleAtomX), multipoleAtomY(multipoleAtomY),
        thole(thole), dampingFactor(dampingFactor), polarity(polarity) {

       covalentInfo.resize(CovalentEnd);

       molecularDipole.resize(3);
       molecularDipole[0]          = inputMolecularDipole[0];
       molecularDipole[1]          = inputMolecularDipole[1];
       molecularDipole[2]          = inputMolecularDipole[2];

       molecularQuadrupole.resize(9);
       molecularQuadrupole[0]      = inputMolecularQuadrupole[0];
       molecularQuadrupole[1]      = inputMolecularQuadrupole[1];
       molecularQuadrupole[2]      = inputMolecularQuadrupole[2];
       molecularQuadrupole[3]      = inputMolecularQuadrupole[3];
       molecularQuadrupole[4]      = inputMolecularQuadrupole[4];
       molecularQuadrupole[5]      = inputMolecularQuadrupole[5];
       molecularQuadrupole[6]      = inputMolecularQuadrupole[6];
       molecularQuadrupole[7]      = inputMolecularQuadrupole[7];
       molecularQuadrupole[8]      = inputMolecularQuadrupole[8];
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_MULTIPOLE_FORCE_H_*/

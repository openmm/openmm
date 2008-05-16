#ifndef OPENMM_GBSAOBCFORCEFIELD_H_
#define OPENMM_GBSAOBCFORCEFIELD_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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
#include <vector>

namespace OpenMM {

/**
 * This class implements an implicit solvation force using the GBSA-OBC model.
 */

class GBSAOBCForceField : public Force {
public:
    /*
     * Create a GBSAOBCForceField.
     * 
     * @param numAtoms    the number of atoms in the system
     */
    GBSAOBCForceField(int numAtoms);
    /**
     * Get the number of atoms in the system.
     */
    int getNumAtoms() const {
        return atoms.size();
    }
    /**
     * Get the force field parameters for an atom.
     * 
     * @param index          the index of the atom for which to get parameters
     * @param charge         the charge of the atom, measured in units of the proton charge
     * @param radius         the GBSA radius of the atom, measured in nm
     * @param scalingFactor  the OBC scaling factor for the atom
     */
    void getAtomParameters(int index, double& charge, double& radius, double& scalingFactor) const;
    /**
     * Set the force field parameters for an atom.
     * 
     * @param index          the index of the atom for which to set parameters
     * @param charge         the charge of the atom, measured in units of the proton charge
     * @param radius         the GBSA radius of the atom, measured in nm
     * @param scalingFactor  the OBC scaling factor for the atom
     */
    void setAtomParameters(int index, double charge, double radius, double scalingFactor);
    /**
     * Get the dielectric constant for the solvent.
     */
    double getSolventDielectric() const {
        return solventDielectric;
    }
    /**
     * Set the dielectric constant for the solvent.
     */
    void setSolventDielectric(double dielectric) {
        solventDielectric = dielectric;
    }
    /**
     * Get the dielectric constant for the solute.
     */
    double getSoluteDielectric() const {
        return soluteDielectric;
    }
    /**
     * Set the dielectric constant for the solute.
     */
    void setSoluteDielectric(double dielectric) {
        soluteDielectric = dielectric;
    }
protected:
    ForceImpl* createImpl();
private:
    class AtomInfo;
    double solventDielectric, soluteDielectric;
    std::vector<AtomInfo> atoms;
};

class GBSAOBCForceField::AtomInfo {
public:
    double charge, radius, scalingFactor;
    AtomInfo() {
        charge = radius = scalingFactor = 0.0;
    }
};

} // namespace OpenMM

#endif /*OPENMM_GBSAOBCFORCEFIELD_H_*/

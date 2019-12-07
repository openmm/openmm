#ifndef OPENMM_NONBONDEDUTILITIES_H_
#define OPENMM_NONBONDEDUTILITIES_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2019 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

namespace OpenMM {

/**
 * This class provides a generic interface for calculating nonbonded interactions.  Clients only need
 * to provide the code for evaluating a single interaction and the list of parameters it depends on.
 * A complete kernel is then synthesized using an appropriate algorithm to evaluate all interactions on
 * all atoms.  Call addInteraction() to define a nonbonded interaction, and addParameter() to define
 * per-particle parameters that the interaction depends on.
 *
 * During each force or energy evaluation, the following sequence of steps takes place:
 *
 * 1. Data structures (e.g. neighbor lists) are calculated to allow nonbonded interactions to be evaluated
 * quickly.
 *
 * 2. calcForcesAndEnergy() is called on each ForceImpl in the System.
 *
 * 3. Finally, the default interaction kernel is invoked to calculate all interactions that were added
 * to it.
 *
 * This sequence means that the default interaction kernel may depend on quantities that were calculated
 * by ForceImpls during calcForcesAndEnergy().
 */

class NonbondedUtilities {
public:
    virtual ~NonbondedUtilities() {
    }
    /**
     * Get whether a cutoff is being used.
     */
    virtual bool getUseCutoff() = 0;
    /**
     * Get whether periodic boundary conditions are being used.
     */
    virtual bool getUsePeriodic() = 0;
    /**
     * Get the maximum cutoff distance used by any interaction.
     */
    virtual double getMaxCutoffDistance() = 0;
};

} // namespace OpenMM

#endif /*OPENMM_NONBONDEDUTILITIES_H_*/

//
// Created by babaid on 05.10.24.
//
#ifndef OPENMM_EXTERNALPUREMDFORCE_H
#define OPENMM_EXTERNALPUREMDFORCE_H

#include "Force.h"
#include "Vec3.h"
#include <map>
#include <vector>
#include "internal/windowsExport.h"


namespace OpenMM {

    class OPENMM_EXPORT ExternalPuremdForce : public Force {
    public:
    /**
     * Create a ExternalPuremdForce.
     */
      ExternalPuremdForce();
    /**
     * Get the number of harmonic bond stretch terms in the potential function
     */
    int getNumAtoms() const {
        return atoms.size();
    }
    /**
     * Add a bond term to the force field.
     *
     * @param particle the index of the particle
     * @return the index of the bond that was added
     */
    int addAtom(int particle);
    /**
     * Get the bonding atom
     *
     * @param index the index of the atoms
     * @particle the particle index is going to be saved here
     */
    void getAtom(int index, int& particle) const;
    /**
     * Update the per-bond parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setBondParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     *
     * The only information this method updates is the values of per-bond parameters.  The set of particles involved
     * in a bond cannot be changed, nor can new bonds be added.
     */
    void updateParametersInContext(Context& context);
    protected:
    ForceImpl* createImpl() const;
    private:
    class AtomInfo;
    std::vector<AtomInfo> atoms;
    bool usePeriodic;
    mutable int numContexts, firstChangedBond, lastChangedBond;
};

/**
 * This is an internal class used to record information about a bond.
 * @private
 */
class ExternalPuremdForce::AtomInfo {
public:
    int particle;
    AtomInfo() {
        particle = -1;
    }
    AtomInfo(int particle) :
            particle(particle) {
    }
};

} // namespace OpenMM


#endif //OPENMM_EXTERNALPUREMDFORCE_H

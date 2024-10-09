//
// Created by babaid on 05.10.24.
//
#ifndef OPENMM_EXTERNALPUREMDFORCE_H_
#define OPENMM_EXTERNALPUREMDFORCE_H_

#include "Force.h"
#include "Vec3.h"
#include <map>
#include <vector>
#include "internal/windowsExport.h"

using namespace OpenMM;
namespace OpenMM {

/**
 * A class that introduces a puremd qmmm force.
 */
    class OPENMM_EXPORT ExternalPuremdForce : public Force {
    public:
    /**
     * Create a puremd force
     */
      ExternalPuremdForce();
    /**
     * Create a ExternalPuremdForce.
     *
     * @param ffieldFile force field file.
     * @param controlFile control file.
     */
      ExternalPuremdForce(const std::string& ffieldfile, const std::string& controlFile);
    /**
     * Get the number of atoms being simulated by puremd
     *
     * @return the number of atoms
     */
    int getNumAtoms() const {
        return atoms.size();
    }
    /**
     * Gets the filenames used by the force field
     *
     * @param ffieldFile Force field file.
     * @param controlFile Control file.
     */
    void getFileNames(std::string& ffieldFile, std::string& controlFile) const
    {
        ffieldFile = ffield_file;
        controlFile = control_file;
    }
    /**
     * Add a bond term to the force field.
     *
     * @param particle the index of the particle
     * @param symbol symbol of the particle
     * @param isQM is it reactive
     * @return the index of the bond that was added
     */
    int addAtom(int particle, char* symbol, bool isQM);
    /**
     * Get the bonding atom
     *
     * @param index the index of the atoms
     * @param particle the particle index is going to be saved here
     * @param symbol1 symbol of the atom
     * @param symbol2 symbol other
     * @param isQM is it reactive
     */
    void getParticleParameters(int index, int& particle, char& symbol1, char& symbol2, int& isQM) const;
    /**
     * Update the per-bond parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setBondParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     *
     * The only information this method updates is the values of per-bond parameters.  The set of particles involved
     * in a bond cannot be changed, nor can new bonds be added.
     *
     * @param context the context
     */
    void updateParametersInContext(Context& context);
    protected:
    ForceImpl* createImpl() const;
    private:
    class AtomInfo;
    //called for each atom contains index, symbol and if its reactive
    std::vector<AtomInfo> atoms;
    //params for puremd
    std::string ffield_file;
    std::string control_file;
    //unused
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
    char symbol1, symbol2;
    bool isQM;

    AtomInfo() {
        particle = -1;
        symbol1 = symbol2 = '\0';
        isQM = false;
    }
    AtomInfo(int particle, std::string symbol, bool isQM) :
            particle(particle), symbol1(symbol[0]), isQM(isQM) {
        if(symbol.size()>1){ symbol2 = symbol[1];}
        else
        {
          symbol2 = '\0';
        }
    }
};

} // namespace OpenMM


#endif //OPENMM_EXTERNALPUREMDFORCE_H_

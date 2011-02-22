#ifndef OPENMM_AMOEBA_VDW_FORCE_H_
#define OPENMM_AMOEBA_VDW_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
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

/**
 * This class implements an interaction between pairs of particles that varies harmonically with the distance
 * between them.  To use it, create a VdwForce object then call addAngle() once for each angle.  After
 * a angle has been added, you can modify its force field parameters by calling setAngleParameters().
 */

 class OPENMM_EXPORT AmoebaVdwForce : public Force {
public:
    /**
     * Create a Amoeba VdwForce.
     */
    AmoebaVdwForce();

    /**
     * Get the number of particles
     */
    int getNumParticles() const {
        return parameters.size();
    }

    /**
     * Set the force field parameters for a vdw particle.
     * 
     * @param particleIndex   the particle index
     * @param ivIndex         the iv index
     * @param classIndex      the class index into the sig-eps table
     * @param sigma           vdw sigma
     * @param epsilon         vdw epsilon
     * @param reductionFactor the reduction factor 
     */
    void setParticleParameters(int particleIndex, int ivIndex, int classIndex,
                               double sigma, double epsilon, double reductionFactor );

    /**
     * Get the force field parameters for a vdw particle.
     * 
     * @param particleIndex   the particle index
     * @param ivIndex         the iv index
     * @param classIndex      the class index into the sig-eps table
     * @param sigma           vdw sigma
     * @param epsilon         vdw epsilon
     * @param reductionFactor the reduction factor 
     */
    void getParticleParameters(int particleIndex, int& ivIndex, int& classIndex,
                               double& sigma, double& epsilon, double& reductionFactor ) const;


    /**
     * Set the force field parameters for a vdw particle.
     * 
     * @param particleIndex   the particle index
     * @param ivIndex         the iv index
     * @param classIndex      the class index into the sig-eps table
     * @param sigma           vdw sigma
     * @param epsilon         vdw epsilon
     * @param reductionFactor the reduction factor
     * @return index of added particle
     */
    int addParticle(int ivIndex, int classIndex, 
                    double sigma, double epsilon, double reductionFactor );

    /**
     * Set sigma combining rule
     * 
     * @param sigmaCombiningRule   sigma combining rule:  'ARITHMETIC', 'GEOMETRIC'. 'CUBIC-MEAN'
     */
    void setSigmaCombiningRule( const std::string& sigmaCombiningRule );

    /**
     * Get sigma combining rule
     * 
     * @return sigmaCombiningRule   sigma combining rule:  'ARITHMETIC', 'GEOMETRIC'. 'CUBIC-MEAN'
     */
    const std::string& getSigmaCombiningRule( void ) const;

    /**
     * Set epsilon combining rule
     * 
     * @param epsilonCombiningRule   epsilon combining rule:  'ARITHMETIC', 'GEOMETRIC'. 'CUBIC-MEAN'
     */
    void setEpsilonCombiningRule( const std::string& epsilonCombiningRule );

    /**
     * Get epsilon combining rule
     * 
     * @return epsilonCombiningRule   epsilon combining rule:  'ARITHMETIC', 'GEOMETRIC'. 'HARMONIC', 'HHG'
     */
    const std::string& getEpsilonCombiningRule( void ) const;

    /**
     * Set exclusions for specified particle
     * 
     * @param particleIndex particle index
     * @param exclusions output vector of exclusions
     */
    void setParticleExclusions( int particleIndex, const std::vector< int >& exclusions );

    /**
     * Get exclusions for specified particle
     * 
     * @param particleIndex particle index
     * @param exclusions output vector of exclusions
     */
    void getParticleExclusions( int particleIndex, std::vector< int >& exclusions ) const;

    /**
     * Set cutoff
     * 
     * @param cutoff cutoff
     */
    void setCutoff( double cutoff );

    /**
     * Get cutoff
     * 
     * @return cutoff
     */
    double getCutoff( void ) const;

    /**
     * Set flag for using neighbor list for vdw ixn
     * 
     * @param neighboristFlag neighbor list flag
     */
    void setUseNeighborList( int neighborListFlag );

    /**
     * Get neighbor list flag for vdw ixn
     * 
     * @return neighbor list flag
     */
    int getUseNeighborList( void ) const;

    /**
     * Set flag for employing periodic boundary conditions
     * 
     * @param pbcFlag if nonozero, use periodic boundary conditions
     */
    void setPBC( int pbcFlag );

    /**
     * Get periodic boundary conditions flag
     * 
     * @return periodic boundary conditions flag (nonzero -> use PBC)
     */
    int getPBC( void ) const;

protected:
    ForceImpl* createImpl();
private:

    class VdwInfo;
    int usePBC;
    int useNeighborList; double cutoff;
    std::string sigmaCombiningRule;
    std::string epsilonCombiningRule;
    std::vector< std::vector<int> > exclusions;

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    std::vector<VdwInfo> parameters;
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

    std::vector< std::vector< std::vector<double> > > sigEpsTable;
};

class AmoebaVdwForce::VdwInfo {
public:
    int ivIndex, classIndex;
    double reductionFactor, sigma, epsilon, cutoff;
    VdwInfo() {
        ivIndex = classIndex = -1;
        reductionFactor      = 0.0;
        sigma                = 1.0;
        epsilon              = 0.0;
    }
    VdwInfo(int ivIndex, int classIndex, double sigma, double epsilon, double  reductionFactor ) :
        ivIndex(ivIndex), classIndex(classIndex), sigma(sigma), epsilon(epsilon), reductionFactor(reductionFactor)  {
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_VDW_FORCE_H_*/

#ifndef OPENMM_AMOEBA_TORSION_FORCE_H_
#define OPENMM_AMOEBA_TORSION_FORCE_H_

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
#include "openmm/OpenMMException.h"
#include "openmm/internal/windowsExport.h"
#include <vector>
#include <sstream>

namespace OpenMM {

/**
 * This class implements an torsion interaction among four particles.
 * To use it, create a TorsionForce object then call addTorsion() once for each angle.  After
 * a angle has been added, you can modify its force field parameters by calling setTorsionParameters().
 */

class OPENMM_EXPORT AmoebaTorsionForce : public Force {

public:

    static const unsigned int ParametersPerTorsion  = 2;

    /**
     * Create a Amoeba TorsionForce.
     */
    AmoebaTorsionForce();

    /**
     * Get the number of torsion terms in the potential function
     */
    int getNumTorsions() const {
        return torsions.size();
    }

    /**
     * Add a torsion term to the force field.
     *
     * @param particle1     the index of the first particle connected by the torsion
     * @param particle2     the index of the second particle connected by the torsion
     * @param particle3     the index of the third particle connected by the torsion
     * @param particle4     the index of the fourth particle connected by the torsion
     * @param torsion1      the vector of torsion params for first index (amplitude, phase, fold)
     * @param torsion2      the vector of torsion params for second index (amplitude, phase, fold)
     * @param torsion3      the vector of torsion params for third index (amplitude, phase, fold)
     * @return the index of the torsion that was added
     */
    int addTorsion(int particle1, int particle2, int particle3, int particle4,
                   const std::vector<double>& torsion1, const std::vector<double>& torsion2, const std::vector<double>& torsion3 );

    /**
     * Get the force field parameters for a torsion term.
     * 
     * @param index         the index of the torsion for which to get parameters
     * @param particle1     the index of the first particle connected by the torsion
     * @param particle2     the index of the second particle connected by the torsion
     * @param particle3     the index of the third particle connected by the torsion
     * @param particle4     the index of the fourth particle connected by the torsion
     * @param torsion1      the vector of torsion params for first index (amplitude, phase, fold)
     * @param torsion2      the vector of torsion params for second index (amplitude, phase, fold)
     * @param torsion3      the vector of torsion params for third index (amplitude, phase, fold)
     */
    void getTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, 
                              std::vector<double>& torsion1, std::vector<double>& torsion2, std::vector<double>& torsion3 ) const;

    /**
     * Set the force field parameters for a torsion term.
     * 
     * @param index         the index of the torsion for which to set parameters
     * @param particle1     the index of the first particle connected by the torsion
     * @param particle2     the index of the second particle connected by the torsion
     * @param particle3     the index of the third particle connected by the torsion
     * @param particle4     the index of the fourth particle connected by the torsion
     * @param torsion1      the vector of torsion params for first index (amplitude, phase, fold)
     * @param torsion2      the vector of torsion params for second index (amplitude, phase, fold)
     * @param torsion3      the vector of torsion params for third index (amplitude, phase, fold)
     */
    void setTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, 
                              const std::vector<double>& torsion1, const std::vector<double>& torsion2, const std::vector<double>& torsion3 );

protected:
    ForceImpl* createImpl();
private:

    class TorsionInfo;

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    std::vector<TorsionInfo> torsions;
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

class AmoebaTorsionForce::TorsionInfo {
private:

    // max number of torsion sets of parameters; in current amoebapro.prm max is 3
    // but code goes up to 6

    static const unsigned int maxTorsions           = 3;

    void _initialize() {

        torsionParameters.resize( maxTorsions );     
        for( unsigned int ii = 0; ii < torsionParameters.size(); ii++ ){
            torsionParameters[ii].resize( AmoebaTorsionForce::ParametersPerTorsion );
            for( unsigned int jj = 0; jj < AmoebaTorsionForce::ParametersPerTorsion; jj++ ){
                torsionParameters[ii][jj] = 0.0;
            }
        }
    }

public:
    int particle1, particle2, particle3, particle4;
    std::vector< std::vector<double> > torsionParameters;
    TorsionInfo() {
        particle1 = particle2  = particle3 = particle4 = -1;
        _initialize();
    }
    TorsionInfo(int particle1, int particle2, int particle3, int particle4, const std::vector<double>& torsion1, const std::vector<double>& torsion2, const std::vector<double>& torsion3 ) :
        particle1(particle1), particle2(particle2), particle3(particle3), particle4(particle4)  {

        if( torsion1.size() != AmoebaTorsionForce::ParametersPerTorsion ){
           std::stringstream buffer;
           buffer << "TorsionInfo::TorsionInfo: torsion1 size(=" << torsion1.size() << ") is not " << AmoebaTorsionForce::ParametersPerTorsion;
           throw OpenMMException( buffer.str() );
        }

        if( torsion2.size() != AmoebaTorsionForce::ParametersPerTorsion ){
           std::stringstream buffer;
           buffer << "TorsionInfo::TorsionInfo: torsion2 size(=" << torsion2.size() << ") is not " << AmoebaTorsionForce::ParametersPerTorsion;
           throw OpenMMException( buffer.str() );
        }

        if( torsion3.size() != AmoebaTorsionForce::ParametersPerTorsion ){
           std::stringstream buffer;
           buffer << "TorsionInfo::TorsionInfo: torsion3 size(=" << torsion3.size() << ") is not " << AmoebaTorsionForce::ParametersPerTorsion;
           throw OpenMMException( buffer.str() );
        }

       _initialize();   
       for( unsigned int ii = 0; ii < AmoebaTorsionForce::ParametersPerTorsion; ii++ ){
           torsionParameters[0][ii] = torsion1[ii];
           torsionParameters[1][ii] = torsion2[ii];
           torsionParameters[2][ii] = torsion3[ii];
       }
    }
    int copyTorsionParameter(int index, const std::vector<double>& torsionParameter ) {

        if( torsionParameter.size() != AmoebaTorsionForce::ParametersPerTorsion ){
           std::stringstream buffer;
           buffer << "TorsionInfo::copyTorsionParameter: input torsionParameter size(=" << torsionParameter.size() << ") is not " << AmoebaTorsionForce::ParametersPerTorsion;
           throw OpenMMException( buffer.str() );
        }

        if( index >= 0 && index < maxTorsions ){
            for( unsigned int ii = 0; ii < AmoebaTorsionForce::ParametersPerTorsion; ii++ ){
                torsionParameters[index][ii] = torsionParameter[ii];
            }
        } else {
           std::stringstream buffer;
           buffer << "TorsionInfo::copyTorsionParameter: input torsionParameter index(=" << index << ") is not in range [0," << maxTorsions << ")!";
           throw OpenMMException( buffer.str() );
        }
        return 0;
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_TORSION_FORCE_H_*/

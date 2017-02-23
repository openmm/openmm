
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef __AmoebaReferenceGeneralizedKirkwoodForce_H__
#define __AmoebaReferenceGeneralizedKirkwoodForce_H__

#include "openmm/Vec3.h"
#include <vector>

using namespace OpenMM;
using namespace std;

namespace OpenMM {

class AmoebaReferenceGeneralizedKirkwoodForce {

public:

    /**
     *  Constructor
     *  
     */
    AmoebaReferenceGeneralizedKirkwoodForce();
 
    /**
     *  Destructor
     *  
     */
    ~AmoebaReferenceGeneralizedKirkwoodForce() {};
 
    /**
     *  Get number of particles 
     *
     *  @return numParticles
     *
     */
    int getNumParticles() const;

    /**
     *  Set numParticles
     *
     *  @param numParticles
     *
     */
    void setNumParticles(int numParticles);

    /**
     * Get includeCavityTerm flag 
     *
     * @return includeCavityTerm
     *
     */
    int getIncludeCavityTerm() const;

    /**
     *  Set includeCavityTerm flag
     *
     *  @param includeCavityTerm flag indicating whether surface area term is to be included
     *
     */
    void setIncludeCavityTerm(int includeCavityTerm);

    /**
     *  Get directPolarization flag 
     *
     *  @return directPolarization
     *
     */
    int getDirectPolarization() const;

    /**
     *  Set directPolarization flag
     * 
     *  @param directPolarization nonzero if direct as opposed to mutual polarization
     */
    void setDirectPolarization(int directPolarization);

    /**
     *  Get solute dielectric
     *
     *  @return soluteDielectric
     */
    double getSoluteDielectric() const;

    /**
     *  Set solute dielectric
     * 
     *  @param soluteDielectric solute dielectric
     * 
     */
    void setSoluteDielectric(double soluteDielectric);

    /**
     *  Get solvent dielectric
     *
     *  @return solventDielectric
     *
     */
    double getSolventDielectric() const;

    /**
     *  Set solvent dielectric 
     *
     *  @param solventDielectric solvent dielectric
     *
     */
    void setSolventDielectric(double solventDielectric);

    /**
     *  Get dielectric offset
     *
     *  @return dielectricOffset
     *
     */
    double getDielectricOffset() const;

    /**
     * Set dielectric offset
     *
     * @param dielectricOffset dielectric offset
     *
     */
    void setDielectricOffset(double dielectricOffset);

    /**
     * Get probeRadius
     *
     * @return probeRadius
     *
     */
    double getProbeRadius() const;

    /**
     * Set probe radius
     *
     * @param probeRadius probe radiue
     *
     */
    void setProbeRadius(double probeRadius);

    /**
     * Get surfaceAreaFactor
     *
     * @return surfaceAreaFactor
     *
     */
    double getSurfaceAreaFactor() const;

    /**
     * Set surface area factor
     *
     * @param surfaceAreaFactor surface area factor
     *
     */
    void setSurfaceAreaFactor(double surfaceAreaFactor);

    /**
     * Set atomic radii
     *
     * @param atomicRadii input vector of atomic radii
     *
     */
    void setAtomicRadii(const vector<double>& atomicRadii);

    /**
     * Get atomic radii
     *
     * @param atomicRadii output vector of atomic radii
     *
     */
    void getAtomicRadii(vector<double>& atomicRadii) const;

    /**
     * Set scale factors
     *
     * @param scaleFactors input vector of scale factors
     *
     */
    void setScaleFactors(const vector<double>& scaleFactors);

    /**
     * Get scale factors
     *
     * @param scaleFactors output vector of scale factors
     *
     */
    void getScaleFactors(vector<double>& scaleFactors) const;

    /**
     * Set charges
     *
     * @param charges input vector of charges
     *
     */
    void setCharges(const vector<double>& charges);

    /**
     * Calculate Grycuk Born radii
     *
     * @param particlePositions particle positions
     *
     */
    void calculateGrycukBornRadii(const vector<Vec3>& particlePositions);
         
    /**
     * Get Grycik Born radii (must have called calculateGrycukBornRadii())
     *
     * @param bornRadii vector of Born radii
     *
     */
    void getGrycukBornRadii(vector<double>& bornRadii) const;     

private:

    int _numParticles;
    int _includeCavityTerm;
    int _directPolarization;

    double _soluteDielectric;
    double _solventDielectric;
    double _dielectricOffset;
    double _probeRadius;
    double _surfaceAreaFactor;

    std::vector<double> _atomicRadii;
    std::vector<double> _scaleFactors;
    std::vector<double> _charges;

    std::vector<double> _bornRadii;

};

} // namespace OpenMM

#endif // _AmoebaReferenceGeneralizedKirkwoodForce___

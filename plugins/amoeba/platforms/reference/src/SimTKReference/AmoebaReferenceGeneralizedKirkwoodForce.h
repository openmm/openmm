
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
     * Get tanhRescaling flag
     *
     * @return tanhRescaling
     */
    int getTanhRescaling() const;

    /**
     * Get Tanh parameters beta0, beta1 and beta2.
    */
    void getTanhParameters(double& b0, double& b1, double& b2) const;

    /**
     * Set the flag signaling whether the solute integral is rescaled by a Tanh function
     * to account for interstitial spaces.
    */
    void setTanhParameters(double b0, double b1, double b2);

    /**
     *  Set tanhRescaling flag
     *
     *  @param tanhRescaling flag indicating whether tanh rescaling of the solute integral is performed.
     */
    void setTanhRescaling(int tanhRescaling);

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
     * Set descreen radii
     *
     * @param descreenRadii input vector of descreen radii
     *
    */
    void setDescreenRadii(const vector<double>& descreenRadii);

    /**
     * Get descreen radii
     *
     * @param descreenRadii output vector of descreen radii
     *
     */
    void getDescreenRadii(vector<double>& descreenRadii) const;

    /**
     * Set neck scale factors
     *
     * @param neckFactors input vector of neck scale factors
     */
    void setNeckFactors(const vector<double>& neckFactors);

    /**
     * Get neck scale factors
     *
     * @param neckFactors output vector of neck scale factors
     */
    void getNeckFactors(vector<double>& neckFactors) const;

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

    /**
     * Get Grycik Solute Integral (must have called calculateGrycukBornRadii())
     *
     * @param bornRadii vector of Born radii
     *
     */
    void getSoluteIntegral(vector<double>& soluteIntegral) const;

    static void getNeckConstants(double rhoDescreened, double rhoDescreening, double constants[]);

private:

    int _numParticles;
    int _includeCavityTerm;
    int _directPolarization;
    int _tanhRescaling;

    double _soluteDielectric;
    double _solventDielectric;
    double _dielectricOffset;
    double _probeRadius;
    double _surfaceAreaFactor;

    double _beta0;
    double _beta1;
    double _beta2;

    std::vector<double> _atomicRadii;
    std::vector<double> _scaleFactors;
    std::vector<double> _charges;
    std::vector<double> _descreenRadii;
    std::vector<double> _neckFactors;

    std::vector<double> _soluteIntegral;
    std::vector<double> _bornRadii;

};

} // namespace OpenMM

#endif // _AmoebaReferenceGeneralizedKirkwoodForce___

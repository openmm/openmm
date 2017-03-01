
/* Portions copyright (c) 2006-2009 Stanford University and Simbios.
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

#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <cstdio>

#include "ReferenceForce.h"
#include "ReferenceObc.h"

using namespace OpenMM;
using namespace std;

/**---------------------------------------------------------------------------------------

    ReferenceObc constructor

    obcParameters      obcParameters object
    
    --------------------------------------------------------------------------------------- */

ReferenceObc::ReferenceObc(ObcParameters* obcParameters) : _obcParameters(obcParameters), _includeAceApproximation(1) {
    _obcChain.resize(_obcParameters->getNumberOfAtoms());
}

/**---------------------------------------------------------------------------------------

    ReferenceObc destructor

    --------------------------------------------------------------------------------------- */

ReferenceObc::~ReferenceObc() {
}

/**---------------------------------------------------------------------------------------

    Get ObcParameters reference

    @return ObcParameters reference

    --------------------------------------------------------------------------------------- */

ObcParameters* ReferenceObc::getObcParameters() const {
    return _obcParameters;
}

/**---------------------------------------------------------------------------------------

    Set ObcParameters reference

    @param ObcParameters reference

    --------------------------------------------------------------------------------------- */

void ReferenceObc::setObcParameters( ObcParameters* obcParameters) {
    _obcParameters = obcParameters;
}

/**---------------------------------------------------------------------------------------

   Return flag signalling whether AceApproximation for nonpolar term is to be included

   @return flag

   --------------------------------------------------------------------------------------- */

int ReferenceObc::includeAceApproximation() const {
    return _includeAceApproximation;
}

/**---------------------------------------------------------------------------------------

   Set flag indicating whether AceApproximation is to be included

   @param includeAceApproximation new includeAceApproximation value

   --------------------------------------------------------------------------------------- */

void ReferenceObc::setIncludeAceApproximation(int includeAceApproximation) {
    _includeAceApproximation = includeAceApproximation;
}

/**---------------------------------------------------------------------------------------

    Return OBC chain derivative: size = _obcParameters->getNumberOfAtoms()

    @return array

    --------------------------------------------------------------------------------------- */

vector<double>& ReferenceObc::getObcChain() {
    return _obcChain;
}

/**---------------------------------------------------------------------------------------

    Get Born radii based on papers:

       J. Phys. Chem. 1996 100, 19824-19839 (HCT paper)
       Proteins: Structure, Function, and Bioinformatcis 55:383-394 (2004) (OBC paper)

    @param atomCoordinates     atomic coordinates
    @param bornRadii           output array of Born radii

    --------------------------------------------------------------------------------------- */

void ReferenceObc::computeBornRadii(const vector<Vec3>& atomCoordinates, vector<double>& bornRadii) {

    ObcParameters* obcParameters              = getObcParameters();

    int numberOfAtoms                         = obcParameters->getNumberOfAtoms();
    const vector<double>& atomicRadii         = obcParameters->getAtomicRadii();
    const vector<double>& scaledRadiusFactor  = obcParameters->getScaledRadiusFactors();
    vector<double>& obcChain                  = getObcChain();

    double dielectricOffset                 = obcParameters->getDielectricOffset();
    double alphaObc                         = obcParameters->getAlphaObc();
    double betaObc                          = obcParameters->getBetaObc();
    double gammaObc                         = obcParameters->getGammaObc();

    // ---------------------------------------------------------------------------------------

    // calculate Born radii

    for (int atomI = 0; atomI < numberOfAtoms; atomI++) {
      
       double radiusI         = atomicRadii[atomI];
       double offsetRadiusI   = radiusI - dielectricOffset;

       double radiusIInverse  = 1.0/offsetRadiusI;
       double sum             = 0.0;

       // HCT code

       for (int atomJ = 0; atomJ < numberOfAtoms; atomJ++) {

          if (atomJ != atomI) {

             double deltaR[ReferenceForce::LastDeltaRIndex];
             if (_obcParameters->getPeriodic())
                 ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomI], atomCoordinates[atomJ], _obcParameters->getPeriodicBox(), deltaR);
             else
                 ReferenceForce::getDeltaR(atomCoordinates[atomI], atomCoordinates[atomJ], deltaR);
             double r               = deltaR[ReferenceForce::RIndex];
             if (_obcParameters->getUseCutoff() && r > _obcParameters->getCutoffDistance())
                 continue;

             double offsetRadiusJ   = atomicRadii[atomJ] - dielectricOffset; 
             double scaledRadiusJ   = offsetRadiusJ*scaledRadiusFactor[atomJ];
             double rScaledRadiusJ  = r + scaledRadiusJ;

             if (offsetRadiusI < rScaledRadiusJ) {
                double rInverse = 1.0/r;
                double l_ij     = offsetRadiusI > fabs(r - scaledRadiusJ) ? offsetRadiusI : fabs(r - scaledRadiusJ);
                       l_ij     = 1.0/l_ij;

                double u_ij     = 1.0/rScaledRadiusJ;

                double l_ij2    = l_ij*l_ij;
                double u_ij2    = u_ij*u_ij;
 
                double ratio    = log((u_ij/l_ij));
                double term     = l_ij - u_ij + 0.25*r*(u_ij2 - l_ij2)  + (0.5*rInverse*ratio) + (0.25*scaledRadiusJ*scaledRadiusJ*rInverse)*(l_ij2 - u_ij2);

                // this case (atom i completely inside atom j) is not considered in the original paper
                // Jay Ponder and the authors of Tinker recognized this and
                // worked out the details

                if (offsetRadiusI < (scaledRadiusJ - r)) {
                   term += 2.0*(radiusIInverse - l_ij);
                }
                sum += term;

             }
          }
       }
 
       // OBC-specific code (Eqs. 6-8 in paper)

       sum              *= 0.5*offsetRadiusI;
       double sum2       = sum*sum;
       double sum3       = sum*sum2;
       double tanhSum    = tanh(alphaObc*sum - betaObc*sum2 + gammaObc*sum3);
       
       bornRadii[atomI]      = 1.0/(1.0/offsetRadiusI - tanhSum/radiusI); 
 
       obcChain[atomI]       = offsetRadiusI*(alphaObc - 2.0*betaObc*sum + 3.0*gammaObc*sum2);
       obcChain[atomI]       = (1.0 - tanhSum*tanhSum)*obcChain[atomI]/radiusI;

    }
}

/**---------------------------------------------------------------------------------------

    Get nonpolar solvation force constribution via ACE approximation

    @param obcParameters parameters
    @param bornRadii                 Born radii
    @param energy                    energy (output): value is incremented from input value 
    @param forces                    forces: values are incremented from input values

    --------------------------------------------------------------------------------------- */

void ReferenceObc::computeAceNonPolarForce(const ObcParameters* obcParameters,
                                      const vector<double>& bornRadii,
                                      double* energy,
                                      vector<double>& forces) const {

    // compute the nonpolar solvation via ACE approximation

    const double probeRadius          = obcParameters->getProbeRadius();
    const double surfaceAreaFactor    = obcParameters->getPi4Asolv();

    const vector<double>& atomicRadii   = obcParameters->getAtomicRadii();
    int numberOfAtoms                     = obcParameters->getNumberOfAtoms();

    // the original ACE equation is based on Eq.2 of

    // M. Schaefer, C. Bartels and M. Karplus, "Solution Conformations
    // and Thermodynamics of Structured Peptides: Molecular Dynamics
    // Simulation with an Implicit Solvation Model", J. Mol. Biol.,
    // 284, 835-848 (1998)  (ACE Method)

    // The original equation includes the factor (atomicRadii[atomI]/bornRadii[atomI]) to the first power,
    // whereas here the ratio is raised to the sixth power: (atomicRadii[atomI]/bornRadii[atomI])**6

    // This modification was made by Jay Ponder who observed it gave better correlations w/
    // observed values. He did not think it was important enough to write up, so there is
    // no paper to cite.

    for (int atomI = 0; atomI < numberOfAtoms; atomI++) {
        if (bornRadii[atomI] > 0.0) {
            double r            = atomicRadii[atomI] + probeRadius;
            double ratio6       = pow(atomicRadii[atomI]/bornRadii[atomI], 6.0);
            double saTerm       = surfaceAreaFactor*r*r*ratio6;
            *energy                += saTerm;
            forces[atomI]          -= 6.0*saTerm/bornRadii[atomI]; 
        }
    }
}

/**---------------------------------------------------------------------------------------

    Get Obc Born energy and forces

    @param atomCoordinates     atomic coordinates
    @param partialCharges      partial charges
    @param forces              forces

    The array bornRadii is also updated and the obcEnergy

    --------------------------------------------------------------------------------------- */

double ReferenceObc::computeBornEnergyForces(const vector<Vec3>& atomCoordinates,
                                           const vector<double>& partialCharges, vector<Vec3>& inputForces) {

    // constants

    const int numberOfAtoms = _obcParameters->getNumberOfAtoms();
    const double dielectricOffset = _obcParameters->getDielectricOffset();
    const double cutoffDistance = _obcParameters->getCutoffDistance();
    const double soluteDielectric = _obcParameters->getSoluteDielectric();
    const double solventDielectric = _obcParameters->getSolventDielectric();
    double preFactor;
    if (soluteDielectric != 0.0 && solventDielectric != 0.0)
        preFactor = 2.0*_obcParameters->getElectricConstant()*((1.0/soluteDielectric) - (1.0/solventDielectric));
    else
        preFactor = 0.0;

    // ---------------------------------------------------------------------------------------

    // compute Born radii

    vector<double> bornRadii(numberOfAtoms);
    computeBornRadii(atomCoordinates, bornRadii);

    // set energy/forces to zero

    double obcEnergy = 0.0;
    vector<double> bornForces(numberOfAtoms, 0.0);

    // ---------------------------------------------------------------------------------------

    // compute the nonpolar solvation via ACE approximation
     
    if (includeAceApproximation()) {
       computeAceNonPolarForce(_obcParameters, bornRadii, &obcEnergy, bornForces);
    }
 
    // ---------------------------------------------------------------------------------------

    // first main loop

    for (int atomI = 0; atomI < numberOfAtoms; atomI++) {
 
       double partialChargeI = preFactor*partialCharges[atomI];
       for (int atomJ = atomI; atomJ < numberOfAtoms; atomJ++) {

          double deltaR[ReferenceForce::LastDeltaRIndex];
          if (_obcParameters->getPeriodic())
              ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomI], atomCoordinates[atomJ], _obcParameters->getPeriodicBox(), deltaR);
          else
              ReferenceForce::getDeltaR(atomCoordinates[atomI], atomCoordinates[atomJ], deltaR);
          if (_obcParameters->getUseCutoff() && deltaR[ReferenceForce::RIndex] > cutoffDistance)
              continue;

          double r2                 = deltaR[ReferenceForce::R2Index];
          double deltaX             = deltaR[ReferenceForce::XIndex];
          double deltaY             = deltaR[ReferenceForce::YIndex];
          double deltaZ             = deltaR[ReferenceForce::ZIndex];

          double alpha2_ij          = bornRadii[atomI]*bornRadii[atomJ];
          double D_ij               = r2/(4.0*alpha2_ij);

          double expTerm            = exp(-D_ij);
          double denominator2       = r2 + alpha2_ij*expTerm; 
          double denominator        = sqrt(denominator2); 
          
          double Gpol               = (partialChargeI*partialCharges[atomJ])/denominator; 
          double dGpol_dr           = -Gpol*(1.0 - 0.25*expTerm)/denominator2;  

          double dGpol_dalpha2_ij   = -0.5*Gpol*expTerm*(1.0 + D_ij)/denominator2;
          
          double energy = Gpol;

          if (atomI != atomJ) {

              if (_obcParameters->getUseCutoff())
                  energy -= partialChargeI*partialCharges[atomJ]/cutoffDistance;
              
              bornForces[atomJ]        += dGpol_dalpha2_ij*bornRadii[atomI];

              deltaX                   *= dGpol_dr;
              deltaY                   *= dGpol_dr;
              deltaZ                   *= dGpol_dr;

              inputForces[atomI][0]    += deltaX;
              inputForces[atomI][1]    += deltaY;
              inputForces[atomI][2]    += deltaZ;

              inputForces[atomJ][0]    -= deltaX;
              inputForces[atomJ][1]    -= deltaY;
              inputForces[atomJ][2]    -= deltaZ;

          } else {
             energy *= 0.5;
          }

          obcEnergy         += energy;
          bornForces[atomI] += dGpol_dalpha2_ij*bornRadii[atomJ];

       }
    }

    // ---------------------------------------------------------------------------------------

    // second main loop

    const vector<double>& obcChain            = getObcChain();
    const vector<double>& atomicRadii         = _obcParameters->getAtomicRadii();

    const double alphaObc                   = _obcParameters->getAlphaObc();
    const double betaObc                    = _obcParameters->getBetaObc();
    const double gammaObc                   = _obcParameters->getGammaObc();
    const vector<double>& scaledRadiusFactor  = _obcParameters->getScaledRadiusFactors();

    // compute factor that depends only on the outer loop index

    for (int atomI = 0; atomI < numberOfAtoms; atomI++) {
       bornForces[atomI] *= bornRadii[atomI]*bornRadii[atomI]*obcChain[atomI];      
    }

    for (int atomI = 0; atomI < numberOfAtoms; atomI++) {
 
       // radius w/ dielectric offset applied

       double radiusI        = atomicRadii[atomI];
       double offsetRadiusI  = radiusI - dielectricOffset;

       for (int atomJ = 0; atomJ < numberOfAtoms; atomJ++) {

          if (atomJ != atomI) {

             double deltaR[ReferenceForce::LastDeltaRIndex];
             if (_obcParameters->getPeriodic())
                ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomI], atomCoordinates[atomJ], _obcParameters->getPeriodicBox(), deltaR);
             else 
                ReferenceForce::getDeltaR(atomCoordinates[atomI], atomCoordinates[atomJ], deltaR);
             if (_obcParameters->getUseCutoff() && deltaR[ReferenceForce::RIndex] > cutoffDistance)
                    continue;
    
             double deltaX             = deltaR[ReferenceForce::XIndex];
             double deltaY             = deltaR[ReferenceForce::YIndex];
             double deltaZ             = deltaR[ReferenceForce::ZIndex];
             double r                  = deltaR[ReferenceForce::RIndex];
 
             // radius w/ dielectric offset applied

             double offsetRadiusJ      = atomicRadii[atomJ] - dielectricOffset;

             double scaledRadiusJ      = offsetRadiusJ*scaledRadiusFactor[atomJ];
             double scaledRadiusJ2     = scaledRadiusJ*scaledRadiusJ;
             double rScaledRadiusJ     = r + scaledRadiusJ;

             // dL/dr & dU/dr are zero (this can be shown analytically)
             // removed from calculation

             if (offsetRadiusI < rScaledRadiusJ) {

                double l_ij          = offsetRadiusI > fabs(r - scaledRadiusJ) ? offsetRadiusI : fabs(r - scaledRadiusJ);
                       l_ij          = 1.0/l_ij;

                double u_ij          = 1.0/rScaledRadiusJ;

                double l_ij2         = l_ij*l_ij;

                double u_ij2         = u_ij*u_ij;
 
                double rInverse      = 1.0/r;
                double r2Inverse     = rInverse*rInverse;

                double t3            = 0.125*(1.0 + scaledRadiusJ2*r2Inverse)*(l_ij2 - u_ij2) + 0.25*log(u_ij/l_ij)*r2Inverse;

                double de            = bornForces[atomI]*t3*rInverse;

                deltaX                  *= de;
                deltaY                  *= de;
                deltaZ                  *= de;
    
                inputForces[atomI][0]   -= deltaX;
                inputForces[atomI][1]   -= deltaY;
                inputForces[atomI][2]   -= deltaZ;
  
                inputForces[atomJ][0]   += deltaX;
                inputForces[atomJ][1]   += deltaY;
                inputForces[atomJ][2]   += deltaZ;
 
             }
          }
       }

    }

    return obcEnergy;
}

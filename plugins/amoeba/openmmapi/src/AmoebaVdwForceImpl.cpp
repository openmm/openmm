/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/AmoebaVdwForceImpl.h"
#include "openmm/amoebaKernels.h"
#include <map>
#include <cmath>

using namespace OpenMM;
using namespace std;

AmoebaVdwForceImpl::AmoebaVdwForceImpl(const AmoebaVdwForce& owner) : owner(owner) {
}

AmoebaVdwForceImpl::~AmoebaVdwForceImpl() {
}

void AmoebaVdwForceImpl::initialize(ContextImpl& context) {
    const System& system = context.getSystem();

    if (owner.getNumParticles() != system.getNumParticles())
        throw OpenMMException("AmoebaVdwForce must have exactly as many particles as the System it belongs to.");

    // check that cutoff < 0.5*boxSize

    if (owner.getNonbondedMethod() == AmoebaVdwForce::CutoffPeriodic) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoffDistance();
        if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
            throw OpenMMException("AmoebaVdwForce: The cutoff distance cannot be greater than half the periodic box size.");
    }   

    kernel = context.getPlatform().createKernel(CalcAmoebaVdwForceKernel::Name(), context);
    kernel.getAs<CalcAmoebaVdwForceKernel>().initialize(context.getSystem(), owner);
}

double AmoebaVdwForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcAmoebaVdwForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

void AmoebaVdwForceImpl::createParameterMatrix(const AmoebaVdwForce& force, vector<int>& type,
        vector<vector<double> >& sigmaMatrix, vector<vector<double> >& epsilonMatrix) {
    int numParticles = force.getNumParticles();
    type.resize(numParticles);
    int numTypes;
    vector<double> typeSigma, typeEpsilon;
    if (force.getUseParticleTypes()) {
        // We get the types directly from the particles.

        double sigma, epsilon, reduction;
        int parent;
        bool isAlchemical;
        for (int i = 0; i < numParticles; i++)
            force.getParticleParameters(i, parent, sigma, epsilon, reduction, isAlchemical, type[i]);
        numTypes = force.getNumParticleTypes();
        typeSigma.resize(numTypes);
        typeEpsilon.resize(numTypes);
        for (int i = 0; i < numTypes; i++)
            force.getParticleTypeParameters(i, typeSigma[i], typeEpsilon[i]);
    }
    else {
        // Identify types by finding every unique sigma/epsilon pair.

        map<pair<double, double>, int> typeForParams;
        for (int i = 0; i < numParticles; i++) {
            double sigma, epsilon, reduction;
            int parent, typeIndex;
            bool isAlchemical;
            force.getParticleParameters(i, parent, sigma, epsilon, reduction, isAlchemical, typeIndex);
            pair<double, double> params = make_pair(sigma, epsilon);
            map<pair<double, double>, int>::iterator entry = typeForParams.find(params);
            if (entry == typeForParams.end()) {
                int index = typeForParams.size();
                typeForParams[params] = index;
            }
            type[i] = typeForParams[params];
        }
        numTypes = typeForParams.size();
        typeSigma.resize(numTypes);
        typeEpsilon.resize(numTypes);
        for (auto params : typeForParams) {
            typeSigma[params.second] = params.first.first;
            typeEpsilon[params.second] = params.first.second;
        }
    }
    
    // Build the matrices by applying combining rules.

    sigmaMatrix.clear();
    epsilonMatrix.clear();
    sigmaMatrix.resize(numTypes, vector<double>(numTypes));
    epsilonMatrix.resize(numTypes, vector<double>(numTypes));
    string sigmaCombiningRule = force.getSigmaCombiningRule();
    string epsilonCombiningRule = force.getEpsilonCombiningRule();
    for (int i = 0; i < numTypes; i++) {
        double iSigma = typeSigma[i];
        double iEpsilon = typeEpsilon[i];
        for (int j = 0; j < numTypes; j++) {
            double jSigma = typeSigma[j];
            double jEpsilon = typeEpsilon[j];
            double sigma, epsilon;
            // ARITHMETIC = 1
            // GEOMETRIC  = 2
            // CUBIC-MEAN = 3
            if (sigmaCombiningRule == "ARITHMETIC")
              sigma = iSigma+jSigma;
            else if (sigmaCombiningRule == "GEOMETRIC")
              sigma = 2*sqrt(iSigma*jSigma);
            else if (sigmaCombiningRule == "CUBIC-MEAN") {
              double iSigma2 = iSigma*iSigma;
              double jSigma2 = jSigma*jSigma;
              if ((iSigma2+jSigma2) != 0.0)
                sigma = 2*(iSigma2*iSigma + jSigma2*jSigma) / (iSigma2+jSigma2);
              else
                sigma = 0.0;
            }
            else
                throw OpenMMException("AmoebaVdwForce: Unknown value for sigma combining rule: "+sigmaCombiningRule);
            sigmaMatrix[i][j] = sigma;
            sigmaMatrix[j][i] = sigma;

            // ARITHMETIC = 1
            // GEOMETRIC  = 2
            // HARMONIC   = 3
            // W-H        = 4
            // HHG        = 5
            if (epsilonCombiningRule == "ARITHMETIC")
              epsilon = 0.5*(iEpsilon+jEpsilon);
            else if (epsilonCombiningRule == "GEOMETRIC")
              epsilon = sqrt(iEpsilon*jEpsilon);
            else if (epsilonCombiningRule == "HARMONIC") {
              if ((iEpsilon+jEpsilon) != 0.0)
                epsilon = 2*(iEpsilon*jEpsilon) / (iEpsilon+jEpsilon);
              else
                epsilon = 0.0;
            }
            else if (epsilonCombiningRule == "W-H") {
              double iSigma3 = iSigma * iSigma * iSigma;
              double jSigma3 = jSigma * jSigma * jSigma;
              double iSigma6 = iSigma3 * iSigma3;
              double jSigma6 = jSigma3 * jSigma3;
              double eps_s = sqrt(iEpsilon*jEpsilon);
              epsilon = (eps_s == 0.0 ? 0.0 : 2*eps_s*iSigma3*jSigma3/(iSigma6+jSigma6));
            }
            else if (epsilonCombiningRule == "HHG") {
              double epsilonS = sqrt(iEpsilon)+sqrt(jEpsilon);
              if (epsilonS != 0.0)
                epsilon = 4*(iEpsilon*jEpsilon) / (epsilonS*epsilonS);
              else
                epsilon = 0.0;
            }
            else
                throw OpenMMException("AmoebaVdwForce: Unknown value for epsilon combining rule: "+epsilonCombiningRule);
            epsilonMatrix[i][j] = epsilon;
            epsilonMatrix[j][i] = epsilon;
        }
    }
    
    // Record any type pairs that override the combining rules.

    if (force.getUseParticleTypes()) {
        for (int i = 0; i < force.getNumTypePairs(); i++) {
            int type1, type2;
            double sigma, epsilon;
            force.getTypePairParameters(i, type1, type2, sigma, epsilon);
            sigmaMatrix[type1][type2] = sigma;
            sigmaMatrix[type2][type1] = sigma;
            epsilonMatrix[type1][type2] = epsilon;
            epsilonMatrix[type2][type1] = epsilon;
        }
    }
}

double AmoebaVdwForceImpl::calcDispersionCorrection(const System& system, const AmoebaVdwForce& force) {

    // Amoeba VdW dispersion correction implemented by LPW
    // There is no dispersion correction if PBC is off or the cutoff is set to the default value of ten billion (AmoebaVdwForce.cpp)
    if (force.getNonbondedMethod() == AmoebaVdwForce::NoCutoff)
        return 0.0;

    // Identify all particle classes (defined by sigma and epsilon), and count the number of
    // particles in each class.

    vector<int> type;
    vector<vector<double> > sigmaMatrix;
    vector<vector<double> > epsilonMatrix;
    createParameterMatrix(force, type, sigmaMatrix, epsilonMatrix);
    int numTypes = sigmaMatrix.size();
    vector<int> typeCounts(numTypes, 0);
    for (int i = 0; i < force.getNumParticles(); i++)
        typeCounts[type[i]]++;

    // Compute the VdW tapering coefficients.  Mostly copied from amoebaCudaGpu.cpp.
    double cutoff = force.getCutoffDistance();
    double vdwTaper = 0.90; // vdwTaper is a scaling factor, it is not a distance.
    double c0 = 0.0;
    double c1 = 0.0;
    double c2 = 0.0;
    double c3 = 0.0;
    double c4 = 0.0;
    double c5 = 0.0;

    double vdwCut = cutoff;
    double vdwTaperCut = vdwTaper*cutoff;

    double vdwCut2 = vdwCut*vdwCut;
    double vdwCut3 = vdwCut2*vdwCut;
    double vdwCut4 = vdwCut2*vdwCut2;
    double vdwCut5 = vdwCut2*vdwCut3;
    double vdwCut6 = vdwCut3*vdwCut3;
    double vdwCut7 = vdwCut3*vdwCut4;

    double vdwTaperCut2 = vdwTaperCut*vdwTaperCut;
    double vdwTaperCut3 = vdwTaperCut2*vdwTaperCut;
    double vdwTaperCut4 = vdwTaperCut2*vdwTaperCut2;
    double vdwTaperCut5 = vdwTaperCut2*vdwTaperCut3;
    double vdwTaperCut6 = vdwTaperCut3*vdwTaperCut3;
    double vdwTaperCut7 = vdwTaperCut3*vdwTaperCut4;

    // get 5th degree multiplicative switching function coefficients;

    double denom = 1.0 / (vdwCut - vdwTaperCut);
    double denom2 = denom*denom;
    denom = denom * denom2*denom2;

    c0 = vdwCut * vdwCut2 * (vdwCut2 - 5.0 * vdwCut * vdwTaperCut + 10.0 * vdwTaperCut2) * denom;
    c1 = -30.0 * vdwCut2 * vdwTaperCut2*denom;
    c2 = 30.0 * (vdwCut2 * vdwTaperCut + vdwCut * vdwTaperCut2) * denom;
    c3 = -10.0 * (vdwCut2 + 4.0 * vdwCut * vdwTaperCut + vdwTaperCut2) * denom;
    c4 = 15.0 * (vdwCut + vdwTaperCut) * denom;
    c5 = -6.0 * denom;

    // Loop over all pairs of types to compute the coefficient.
    // Copied over from TINKER - numerical integration.
    double range = 20.0;
    double cut = vdwTaperCut; // This is where tapering BEGINS
    double off = vdwCut; // This is where tapering ENDS
    int nstep = 200;
    int ndelta = int(double(nstep) * (range - cut));
    double rdelta = (range - cut) / double(ndelta);
    double offset = cut - 0.5 * rdelta;

    // Buffered-14-7 buffering constants
    double dhal = 0.07; 
    double ghal = 0.12;

    double elrc = 0.0; // This number is incremented and passed out at the end
    double e = 0.0;

    // Double loop over different atom types.
    
    for (int i = 0; i < numTypes; i++) {
        for (int j = 0; j < numTypes; j++) {
            double sigma = sigmaMatrix[i][j];
            double epsilon = epsilonMatrix[i][j];
            int count = typeCounts[i]*typeCounts[j];
            // Below is an exact copy of stuff from the previous block.
            double rv = sigma;
            double termik = 2.0 * M_PI * count; // termik is equivalent to 2 * pi * count.
            double rv2 = rv * rv;
            double rv6 = rv2 * rv2 * rv2;
            double rv7 = rv6 * rv;
            double etot = 0.0;
            double r2 = 0.0;
            for (int j = 1; j <= ndelta; j++) {
                double r = offset + double(j) * rdelta;
                r2 = r*r;
                double r3 = r2 * r;
                double r6 = r3 * r3;
                if (force.getPotentialFunction() == AmoebaVdwForce::LennardJones) {
                    double p6 = rv6 / r6;
                    double p12 = p6 * p6;
                    e = 4 * epsilon * (p12 - p6);
                }
                else {
                    double r7 = r6 * r;
                    double rho = r7 + ghal * rv7;
                    double tau = (dhal+1.0) / (r+dhal*rv);
                    double tau7 = pow(tau, 7);
                    e = epsilon * rv7 * tau7 * ((ghal+1.0)*rv7/rho-2.0);
                }
                double taper = 0.0;
                if (r < off) {
                    double r4 = r2 * r2;
                    double r5 = r2 * r3;
                    taper = c5 * r5 + c4 * r4 + c3 * r3 + c2 * r2 + c1 * r + c0;
                    e = e * (1.0 - taper);
                }
                etot = etot + e * rdelta * r2;
            }
            elrc = elrc + termik * etot;
        }
    }
    return elrc;
}

std::vector<std::string> AmoebaVdwForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcAmoebaVdwForceKernel::Name());
    return names;
}

void AmoebaVdwForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcAmoebaVdwForceKernel>().copyParametersToContext(context, owner);
}



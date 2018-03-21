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

using std::pair;
using std::vector;
using std::set;

AmoebaVdwForceImpl::AmoebaVdwForceImpl(const AmoebaVdwForce& owner)
    : owner(owner) {}

AmoebaVdwForceImpl::~AmoebaVdwForceImpl() {}

void AmoebaVdwForceImpl::initialize(ContextImpl& context) {
    const System& system = context.getSystem();

    if (owner.getNumParticles() != system.getNumParticles())
        throw OpenMMException("AmoebaVdwForce must have exactly as many particles as the System it belongs to.");

    // check that cutoff < 0.5*boxSize

    if (owner.getNonbondedMethod() == AmoebaVdwForce::CutoffPeriodic) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoffDistance();
        if (cutoff > 0.5 * boxVectors[0][0] || cutoff > 0.5 * boxVectors[1][1] || cutoff > 0.5 * boxVectors[2][2])
            throw OpenMMException("AmoebaVdwForce: The cutoff distance cannot be greater than half the periodic box size.");
    }

    kernel = context.getPlatform().createKernel(CalcAmoebaVdwForceKernel::Name(), context);
    kernel.getAs<CalcAmoebaVdwForceKernel>().initialize(context.getSystem(), owner);
}

double AmoebaVdwForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups & (1 << owner.getForceGroup())) != 0)
        return kernel.getAs<CalcAmoebaVdwForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

double AmoebaVdwForceImpl::calcDispersionCorrection(const System& system, const AmoebaVdwForce& force) {
    // Amoeba VdW dispersion correction implemented by LPW
    // There is no dispersion correction if PBC is off or the cutoff is set to the default value of ten billion (AmoebaVdwForce.cpp)
    if (force.getNonbondedMethod() == AmoebaVdwForce::NoCutoff) {
        return 0.0;
    }
    // Identify all particle classes (defined by sigma and epsilon and reduction), and count the number of
    // particles in each class.
    double V1lambda = 1.0;
    map<int, int> NonLigclassCounts;
    map<int, int> LigclassCounts;
    for (int i = 0; i < force.getNumParticles(); i++) {
        double sigma, epsilon, reduction, vlambda;
        // The variables reduction, ivindex are not used.
        int ivindex;
        int key; // key == vdwType
        // Get the sigma and epsilon parameters, ignoring everything else.
        force.getParticleParameters(i, ivindex, key, sigma, epsilon, reduction, vlambda);
        if (vlambda < 1.0) {
            V1lambda = vlambda;
            map<int, int>::iterator entry = LigclassCounts.find(key);
            if (entry == LigclassCounts.end()) {
                LigclassCounts[key] = 1;
                NonLigclassCounts[key] = 0;
            }
            else {
                entry->second++;
            }
        }
        else {
            map<int, int>::iterator entry = NonLigclassCounts.find(key);
            if (entry == NonLigclassCounts.end()) {
                LigclassCounts[key] = 0;
                NonLigclassCounts[key] = 1;
            }
            else {
                entry->second++;
            }
        }
    }

    // Compute the VdW tapering coefficients.
    double cutoff = force.getCutoffDistance();
    double vdwTaper = 0.90; // vdwTaper is a scaling factor, it is not a distance.
    double c0 = 0.0;
    double c1 = 0.0;
    double c2 = 0.0;
    double c3 = 0.0;
    double c4 = 0.0;
    double c5 = 0.0;

    double vdwCut = cutoff;
    double vdwTaperCut = vdwTaper * cutoff;

    double vdwCut2 = vdwCut * vdwCut;
    double vdwCut3 = vdwCut2 * vdwCut;
    double vdwCut4 = vdwCut2 * vdwCut2;
    double vdwCut5 = vdwCut2 * vdwCut3;
    double vdwCut6 = vdwCut3 * vdwCut3;
    double vdwCut7 = vdwCut3 * vdwCut4;

    double vdwTaperCut2 = vdwTaperCut * vdwTaperCut;
    double vdwTaperCut3 = vdwTaperCut2 * vdwTaperCut;
    double vdwTaperCut4 = vdwTaperCut2 * vdwTaperCut2;
    double vdwTaperCut5 = vdwTaperCut2 * vdwTaperCut3;
    double vdwTaperCut6 = vdwTaperCut3 * vdwTaperCut3;
    double vdwTaperCut7 = vdwTaperCut3 * vdwTaperCut4;

    // get 5th degree multiplicative switching function coefficients;

    double denom = 1.0 / (vdwCut - vdwTaperCut);
    double denom2 = denom * denom;
    denom = denom * denom2 * denom2;

    c0 = vdwCut * vdwCut2 * (vdwCut2 - 5.0 * vdwCut * vdwTaperCut + 10.0 * vdwTaperCut2) * denom;
    c1 = -30.0 * vdwCut2 * vdwTaperCut2 * denom;
    c2 = 30.0 * (vdwCut2 * vdwTaperCut + vdwCut * vdwTaperCut2) * denom;
    c3 = -10.0 * (vdwCut2 + 4.0 * vdwCut * vdwTaperCut + vdwTaperCut2) * denom;
    c4 = 15.0 * (vdwCut + vdwTaperCut) * denom;
    c5 = -6.0 * denom;

    // Loop over all pairs of classes to compute the coefficient.
    // Copied over from TINKER - numerical integration.
    double range = 10;
    double cut = vdwTaperCut; // This is where tapering BEGINS
    double off = vdwCut; // This is where tapering ENDS
    int nstep = 20;
    int ndelta = int(double(nstep) * (range - cut));
    double rdelta = (range - cut) / double(ndelta);
    double offset = cut - 0.5 * rdelta;
    double dhal = 0.07; // This magic number also appears in kCalculateAmoebaCudaVdw14_7.cu
    double ghal = 0.12; // This magic number also appears in kCalculateAmoebaCudaVdw14_7.cu
    double elrc = 0.0; // This number is incremented and passed out at the end
    double e = 0.0;
    double sigma, epsilon; // The pairwise sigma and epsilon parameters.
    int i = 0, k = 0; // Loop counters.

    std::string functionalForm = force.getFunctionalForm();

    // Double loop over different atom types.
    for (map<int, int>::iterator class1 = NonLigclassCounts.begin(); class1 != NonLigclassCounts.end(); ++class1) {
        k = 0;
        for (map<int, int>::const_iterator class2 = NonLigclassCounts.begin(); class2 != NonLigclassCounts.end(); ++class2) {
            int key1 = class1->first;
            int key2 = class2->first;
            force.getVdwprParameters(key1, key2, sigma, epsilon);

            int countNon1 = class1->second;
            int countNon2 = class2->second;
            map<int, int>::iterator Ligclass1 = LigclassCounts.find(key1);
            int countLig1 = Ligclass1->second;
            map<int, int>::iterator Ligclass2 = LigclassCounts.find(key2);
            int countLig2 = Ligclass2->second;
            double rv = sigma;
            int total1 = countNon1 + countLig1;
            int total2 = countNon2 + countLig2;

            double termik = 0.0;
            if (force.getCoupleMethod() == AmoebaVdwForce::Decouple) {
                termik = total1 * total2 - (1.0 - V1lambda) * (countLig1 * countNon2 + countNon1 * countLig2); // termik is equivalent to 2*pi*count.
            } else if (force.getCoupleMethod() == AmoebaVdwForce::Annihilate) {
                termik = V1lambda * total1 * total2 + (1.0 - V1lambda) * countNon1 * countNon2;
            }

            termik *= (2.0 * M_PI);
            double rv2 = rv * rv;
            double rv6 = rv2 * rv2 * rv2;
            double rv7 = rv6 * rv;
            double etot = 0.0;
            double r2 = 0.0;
            for (int j = 1; j <= ndelta; j++) {
                double r = offset + double(j) * rdelta;
                r2 = r * r;
                double r3 = r2 * r;
                double r6 = r3 * r3;
                double r7 = r6 * r;
                double rho = r / rv;
                if (functionalForm == "BUFFERED-14-7") {
                    double term1 = pow(((dhal + 1.0) / (dhal + rho)), 7);
                    double term2 = ((ghal + 1.0) / (ghal + pow(rho, 7))) - 2.0;
                    e = epsilon * term1 * term2;
                } else if (functionalForm == "LENNARD-JONES") {
                    double p6 = rv6 / r6;
                    double p12 = p6 * p6;
                    e = epsilon * (p12 - 2.0 * p6);
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
            k++;
        }
        i++;
    }
    std::map<int, int>::iterator it;
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

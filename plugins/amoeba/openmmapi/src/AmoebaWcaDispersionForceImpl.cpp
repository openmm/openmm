/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
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
#include "openmm/internal/AmoebaWcaDispersionForceImpl.h"
#include "openmm/amoebaKernels.h"
#include <cmath>

using namespace OpenMM;

using std::pair;
using std::vector;
using std::set;

AmoebaWcaDispersionForceImpl::AmoebaWcaDispersionForceImpl(const AmoebaWcaDispersionForce& owner) : owner(owner) {
}

AmoebaWcaDispersionForceImpl::~AmoebaWcaDispersionForceImpl() {
}

void AmoebaWcaDispersionForceImpl::initialize(ContextImpl& context) {
 
    const System& system = context.getSystem();
    if (owner.getNumParticles() != system.getNumParticles())
        throw OpenMMException("AmoebaWcaDispersionForce must have exactly as many particles as the System it belongs to.");

    kernel = context.getPlatform().createKernel(CalcAmoebaWcaDispersionForceKernel::Name(), context);
    kernel.getAs<CalcAmoebaWcaDispersionForceKernel>().initialize(context.getSystem(), owner);
}

double AmoebaWcaDispersionForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcAmoebaWcaDispersionForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

static double tailCorrection(double ri, double eps, double rmin) {
    if (ri < rmin) {
        double r3 = ri * ri * ri;
        double rmin3 = rmin * rmin * rmin;
        return -eps * 4.0 * M_PI * (rmin3 - r3) / 3.0 - eps * M_PI * 18.0 * rmin3 / 11.0;
    }
    else {
        double ri2 = ri * ri;
        double ri4 = ri2 * ri2;
        double ri7 = ri * ri2 * ri4;
        double ri11 = ri7 * ri4;
        double rmin2 = rmin * rmin;
        double rmin4 = rmin2 * rmin2;
        double rmin7 = rmin * rmin2 * rmin4;
        return 2.0 * M_PI * eps * rmin7 * (2.0 * rmin7 - 11.0 * ri7) / (11.0 * ri11);
    }
}

void AmoebaWcaDispersionForceImpl::getMaximumDispersionEnergy(const AmoebaWcaDispersionForce &force, int particleIndex,
                                                              double &maxDispersionEnergy) {

    // Atom i parameters
    double rmini, epsi;
    force.getParticleParameters(particleIndex, rmini, epsi);
    if (epsi <= 0.0 || rmini <= 0.0) {
        maxDispersionEnergy = 0.0;
        return;
    }
    double rmini2 = rmini * rmini;
    double sepsi = std::sqrt(epsi);

    // Atom i with water oxygen.
    double epso = force.getEpso();
    double emixo = std::sqrt(epso) + sepsi;
    emixo = 4.0 * epso * epsi / (emixo * emixo);

    double rmino = force.getRmino();
    double rmino2 = rmino * rmino;
    double rmixo = 2.0 * (rmino2 * rmino + rmini2 * rmini) / (rmino2 + rmini2);

    // Atom i with water hydrogen.
    double epsh = force.getEpsh();
    double emixh = std::sqrt(epsh) + sepsi;
    emixh = 4.0 * epsh * epsi / (emixh * emixh);

    double rminh = force.getRminh();
    double rminh2 = rminh * rminh;
    double rmixh = 2.0 * (rminh2 * rminh + rmini2 * rmini) / (rminh2 + rmini2);

    // Compute the tail correction (i.e. the dispersion integral in the absence of any other solute atoms).
    double dispersionOffest = force.getDispoff();
    double riO = rmixo / 2.0 + dispersionOffest;
    double cdisp = tailCorrection(riO, emixo, rmixo);
    double riH = rmixh / 2.0 + dispersionOffest;
    cdisp += 2.0 * tailCorrection(riH, emixh, rmixh);

    maxDispersionEnergy = force.getSlevy() * force.getAwater() * cdisp;
}

double AmoebaWcaDispersionForceImpl::getTotalMaximumDispersionEnergy(const AmoebaWcaDispersionForce &force) {
    double totalMaximumDispersionEnergy = 0.0;
    for (int ii = 0; ii < force.getNumParticles(); ii++) {
        double maximumDispersionEnergy;
        getMaximumDispersionEnergy(force, ii, maximumDispersionEnergy);
        totalMaximumDispersionEnergy += maximumDispersionEnergy;
    }
    return totalMaximumDispersionEnergy;
}

std::vector<std::string> AmoebaWcaDispersionForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcAmoebaWcaDispersionForceKernel::Name());
    return names;
}

void AmoebaWcaDispersionForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcAmoebaWcaDispersionForceKernel>().copyParametersToContext(context, owner);
}

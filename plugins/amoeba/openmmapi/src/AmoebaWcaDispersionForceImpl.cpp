/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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
void AmoebaWcaDispersionForceImpl::getMaximumDispersionEnergy(const AmoebaWcaDispersionForce& force, int particleIndex, double& maxDispersionEnergy) {

    const double pi  = 3.1415926535897;

    // from last loop in subroutine knp in ksolv.f

    double rdisp, epsi;
    force.getParticleParameters(particleIndex, rdisp, epsi);
    if (epsi <= 0.0 || rdisp <= 0.0) {
        maxDispersionEnergy = 0.0;
        return;
    }
    double rmini     = rdisp;
    rdisp           += force.getDispoff();
    double epso      = force.getEpso();
    double emixo     = std::sqrt(epso) + std::sqrt(epsi);
           emixo     = 4.0*epso*epsi/(emixo*emixo);

    double rmino     = force.getRmino();
    double rmino2    = rmino*rmino;
    double rmini2    = rmini*rmini;
    double rmixo     = 2.0*(rmino2*rmino + rmini2*rmini) / (rmino2 + rmini2);

    double rmixo3    = rmixo*rmixo*rmixo;
    double rmixo7    = rmixo*rmixo3*rmixo3;
    double ao        = emixo*rmixo7;

    double epsh      = force.getEpsh();
    double emixh     = std::sqrt(epsh) + std::sqrt(epsi);
           emixh     = 4.0*epsh*epsi/(emixh*emixh);
    double rminh     = force.getRminh();
    double rminh2    = rminh*rminh;
    double rmixh     = rminh*rminh + rmini2;
           rmixh     = 2.0 * (rminh2*rminh + rmini2*rmini) / (rminh2 + rmini2);
    double rmixh3    = rmixh*rmixh*rmixh;
    double rmixh7    = rmixh3*rmixh3*rmixh;
    double ah        = emixh*rmixh7;

    double rdisp3    = rdisp*rdisp*rdisp;
    double rdisp7    = rdisp*rdisp3*rdisp3;
    double rdisp11   = rdisp7*rdisp3*rdisp;

    double cdisp;
    if (rdisp < rmixh) {
        cdisp = -4.0*pi*emixh*(rmixh3-rdisp3)/3.0 - emixh*18.0/11.0*rmixh3*pi;
    } else {
        cdisp = 2.0*pi*(2.0*rmixh7-11.0*rdisp7)*ah/ (11.0*rdisp11);
    }
    cdisp *= 2.0;
    if (rdisp < rmixo) {
        cdisp -= 4.0*pi*emixo*(rmixo3-rdisp3)/3.0;
        cdisp -= emixo*18.0/11.0*rmixo3*pi;
    } else {
        cdisp += 2.0*pi*(2.0*rmixo7-11.0*rdisp7) * ao/(11.0*rdisp11);
    }
    maxDispersionEnergy = force.getSlevy()*force.getAwater()*cdisp;
}

double AmoebaWcaDispersionForceImpl::getTotalMaximumDispersionEnergy(const AmoebaWcaDispersionForce& force) {

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

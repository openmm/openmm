/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "openmm/internal/QTBIntegratorUtilities.h"
#include "SimTKOpenMMRealType.h"
#include <map>

using namespace OpenMM;
using namespace std;

void QTBIntegratorUtilities::findTypes(const System& system, const QTBIntegrator& integrator, vector<int>& particleType,
        vector<std::vector<int> >& typeParticles, vector<double>& typeMass, vector<double>& typeAdaptationRate) {
    // Record information about groups defined by particle types.

    particleType.resize(system.getNumParticles());
    map<int, int> typeIndex;
    map<int, double> massTable;
    const auto& types = integrator.getParticleTypes();
    double defaultAdaptationRate = integrator.getDefaultAdaptationRate();
    for (auto particle : types) {
        int type = particle.second;
        double mass = system.getParticleMass(particle.first);
        if (typeIndex.find(type) == typeIndex.end()) {
            int index = typeIndex.size();
            typeIndex[type] = index;
            double rate = defaultAdaptationRate;
            const auto& typeRates = integrator.getTypeAdaptationRates();
            if (typeRates.find(type) != typeRates.end())
                rate = typeRates.at(type);
            typeAdaptationRate.push_back(rate);
            typeParticles.push_back(vector<int>());
            typeMass.push_back(mass);
            massTable[type] = mass;
        }
        if (mass != massTable[type])
            throw OpenMMException("QTBIntegrator: All particles of the same type must have the same mass");
        particleType[particle.first] = typeIndex[type];
        typeParticles[typeIndex[type]].push_back(particle.first);
    }
    for (int i = 0; i < system.getNumParticles(); i++)
        if (types.find(i) == types.end()) {
            // This particle's type isn't set, so define a new type for it.
            particleType[i] = typeParticles.size();
            typeAdaptationRate.push_back(defaultAdaptationRate);
            typeParticles.push_back({i});
            typeMass.push_back(system.getParticleMass(i));
        }
}

void QTBIntegratorUtilities::calculateSpectrum(double temperature, double friction, double dt, int numFreq, vector<double>& theta, vector<double>& thetad, ThreadPool& threads) {
    // Compute the standard spectrum.

    double hbar = 1.054571628e-34*AVOGADRO/(1000*1e-12);
    double kT = BOLTZ*temperature;
    theta.resize(numFreq);
    theta[0] = kT;
    for (int i = 1; i < numFreq; i++) {
        double w = M_PI*i/(numFreq*dt);
        theta[i] = hbar*w*(0.5+1/(exp(hbar*w/kT)-1));
    }

    // Compute the deconvolved version.  The algorithm is described in the supplementary
    // information in https://doi.org/10.1021/acs.jpclett.1c01722.

    auto C = [&](double w0, double w) {
        double t = w*w-w0*w0;
        return (friction/M_PI)*w0*w0/(t*t+friction*friction*w*w);
    };
    double dw = M_PI/(numFreq*dt);

    // Normalize the kernel to reduce error at low frequencies.

    vector<double> sum(numFreq, 0.0), scale(numFreq);
    for (int i = 0; i < numFreq; i++) {
        double wi = M_PI*(i+0.5)/(numFreq*dt);
        for (int j = 0; j < numFreq; j++) {
            double wj = M_PI*(j+0.5)/(numFreq*dt);
            sum[i] += C(wi, wj)*dw;
        }
    }
    for (int i = 0; i < numFreq; i++)
        scale[i] = 0.5/sum[i];

    // Compute intermediate quantities.

    vector<vector<double> > D(numFreq, vector<double>(numFreq));
    vector<double> h(numFreq);
    vector<double> fcurrent(numFreq), fnext(numFreq);
    for (int i = 0; i < numFreq; i++)
        fcurrent[i] = 0.5*theta[i];
    threads.execute([&] (ThreadPool& threads, int threadIndex) {
        for (int i = threadIndex; i < numFreq; i += threads.getNumThreads()) {
            double wi = M_PI*(i+0.5)/(numFreq*dt);
            h[i] = 0.0;
            for (int j = 0; j < numFreq; j++) {
                double wj = M_PI*(j+0.5)/(numFreq*dt);
                h[i] += dw*C(wj, wi)*fcurrent[j]*scale[j];
                D[i][j] = 0.0;
                for (int k = 0; k < numFreq; k++) {
                    double wk = M_PI*(k+0.5)/(numFreq*dt);
                    D[i][j] += dw*C(wk, wi)*C(wk, wj)*scale[k]*scale[k];
                }
            }
        }
    });
    threads.waitForThreads();

    // Perform the iteration.

    for (int iteration = 0; iteration < 20; iteration++) {
        for (int i = 0; i < numFreq; i++) {
            double denom = 0.0;
            for (int j = 0; j < numFreq; j++)
                denom += dw*D[i][j]*fcurrent[j];
            fnext[i] = fcurrent[i]*h[i]/denom;
        }
        fcurrent = fnext;
    }
    thetad = fnext;
}

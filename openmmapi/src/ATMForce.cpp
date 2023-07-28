/* -------------------------------------------------------------------------- *
 *                    OpenMM's Alchemical Transfer Force                      *
 * -------------------------------------------------------------------------- *
 * This is a Force of the OpenMM molecular simulation toolkit                 *
 * that implements the Alchemical Transfer Potential                          *
 * for absolute and relative binding free energy estimation                   *
 * (https://doi.org/10.1021/acs.jcim.1c01129). The code is derived from the   *
 * ATMMetaForce plugin                                                        *
 * https://github.com/Gallicchio-Lab/openmm-atmmetaforce-plugin               *
 * with support from the National Science Foundation CAREER 1750511           *
 *                                                                            *
 * Portions copyright (c) 2021-2023 by the Authors                            *
 * Authors: Emilio Gallicchio                                                 *
 * Contributors: Peter Eastman                                                *
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

#include "openmm/ATMForce.h"
#include "openmm/Force.h"
#include "openmm/serialization/XmlSerializer.h"
#include "openmm/internal/ATMForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

ATMForce::ATMForce(double lambda1, double lambda2, double alpha, double u0, double w0, double umax, double ubcore, double acore, double direction):
    defaultLambda1(lambda1), defaultLambda2(lambda2), defaultAlpha(alpha), defaultU0(u0), defaultW0(w0),
    defaultUmax(umax), defaultUbcore(ubcore), defaultAcore(acore), defaultDirection(direction) {
}

ATMForce::~ATMForce() {
    for (Force* force : forces)
        delete force;
}

int ATMForce::addParticle(const std::vector<double>& displacement) {
    double dx, dy, dz, dx0, dy0, dz0;
    if ( displacement.size() == 3 ) {
      dx = displacement[0];
      dy = displacement[1];
      dz = displacement[2];
      dx0 = 0;
      dy0 = 0;
      dz0 = 0;
    } else if ( displacement.size() == 6 ){
      dx = displacement[0];
      dy = displacement[1];
      dz = displacement[2];
      dx0 = displacement[3];
      dy0 = displacement[4];
      dz0 = displacement[5];
    } else {
      throw OpenMMException("ATMForce::addParticle(): the displacement vector must have either 3 or 6 elements");
    }
    particles.push_back(ParticleInfo(particles.size(), dx, dy, dz, dx0, dy0, dz0 ));
    return particles.size()-1;
}

void ATMForce::getParticleParameters(int index, std::vector<double>& displacement) const {
    ASSERT_VALID_INDEX(index, particles);
    displacement.resize(6);
    displacement[0] = particles[index].dx;
    displacement[1] = particles[index].dy;
    displacement[2] = particles[index].dz;
    displacement[3] = particles[index].dx0;
    displacement[4] = particles[index].dy0;
    displacement[5] = particles[index].dz0;
}

void ATMForce::setParticleParameters(int index, const std::vector<double>& displacement) {
    ASSERT_VALID_INDEX(index, particles);
    double dx, dy, dz, dx0, dy0, dz0;
    if ( displacement.size() == 3 ) {
      dx = displacement[0];
      dy = displacement[1];
      dz = displacement[2];
      dx0 = 0;
      dy0 = 0;
      dz0 = 0;
    } else if ( displacement.size() == 6 ){
      dx = displacement[0];
      dy = displacement[1];
      dz = displacement[2];
      dx0 = displacement[3];
      dy0 = displacement[4];
      dz0 = displacement[5];
    } else {
      throw OpenMMException("ATMForce::setParticleParameters(): the displacement vector must have either 3 or 6 elements");
    }
    particles[index].dx = dx;
    particles[index].dy = dy;
    particles[index].dz = dz;
    particles[index].dx0 = dx0;
    particles[index].dy0 = dy0;
    particles[index].dz0 = dz0;
}

int ATMForce::addForce(Force* force) {
    forces.push_back(force);
    return forces.size()-1;
}

Force& ATMForce::getForce(int index) const {
    ASSERT_VALID_INDEX(index, forces);
    return *forces[index];
}

ForceImpl* ATMForce::createImpl() const {
    return new ATMForceImpl(*this);
}

void ATMForce::updateParametersInContext(OpenMM::Context& context) {
    dynamic_cast<ATMForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

double ATMForce::getPerturbationEnergy(const OpenMM::Context& context) const {
    return dynamic_cast<const ATMForceImpl&>(getImplInContext(context)).getPerturbationEnergy();
}


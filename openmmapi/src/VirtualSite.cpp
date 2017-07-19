/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012-2017 Stanford University and the Authors.      *
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

#include "openmm/VirtualSite.h"
#include "openmm/OpenMMException.h"
#include <cmath>
#include <vector>

using namespace OpenMM;
using namespace std;

void VirtualSite::setParticles(const vector<int>& particleIndices) {
    particles = particleIndices;
}

int VirtualSite::getNumParticles() const {
    return particles.size();
}

int VirtualSite::getParticle(int particle) const {
    return particles[particle];
}

TwoParticleAverageSite::TwoParticleAverageSite(int particle1, int particle2, double weight1, double weight2) :
        weight1(weight1), weight2(weight2) {
    vector<int> particles(2);
    particles[0] = particle1;
    particles[1] = particle2;
    setParticles(particles);
}

double TwoParticleAverageSite::getWeight(int particle) const {
    if (particle == 0)
        return weight1;
    if (particle == 1)
        return weight2;
    throw OpenMMException("Illegal index for particle");
}

ThreeParticleAverageSite::ThreeParticleAverageSite(int particle1, int particle2, int particle3, double weight1, double weight2, double weight3) :
        weight1(weight1), weight2(weight2), weight3(weight3) {
    vector<int> particles(3);
    particles[0] = particle1;
    particles[1] = particle2;
    particles[2] = particle3;
    setParticles(particles);
}

double ThreeParticleAverageSite::getWeight(int particle) const {
    if (particle == 0)
        return weight1;
    if (particle == 1)
        return weight2;
    if (particle == 2)
        return weight3;
    throw OpenMMException("Illegal index for particle");
}

OutOfPlaneSite::OutOfPlaneSite(int particle1, int particle2, int particle3, double weight12, double weight13, double weightCross) :
        weight12(weight12), weight13(weight13), weightCross(weightCross) {
    vector<int> particles(3);
    particles[0] = particle1;
    particles[1] = particle2;
    particles[2] = particle3;
    setParticles(particles);
}

double OutOfPlaneSite::getWeight12() const {
    return weight12;
}

double OutOfPlaneSite::getWeight13() const {
    return weight13;
}

double OutOfPlaneSite::getWeightCross() const {
    return weightCross;
}

LocalCoordinatesSite::LocalCoordinatesSite(const vector<int>& particles, const vector<double>& originWeights, const vector<double>& xWeights, const vector<double>& yWeights, const Vec3& localPosition) :
        originWeights(originWeights), xWeights(xWeights), yWeights(yWeights), localPosition(localPosition) {
    int numParticles = particles.size();
    if (numParticles < 2)
        throw OpenMMException("LocalCoordinatesSite: Must depend on at least two other particles");
    if (originWeights.size() != numParticles || xWeights.size() != numParticles || yWeights.size() != numParticles)
        throw OpenMMException("LocalCoordinatesSite: Number of weights does not match number of particles");
    double originSum = 0, xSum = 0, ySum = 0;
    for (int i = 0; i < numParticles; i++) {
        originSum += originWeights[i];
        xSum += xWeights[i];
        ySum += yWeights[i];
    }
    if (fabs(originSum-1.0) > 1e-6)
        throw OpenMMException("LocalCoordinatesSite: Weights for computing origin must add to 1");
    if (fabs(xSum) > 1e-6)
        throw OpenMMException("LocalCoordinatesSite: Weights for computing x axis must add to 0");
    if (fabs(ySum) > 1e-6)
        throw OpenMMException("LocalCoordinatesSite: Weights for computing y axis must add to 0");
    setParticles(particles);
}

LocalCoordinatesSite::LocalCoordinatesSite(int particle1, int particle2, int particle3, const Vec3& originWeights, const Vec3& xWeights, const Vec3& yWeights, const Vec3& localPosition) :
        localPosition(localPosition) {
    if (fabs(originWeights[0]+originWeights[1]+originWeights[2]-1.0) > 1e-6)
        throw OpenMMException("LocalCoordinatesSite: Weights for computing origin must add to 1");
    if (fabs(xWeights[0]+xWeights[1]+xWeights[2]) > 1e-6)
        throw OpenMMException("LocalCoordinatesSite: Weights for computing x axis must add to 0");
    if (fabs(yWeights[0]+yWeights[1]+yWeights[2]) > 1e-6)
        throw OpenMMException("LocalCoordinatesSite: Weights for computing y axis must add to 0");
    vector<int> particles(3);
    particles[0] = particle1;
    particles[1] = particle2;
    particles[2] = particle3;
    setParticles(particles);
    this->originWeights.push_back(originWeights[0]);
    this->originWeights.push_back(originWeights[1]);
    this->originWeights.push_back(originWeights[2]);
    this->xWeights.push_back(xWeights[0]);
    this->xWeights.push_back(xWeights[1]);
    this->xWeights.push_back(xWeights[2]);
    this->yWeights.push_back(yWeights[0]);
    this->yWeights.push_back(yWeights[1]);
    this->yWeights.push_back(yWeights[2]);
}

void LocalCoordinatesSite::getOriginWeights(vector<double>& weights) const {
    weights = originWeights;
}

Vec3 LocalCoordinatesSite::getOriginWeights() const {
    if (originWeights.size() != 3)
        throw OpenMMException("LocalCoordinatesSite: This version of getOriginWeights() requires the site to depend on three particles");
    return Vec3(originWeights[0], originWeights[1], originWeights[2]);
}

void LocalCoordinatesSite::getXWeights(vector<double>& weights) const {
    weights = xWeights;
}

Vec3 LocalCoordinatesSite::getXWeights() const {
    if (xWeights.size() != 3)
        throw OpenMMException("LocalCoordinatesSite: This version of getXWeights() requires the site to depend on three particles");
    return Vec3(xWeights[0], xWeights[1], xWeights[2]);
}

void LocalCoordinatesSite::getYWeights(vector<double>& weights) const {
    weights = yWeights;
}

Vec3 LocalCoordinatesSite::getYWeights() const {
    if (yWeights.size() != 3)
        throw OpenMMException("LocalCoordinatesSite: This version of getYWeights() requires the site to depend on three particles");
    return Vec3(yWeights[0], yWeights[1], yWeights[2]);
}

const Vec3& LocalCoordinatesSite::getLocalPosition() const {
    return localPosition;
}

#ifndef OPENMM_VIRTUALSITE_H_
#define OPENMM_VIRTUALSITE_H_

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

#include <openmm/Vec3.h>
#include "internal/windowsExport.h"
#include <vector>

namespace OpenMM {

/**
 * A VirtualSite describes the rules for computing a particle's position based on
 * other particles.  This is an abstract class.  Subclasses define particular rules.
 * To define a virtual site, create an instance of a VirtualSite subclass and then
 * call setVirtualSite() on the System.
 */

class OPENMM_EXPORT VirtualSite {
public:
    virtual ~VirtualSite() {
    }
    /**
     * Get the number of particles this virtual site depends on.
     */
    int getNumParticles() const;
    /**
     * Get the index of a particle this virtual site depends on.
     * 
     * @param particle    the particle to get (between 0 and getNumParticles())
     * @return the index of the particle in the System
     */
    int getParticle(int particle) const;
protected:
    VirtualSite() {
    }
    void setParticles(const std::vector<int>& particleIndices);
private:
    std::vector<int> particles;
};

/**
 * This is a VirtualSite that computes the particle location as a weighted average
 * of two other particle's locations.  Assuming the weights add up to 1, this means
 * the virtual site is on the line passing through the two particles.
 */
class OPENMM_EXPORT TwoParticleAverageSite : public VirtualSite {
public:
    /**
     * Create a new TwoParticleAverageSite virtual site.  Normally weight1 and weight2
     * should add up to 1, although this is not strictly required.
     * 
     * @param particle1    the index of the first particle
     * @param particle2    the index of the second particle
     * @param weight1      the weight factor (between 0 and 1) for the first particle
     * @param weight2      the weight factor (between 0 and 1) for the second particle
     */
    TwoParticleAverageSite(int particle1, int particle2, double weight1, double weight2);
    /**
     * Get the weight factor used for a particle this virtual site depends on.
     * 
     * @param particle    the particle to get (between 0 and getNumParticles())
     * @return the weight factor used for that particle
     */
    double getWeight(int particle) const;
private:
    double weight1, weight2;
};

/**
 * This is a VirtualSite that computes the particle location as a weighted average
 * of three other particle's locations.  Assuming the weights add up to 1, this means
 * the virtual site is in the plane of the three particles.
 */
class OPENMM_EXPORT ThreeParticleAverageSite : public VirtualSite {
public:
    /**
     * Create a new ThreeParticleAverageSite virtual site.  Normally the weights
     * should add up to 1, although this is not strictly required.
     * 
     * @param particle1    the index of the first particle
     * @param particle2    the index of the second particle
     * @param particle3    the index of the third particle
     * @param weight1      the weight factor (between 0 and 1) for the first particle
     * @param weight2      the weight factor (between 0 and 1) for the second particle
     * @param weight2      the weight factor (between 0 and 1) for the third particle
     */
    ThreeParticleAverageSite(int particle1, int particle2, int particle3, double weight1, double weight2, double weight3);
    /**
     * Get the weight factor used for a particle this virtual site depends on.
     * 
     * @param particle    the particle to get (between 0 and getNumParticles())
     * @return the weight factor used for that particle
     */
    double getWeight(int particle) const;
private:
    double weight1, weight2, weight3;
};

/**
 * This is a VirtualSite that computes the particle location based on three other
 * particles' locations.  If r<sub>1</sub> is the location of particle 1,
 * r<sub>12</sub> is the vector from particle 1 to particle 2, and
 * r<sub>13</sub> is the vector from particle 1 to particle 3, then the virtual
 * site location is given by
 * 
 * r<sub>1</sub> + w<sub>12</sub>r<sub>12</sub> + w<sub>13</sub>r<sub>13</sub> + w<sub>cross</sub>(r<sub>12</sub> x r<sub>13</sub>)
 * 
 * The three weight factors are user-specified.  This allows the virtual site location
 * to be out of the plane of the three particles.
 */
class OPENMM_EXPORT OutOfPlaneSite : public VirtualSite {
public:
    /**
     * Create a new OutOfPlaneSite virtual site.
     * 
     * @param particle1    the index of the first particle
     * @param particle2    the index of the second particle
     * @param particle3    the index of the third particle
     * @param weight12     the weight factor for the vector from particle1 to particle2
     * @param weight13     the weight factor for the vector from particle1 to particle3
     * @param weightCross  the weight factor for the cross product
     */
    OutOfPlaneSite(int particle1, int particle2, int particle3, double weight12, double weight13, double weightCross);
    /**
     * Get the weight factor for the vector from particle1 to particle2.
     */
    double getWeight12() const;
    /**
     * Get the weight factor for the vector from particle1 to particle3.
     */
    double getWeight13() const;
    /**
     * Get the weight factor for the cross product.
     */
    double getWeightCross() const;
private:
    double weight12, weight13, weightCross;
};

/**
 * This is a VirtualSite that uses the locations of several other particles to compute a local
 * coordinate system, then places the virtual site at a fixed location in that coordinate
 * system.  The origin of the coordinate system and the directions of its x and y axes
 * are each specified as a weighted sum of the locations of the other particles:
 * 
 * origin = w<sup>o</sup><sub>1</sub>r<sub>1</sub> +  w<sup>o</sup><sub>2</sub>r<sub>2</sub> +  w<sup>o</sup><sub>3</sub>r<sub>3</sub> + ...
 * 
 * xdir = w<sup>x</sup><sub>1</sub>r<sub>1</sub> +  w<sup>x</sup><sub>2</sub>r<sub>2</sub> +  w<sup>x</sup><sub>3</sub>r<sub>3</sub> + ...
 * 
 * ydir = w<sup>y</sup><sub>1</sub>r<sub>1</sub> +  w<sup>y</sup><sub>2</sub>r<sub>2</sub> +  w<sup>y</sup><sub>3</sub>r<sub>3</sub> + ...
 * 
 * For the origin, the weights must add to one.  For example if
 * (w<sup>o</sup><sub>1</sub>, w<sup>o</sup><sub>2</sub>, w<sup>o</sup><sub>3</sub>) = (1.0, 0.0, 0.0),
 * the origin of the local coordinate system is at the location of particle 1.  For xdir and ydir,
 * the weights must add to zero.  For example, if
 * (w<sup>x</sup><sub>1</sub>, w<sup>x</sup><sub>2</sub>, w<sup>x</sup><sub>3</sub>) = (-1.0, 0.5, 0.5),
 * the x axis points from particle 1 toward the midpoint between particles 2 and 3.
 * 
 * The z direction is computed as zdir = xdir x ydir.  To ensure the axes are all orthogonal,
 * ydir is then recomputed as ydir = zdir x xdir.  All three axis vectors are then normalized, and
 * the virtual site location is set to
 * 
 * origin + x*xdir + y*ydir + z*zdir
 */
class OPENMM_EXPORT LocalCoordinatesSite : public VirtualSite {
public:
    /**
     * Create a new LocalCoordinatesSite virtual site.
     * 
     * @param particles      the indices of the particle the site depends on
     * @param originWeights  the weight factors for the particles when computing the origin location
     * @param xWeights       the weight factors for the particles when computing xdir
     * @param yWeights       the weight factors for the particles when computing ydir
     * @param localPosition  the position of the virtual site in the local coordinate system
     */
    LocalCoordinatesSite(const std::vector<int>& particles, const std::vector<double>& originWeights, const std::vector<double>& xWeights, const std::vector<double>& yWeights, const Vec3& localPosition);
    /**
     * Create a new LocalCoordinatesSite virtual site.  This constructor assumes the site depends on
     * exactly three other particles.
     * 
     * @param particle1      the index of the first particle
     * @param particle2      the index of the second particle
     * @param particle3      the index of the third particle
     * @param originWeights  the weight factors for the three particles when computing the origin location
     * @param xWeights       the weight factors for the three particles when computing xdir
     * @param yWeights       the weight factors for the three particles when computing ydir
     * @param localPosition  the position of the virtual site in the local coordinate system
     */
    LocalCoordinatesSite(int particle1, int particle2, int particle3, const Vec3& originWeights, const Vec3& xWeights, const Vec3& yWeights, const Vec3& localPosition);
    /**
     * Get the weight factors for the particles when computing the origin location.
     */
    void getOriginWeights(std::vector<double>& weights) const;
    /**
     * Get the weight factors for the particles when computing the origin location.
     * This version assumes the site depends on exactly three other particles, and
     * throws an exception if that is not the case.
     */
    Vec3 getOriginWeights() const;
    /**
     * Get the weight factors for the particles when computing xdir.
     */
    void getXWeights(std::vector<double>& weights) const;
    /**
     * Get the weight factors for the particles when computing xdir.
     * This version assumes the site depends on exactly three other particles, and
     * throws an exception if that is not the case.
     */
    Vec3 getXWeights() const;
    /**
     * Get the weight factors for the particles when computing ydir.
     */
    void getYWeights(std::vector<double>& weights) const;
    /**
     * Get the weight factors for the particles when computing ydir.
     * This version assumes the site depends on exactly three other particles, and
     * throws an exception if that is not the case.
     */
    Vec3 getYWeights() const;
    /**
     * Get the position of the virtual site in the local coordinate system.
     */
    const Vec3& getLocalPosition() const;
private:
    std::vector<double> originWeights, xWeights, yWeights;
    Vec3 localPosition;
};

} // namespace OpenMM

#endif /*OPENMM_VIRTUALSITE_H_*/

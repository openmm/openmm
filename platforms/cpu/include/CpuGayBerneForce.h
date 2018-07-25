/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2016-2018 Stanford University and the Authors.      *
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

#ifndef OPENMM_CPU_GAYBERNEFORCE_H__
#define OPENMM_CPU_GAYBERNEFORCE_H__

#include "openmm/GayBerneForce.h"
#include "openmm/internal/ThreadPool.h"
#include "CpuNeighborList.h"
#include "CpuPlatform.h"
#include "openmm/Vec3.h"
#include <set>
#include <utility>

namespace OpenMM {

class CpuGayBerneForce {
public:
    struct Matrix;

    /**
     * Constructor.
     */
    CpuGayBerneForce(const GayBerneForce& force);

    /**
     * Compute the interaction.
     *
     * @param positions     the positions of the atoms
     * @param forces        forces will be added to this vector
     * @param threadForce   individual threads can add their forces to this vector
     * @param boxVectors    the periodic box vectors
     * @param data          the platform data for the current context
     * @return the energy of the interaction
     */
    double calculateForce(const std::vector<Vec3>& positions, std::vector<Vec3>& forces, std::vector<AlignedArray<float> >& threadForce, Vec3* boxVectors, CpuPlatform::PlatformData& data);

    /**
     * This routine contains the code executed by each thread.
     */
    void threadComputeForce(ThreadPool& threads, int threadIndex, CpuNeighborList* neighborList);
    
    /**
     * Get the exclusions being used by the force.
     */
    const std::vector<std::set<int> >& getExclusions() const;

private:
    struct ParticleInfo;
    struct ExceptionInfo;
    std::vector<ParticleInfo> particles;
    std::vector<ExceptionInfo> exceptions;
    std::set<std::pair<int, int> > exclusions;
    std::vector<std::set<int> > particleExclusions;
    GayBerneForce::NonbondedMethod nonbondedMethod;
    double cutoffDistance, switchingDistance;
    bool useSwitchingFunction;
    std::vector<double> s;
    std::vector<Matrix> A, B, G;
    std::vector<double> threadEnergy;
    std::vector<std::vector<Vec3> > threadTorque;
    // The following variables are used to make information accessible to the individual threads.
    Vec3 const* positions;
    std::vector<AlignedArray<float> >* threadForce;
    Vec3* boxVectors;
    std::atomic<int> atomicCounter;

    void computeEllipsoidFrames(const std::vector<Vec3>& positions);
    
    void applyTorques(const std::vector<Vec3>& positions, std::vector<Vec3>& forces);

    double computeOneInteraction(int particle1, int particle2, double sigma, double epsilon, const Vec3* positions,
            float* forces, std::vector<Vec3>& torques, const Vec3* boxVectors);
};

struct CpuGayBerneForce::ParticleInfo {
    int xparticle, yparticle;
    double sigmaOver2, sqrtEpsilon, rx, ry, rz, ex, ey, ez;
    bool isPointParticle;
};

struct CpuGayBerneForce::ExceptionInfo {
    int particle1, particle2;
    double sigma, epsilon;
};

struct CpuGayBerneForce::Matrix {
    double v[3][3];
    Vec3 operator*(const Vec3& r) {
        return Vec3(v[0][0]*r[0] + v[0][1]*r[1] + v[0][2]*r[2],
                    v[1][0]*r[0] + v[1][1]*r[1] + v[1][2]*r[2],
                    v[2][0]*r[0] + v[2][1]*r[1] + v[2][2]*r[2]);
    }

    Matrix operator+(const Matrix& m) {
        Matrix result;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                result.v[i][j] = v[i][j]+m.v[i][j];
        return result;
    }

    double determinant() {
        return (v[0][0]*v[1][1]*v[2][2] + v[0][1]*v[1][2]*v[2][0] + v[0][2]*v[1][0]*v[2][1] -
                v[0][0]*v[1][2]*v[2][1] - v[0][1]*v[1][0]*v[2][2] - v[0][2]*v[1][1]*v[2][0]);
    }
    
    Matrix inverse() {
        double invDet = 1/determinant();
        Matrix result;
        result.v[0][0] = invDet*(v[1][1]*v[2][2] - v[1][2]*v[2][1]);
        result.v[1][0] = -invDet*(v[1][0]*v[2][2] - v[1][2]*v[2][0]);
        result.v[2][0] = invDet*(v[1][0]*v[2][1] - v[1][1]*v[2][0]);
        result.v[0][1] = -invDet*(v[0][1]*v[2][2] - v[0][2]*v[2][1]);
        result.v[1][1] = invDet*(v[0][0]*v[2][2] - v[0][2]*v[2][0]);
        result.v[2][1] = -invDet*(v[0][0]*v[2][1] - v[0][1]*v[2][0]);
        result.v[0][2] = invDet*(v[0][1]*v[1][2] - v[0][2]*v[1][1]);
        result.v[1][2] = -invDet*(v[0][0]*v[1][2] - v[0][2]*v[1][0]);
        result.v[2][2] = invDet*(v[0][0]*v[1][1] - v[0][1]*v[1][0]);
        return result;
    }
};

static Vec3 operator*(const Vec3& r, CpuGayBerneForce::Matrix& m) {
    return Vec3(m.v[0][0]*r[0] + m.v[1][0]*r[1] + m.v[2][0]*r[2],
                m.v[0][1]*r[0] + m.v[1][1]*r[1] + m.v[2][1]*r[2],
                m.v[0][2]*r[0] + m.v[1][2]*r[1] + m.v[2][2]*r[2]);
}

} // namespace OpenMM

#endif // OPENMM_CPU_GAYBERNEFORCE_H__

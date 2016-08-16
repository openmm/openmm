/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2016 Stanford University and the Authors.           *
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

#ifndef __ReferenceGayBerneForce_H__
#define __ReferenceGayBerneForce_H__

#include "openmm/GayBerneForce.h"
#include "RealVec.h"
#include <set>
#include <utility>

namespace OpenMM {

class OPENMM_EXPORT ReferenceGayBerneForce {
public:
    struct Matrix;
    /**
     * Constructor.
     */
    ReferenceGayBerneForce(const GayBerneForce& force);

    /**
     * Compute the interaction.
     *
     * @param positions     the positions of the atoms
     * @param forces        forces will be added to this vector
     * @param boxVectors    the periodic box vectors
     * @return the energy of the interaction
     */
    RealOpenMM calculateForce(const std::vector<RealVec>& positions, std::vector<RealVec>& forces, const RealVec* boxVectors);

private:
    struct ParticleInfo;
    struct ExceptionInfo;
    std::vector<ParticleInfo> particles;
    std::vector<ExceptionInfo> exceptions;
    std::set<std::pair<int, int> > exclusions;
    GayBerneForce::NonbondedMethod nonbondedMethod;
    RealOpenMM cutoffDistance, switchingDistance;
    bool useSwitchingFunction;
    std::vector<RealOpenMM> s;
    std::vector<Matrix> A, B, G;

    void computeEllipsoidFrames(const std::vector<RealVec>& positions);
    
    void applyTorques(const std::vector<RealVec>& positions, std::vector<RealVec>& forces, const std::vector<RealVec>& torques);

    RealOpenMM computeOneInteraction(int particle1, int particle2, RealOpenMM sigma, RealOpenMM epsilon, const std::vector<RealVec>& positions,
            std::vector<RealVec>& forces, std::vector<RealVec>& torques, const RealVec* boxVectors);
};

struct ReferenceGayBerneForce::ParticleInfo {
    int xparticle, yparticle;
    RealOpenMM sigma, epsilon, rx, ry, rz, ex, ey, ez;
};

struct ReferenceGayBerneForce::ExceptionInfo {
    int particle1, particle2;
    RealOpenMM sigma, epsilon;
};

struct ReferenceGayBerneForce::Matrix {
    RealOpenMM v[3][3];
    RealVec operator*(const RealVec& r) {
        return RealVec(v[0][0]*r[0] + v[0][1]*r[1] + v[0][2]*r[2],
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

    RealOpenMM determinant() {
        return (v[0][0]*v[1][1]*v[2][2] + v[0][1]*v[1][2]*v[2][0] + v[0][2]*v[1][0]*v[2][1] -
                v[0][0]*v[1][2]*v[2][1] - v[0][1]*v[1][0]*v[2][2] - v[0][2]*v[1][1]*v[2][0]);
    }
    
    Matrix inverse() {
        RealOpenMM invDet = 1/determinant();
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

static RealVec operator*(const RealVec& r, ReferenceGayBerneForce::Matrix& m) {
    return RealVec(m.v[0][0]*r[0] + m.v[1][0]*r[1] + m.v[2][0]*r[2],
                   m.v[0][1]*r[0] + m.v[1][1]*r[1] + m.v[2][1]*r[2],
                   m.v[0][2]*r[0] + m.v[1][2]*r[1] + m.v[2][2]*r[2]);
}

} // namespace OpenMM

#endif // __ReferenceGayBerneForce_H__

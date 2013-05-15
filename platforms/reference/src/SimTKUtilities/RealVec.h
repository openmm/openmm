#ifndef OPENMM_REALVEC_H_
#define OPENMM_REALVEC_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2013 Stanford University and the Authors.      *
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

#include "SimTKOpenMMRealType.h"
#include "openmm/Vec3.h"
#include <cassert>
#include <iosfwd>

namespace OpenMM {

/**
 * This is identical to Vec3, except that the components are of type RealOpenMM, so
 * it can be compiled in either single or double precision.  Automatic conversion
 * between this class and Vec3 is supported.
 */

class RealVec {
public:
    /**
     * Create a RealVec whose elements are all 0.
     */
    RealVec() {
        data[0] = data[1] = data[2] = 0.0;
    }
    /**
     * Create a RealVec with specified x, y, and z components.
     */
    RealVec(RealOpenMM x, RealOpenMM y, RealOpenMM z) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    /**
     * Create a RealVec from a Vec3.
     */
    RealVec(Vec3 v) {
        data[0] = v[0];
        data[1] = v[1];
        data[2] = v[2];
    }
    /**
     * Create a Vec3 from a RealVec.
     */
    operator Vec3() const {
        return Vec3(data[0], data[1], data[2]);
    }
    RealOpenMM operator[](int index) const {
        assert(index >= 0 && index < 3);
        return data[index];
    }
    RealOpenMM& operator[](int index) {
        assert(index >= 0 && index < 3);
        return data[index];
    }

    // Arithmetic operators

    // unary plus
    RealVec operator+() const {
        return RealVec(*this);
    }

    // plus
    RealVec operator+(const RealVec& rhs) const {
        const RealVec& lhs = *this;
        return RealVec(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]);
    }

    RealVec& operator+=(const RealVec& rhs) {
        data[0] += rhs[0];
        data[1] += rhs[1];
        data[2] += rhs[2];
        return *this;
    }

    // unary minus
    RealVec operator-() const {
        const RealVec& lhs = *this;
        return RealVec(-lhs[0], -lhs[1], -lhs[2]);
    }

    // minus
    RealVec operator-(const RealVec& rhs) const {
        const RealVec& lhs = *this;
        return RealVec(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]);
    }

    RealVec& operator-=(const RealVec& rhs) {
        data[0] -= rhs[0];
        data[1] -= rhs[1];
        data[2] -= rhs[2];
        return *this;
    }

    // scalar product
    RealVec operator*(RealOpenMM rhs) const {
        const RealVec& lhs = *this;
        return RealVec(lhs[0]*rhs, lhs[1]*rhs, lhs[2]*rhs);
    }

    RealVec& operator*=(RealOpenMM rhs) {
        data[0] *= rhs;
        data[1] *= rhs;
        data[2] *= rhs;
        return *this;
    }

    // scalar division
    RealVec operator/(double rhs) const {
        const RealVec& lhs = *this;
        double scale = 1.0/rhs;
        return RealVec(lhs[0]*scale, lhs[1]*scale, lhs[2]*scale);
    }

    RealVec& operator/=(double rhs) {
        double scale = 1.0/rhs;
        data[0] *= scale;
        data[1] *= scale;
        data[2] *= scale;
        return *this;
    }

    // dot product
    RealOpenMM dot(const RealVec& rhs) const {
        const RealVec& lhs = *this;
        return lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2];
    }

    // cross product
    RealVec cross(const RealVec& rhs) const {
        return RealVec(data[1]*rhs[2]-data[2]*rhs[1], data[2]*rhs[0]-data[0]*rhs[2], data[0]*rhs[1]-data[1]*rhs[0]);
    }

private:
    RealOpenMM data[3];
};

template <class CHAR, class TRAITS>
std::basic_ostream<CHAR,TRAITS>& operator<<(std::basic_ostream<CHAR,TRAITS>& o, const RealVec& v) {
    o<<'['<<v[0]<<", "<<v[1]<<", "<<v[2]<<']';
    return o;
}

} // namespace OpenMM

#endif /*OPENMM_REALVEC_H_*/

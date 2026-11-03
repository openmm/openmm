#ifndef OPENMM_FIXED_H_
#define OPENMM_FIXED_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2026 Stanford University and the Authors.           *
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

namespace OpenMM {

/**
 * This class represents a 64 bit fixed point value.  The template parameter is the number of bits to use for the
 * fractional part.  This class is used to accumulate sums in a way that is deterministic and independent of order.
 * It includes only the features that are needed for that purpose.
 */
template <int BITS>
class Fixed {
public:
    Fixed() : value(0) {
    }
    Fixed(double d) : value((int64_t) (d*(1L<<BITS))) {
    }
    operator double() const {
        return value/(double) (1L<<BITS);
    }
    Fixed<BITS> operator+(Fixed<BITS> v) {
        return Fixed<BITS>(value+v.value);
    }
    Fixed<BITS> operator-(Fixed<BITS> v) {
        return Fixed<BITS>(value-v.value);
    }
    Fixed<BITS>& operator+=(Fixed<BITS> v) {
        value += v.value;
        return *this;
    }
    Fixed<BITS>& operator-=(Fixed<BITS> v) {
        value -= v.value;
        return *this;
    }
private:
    int64_t value;
    Fixed(int64_t i) : value(i) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_FIXED_H_*/

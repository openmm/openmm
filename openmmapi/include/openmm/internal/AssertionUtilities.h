#ifndef OPENMM_ASSERTIONUTILITIES_H_
#define OPENMM_ASSERTIONUTILITIES_H_

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

/**
 * This file provides a variety of macros useful in test cases.
 */

#include "windowsExport.h"
#include <cmath>
#include <sstream>

namespace OpenMM {

void OPENMM_EXPORT throwException(const char* file, int line, const std::string& details);

} // namespace OpenMM

#define ASSERT(cond) {if (!(cond)) throwException(__FILE__, __LINE__, "");};

#define ASSERT_EQUAL(expected, found) {if (!((expected) == (found))) {std::stringstream details; details << "Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_CONTAINERS(expected, found) {if (!((expected) == (found))) {std::stringstream details; details << "Containers not equal"; throwException(__FILE__, __LINE__, details.str());}};
 
#define ASSERT_EQUAL_TOL(expected, found, tol) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << "Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC(expected, found, tol) {double _norm_ = std::sqrt((expected).dot(expected)); double _scale_ = _norm_ > 1.0 ? _norm_ : 1.0; if ((std::abs(((expected)[0])-((found)[0]))/_scale_ > (tol)) || (std::abs(((expected)[1])-((found)[1]))/_scale_ > (tol)) || (std::abs(((expected)[2])-((found)[2]))/_scale_ > (tol))) {std::stringstream details; details << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_USUALLY_TRUE(cond) {if (!(cond)) throwException(__FILE__, __LINE__, "(This test is stochastic and may occasionally fail)");};

#define ASSERT_USUALLY_EQUAL_TOL(expected, found, tol) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << "Expected "<<(expected)<<", found "<<(found)<<" (This test is stochastic and may occasionally fail)"; throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_VALID_INDEX(index, vector) {if (index < 0 || index >= (int) vector.size()) throwException(__FILE__, __LINE__, "Index out of range");};

#endif /*OPENMM_ASSERTIONUTILITIES_H_*/

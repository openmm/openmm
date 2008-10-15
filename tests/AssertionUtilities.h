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
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "OpenMMException.h"
#include <cmath>
#include <string>
#include <sstream>

void throwException(const char* file, int line, const std::string& details) {
    std::string fn(file);
    std::string::size_type pos = fn.find_last_of("/\\");
    if (pos+1>=fn.size())
        pos=0;
    std::string filename(fn,(int)(pos+1),(int)(fn.size()-(pos+1)));
    std::stringstream message;
    message << "Assertion failure at "<<filename<<":"<<line;
    if (details.size() > 0)
        message << ".  "<<details;
    throw OpenMM::OpenMMException(message.str());
}

#define ASSERT(cond) {if (!(cond)) throwException(__FILE__, __LINE__, "");};

#define ASSERT_EQUAL(expected, found) {if (!((expected) == (found))) {std::stringstream details; details << "Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_TOL(expected, found, tol) {double _scale_ = std::fabs(expected) > 1.0 ? std::fabs(expected) : 1.0; if (!(std::fabs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << "Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC(expected, found, tol) {ASSERT_EQUAL_TOL((expected)[0], (found)[0], (tol)); ASSERT_EQUAL_TOL((expected)[1], (found)[1], (tol)); ASSERT_EQUAL_TOL((expected)[2], (found)[2], (tol));};

#endif /*OPENMM_ASSERTIONUTILITIES_H_*/

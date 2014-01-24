/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

#include "openmm/TabulatedFunction.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;
using namespace std;

Continuous1DFunction::Continuous1DFunction(const std::vector<double>& values, double min, double max) {
    if (max <= min)
        throw OpenMMException("Continuous1DFunction: max <= min for a tabulated function.");
    if (values.size() < 2)
        throw OpenMMException("Continuous1DFunction: a tabulated function must have at least two points");
    this->values = values;
    this->min = min;
    this->max = max;
}

void Continuous1DFunction::getFunctionParameters(std::vector<double>& values, double& min, double& max) const {
    values = this->values;
    min = this->min;
    max = this->max;
}

void Continuous1DFunction::setFunctionParameters(const std::vector<double>& values, double min, double max) {
    if (max <= min)
        throw OpenMMException("Continuous1DFunction: max <= min for a tabulated function.");
    if (values.size() < 2)
        throw OpenMMException("Continuous1DFunction: a tabulated function must have at least two points");
    this->values = values;
    this->min = min;
    this->max = max;
}

Discrete1DFunction::Discrete1DFunction(const std::vector<double>& values) {
    this->values = values;
}

void Discrete1DFunction::getFunctionParameters(std::vector<double>& values) const {
    values = this->values;
}

void Discrete1DFunction::setFunctionParameters(const std::vector<double>& values) {
    this->values = values;
}

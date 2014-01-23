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

#include "ReferenceTabulatedFunction.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/SplineFitter.h"

using namespace OpenMM;
using namespace std;
using Lepton::CustomFunction;

extern "C" CustomFunction* createReferenceTabulatedFunction(const TabulatedFunction& function) {
    if (dynamic_cast<const Continuous1DFunction*>(&function) != NULL)
        return new ReferenceContinuous1DFunction(dynamic_cast<const Continuous1DFunction&>(function));
    throw OpenMMException("createReferenceTabulatedFunction: Unknown function type");
}

ReferenceContinuous1DFunction::ReferenceContinuous1DFunction(const Continuous1DFunction& function) : function(function) {
    function.getFunctionParameters(values, min, max);
    int numValues = values.size();
    x.resize(numValues);
    for (int i = 0; i < numValues; i++)
        x[i] = min+i*(max-min)/(numValues-1);
    SplineFitter::createNaturalSpline(x, values, derivs);
}

int ReferenceContinuous1DFunction::getNumArguments() const {
    return 1;
}

double ReferenceContinuous1DFunction::evaluate(const double* arguments) const {
    double t = arguments[0];
    if (t < min || t > max)
        return 0.0;
    return SplineFitter::evaluateSpline(x, values, derivs, t);
}

double ReferenceContinuous1DFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    double t = arguments[0];
    if (t < min || t > max)
        return 0.0;
    return SplineFitter::evaluateSplineDerivative(x, values, derivs, t);
}

CustomFunction* ReferenceContinuous1DFunction::clone() const {
    return new ReferenceContinuous1DFunction(function);
}

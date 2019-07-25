/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014-2016 Stanford University and the Authors.      *
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

#ifdef _MSC_VER

#if _MSC_VER < 1800
/**
 * We need to define this ourselves, since Visual Studio is missing round() from cmath.
 */
static int round(double x) {
    return (int) (x+0.5);
}
#else
#include <cmath>
#endif  // MSC_VER < 1800


#else
#include <cmath>
#endif

using namespace OpenMM;
using namespace std;
using Lepton::CustomFunction;

extern "C" OPENMM_EXPORT CustomFunction* createReferenceTabulatedFunction(const TabulatedFunction& function) {
    if (dynamic_cast<const Continuous1DFunction*>(&function) != NULL)
        return new ReferenceContinuous1DFunction(dynamic_cast<const Continuous1DFunction&>(function));
    if (dynamic_cast<const Continuous2DFunction*>(&function) != NULL)
        return new ReferenceContinuous2DFunction(dynamic_cast<const Continuous2DFunction&>(function));
    if (dynamic_cast<const Continuous3DFunction*>(&function) != NULL)
        return new ReferenceContinuous3DFunction(dynamic_cast<const Continuous3DFunction&>(function));
    if (dynamic_cast<const Discrete1DFunction*>(&function) != NULL)
        return new ReferenceDiscrete1DFunction(dynamic_cast<const Discrete1DFunction&>(function));
    if (dynamic_cast<const Discrete2DFunction*>(&function) != NULL)
        return new ReferenceDiscrete2DFunction(dynamic_cast<const Discrete2DFunction&>(function));
    if (dynamic_cast<const Discrete3DFunction*>(&function) != NULL)
        return new ReferenceDiscrete3DFunction(dynamic_cast<const Discrete3DFunction&>(function));
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

ReferenceContinuous1DFunction::ReferenceContinuous1DFunction(const ReferenceContinuous1DFunction& other) : function(other.function) {
    function.getFunctionParameters(values, min, max);
    x = other.x;
    values = other.values;
    derivs = other.derivs;
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
    return new ReferenceContinuous1DFunction(*this);
}

ReferenceContinuous2DFunction::ReferenceContinuous2DFunction(const Continuous2DFunction& function) : function(function) {
    function.getFunctionParameters(xsize, ysize, values, xmin, xmax, ymin, ymax);
    x.resize(xsize);
    y.resize(ysize);
    for (int i = 0; i < xsize; i++)
        x[i] = xmin+i*(xmax-xmin)/(xsize-1);
    for (int i = 0; i < ysize; i++)
        y[i] = ymin+i*(ymax-ymin)/(ysize-1);
    SplineFitter::create2DNaturalSpline(x, y, values, c);
}

ReferenceContinuous2DFunction::ReferenceContinuous2DFunction(const ReferenceContinuous2DFunction& other) : function(other.function) {
    function.getFunctionParameters(xsize, ysize, values, xmin, xmax, ymin, ymax);
    x = other.x;
    y = other.y;
    values = other.values;
    c = other.c;
}

int ReferenceContinuous2DFunction::getNumArguments() const {
    return 2;
}

double ReferenceContinuous2DFunction::evaluate(const double* arguments) const {
    double u = arguments[0];
    if (u < xmin || u > xmax)
        return 0.0;
    double v = arguments[1];
    if (v < ymin || v > ymax)
        return 0.0;
    return SplineFitter::evaluate2DSpline(x, y, values, c, u, v);
}

double ReferenceContinuous2DFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    double u = arguments[0];
    if (u < xmin || u > xmax)
        return 0.0;
    double v = arguments[1];
    if (v < ymin || v > ymax)
        return 0.0;
    double dx, dy;
    SplineFitter::evaluate2DSplineDerivatives(x, y, values, c, u, v, dx, dy);
    if (derivOrder[0] == 1 && derivOrder[1] == 0)
        return dx;
    if (derivOrder[0] == 0 && derivOrder[1] == 1)
        return dy;
    throw OpenMMException("ReferenceContinuous2DFunction: Unsupported derivative order");
}

CustomFunction* ReferenceContinuous2DFunction::clone() const {
    return new ReferenceContinuous2DFunction(*this);
}

ReferenceContinuous3DFunction::ReferenceContinuous3DFunction(const Continuous3DFunction& function) : function(function) {
    function.getFunctionParameters(xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax);
    x.resize(xsize);
    y.resize(ysize);
    z.resize(zsize);
    for (int i = 0; i < xsize; i++)
        x[i] = xmin+i*(xmax-xmin)/(xsize-1);
    for (int i = 0; i < ysize; i++)
        y[i] = ymin+i*(ymax-ymin)/(ysize-1);
    for (int i = 0; i < zsize; i++)
        z[i] = zmin+i*(zmax-zmin)/(zsize-1);
    SplineFitter::create3DNaturalSpline(x, y, z, values, c);
}

ReferenceContinuous3DFunction::ReferenceContinuous3DFunction(const ReferenceContinuous3DFunction& other) : function(other.function) {
    function.getFunctionParameters(xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax);
    x = other.x;
    y = other.y;
    z = other.z;
    values = other.values;
    c = other.c;
}

int ReferenceContinuous3DFunction::getNumArguments() const {
    return 3;
}

double ReferenceContinuous3DFunction::evaluate(const double* arguments) const {
    double u = arguments[0];
    if (u < xmin || u > xmax)
        return 0.0;
    double v = arguments[1];
    if (v < ymin || v > ymax)
        return 0.0;
    double w = arguments[2];
    if (w < zmin || w > zmax)
        return 0.0;
    return SplineFitter::evaluate3DSpline(x, y, z, values, c, u, v, w);
}

double ReferenceContinuous3DFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    double u = arguments[0];
    if (u < xmin || u > xmax)
        return 0.0;
    double v = arguments[1];
    if (v < ymin || v > ymax)
        return 0.0;
    double w = arguments[2];
    if (w < zmin || w > zmax)
        return 0.0;
    double dx, dy, dz;
    SplineFitter::evaluate3DSplineDerivatives(x, y, z, values, c, u, v, w, dx, dy, dz);
    if (derivOrder[0] == 1 && derivOrder[1] == 0 && derivOrder[2] == 0)
        return dx;
    if (derivOrder[0] == 0 && derivOrder[1] == 1 && derivOrder[2] == 0)
        return dy;
    if (derivOrder[0] == 0 && derivOrder[1] == 0 && derivOrder[2] == 1)
        return dz;
    throw OpenMMException("ReferenceContinuous3DFunction: Unsupported derivative order");
}

CustomFunction* ReferenceContinuous3DFunction::clone() const {
    return new ReferenceContinuous3DFunction(*this);
}

ReferenceDiscrete1DFunction::ReferenceDiscrete1DFunction(const Discrete1DFunction& function) : function(function) {
    function.getFunctionParameters(values);
}

int ReferenceDiscrete1DFunction::getNumArguments() const {
    return 1;
}

double ReferenceDiscrete1DFunction::evaluate(const double* arguments) const {
    int i = (int) round(arguments[0]);
    if (i < 0 || i >= values.size())
        throw OpenMMException("ReferenceDiscrete1DFunction: argument out of range");
    return values[i];
}

double ReferenceDiscrete1DFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    return 0.0;
}

CustomFunction* ReferenceDiscrete1DFunction::clone() const {
    return new ReferenceDiscrete1DFunction(function);
}

ReferenceDiscrete2DFunction::ReferenceDiscrete2DFunction(const Discrete2DFunction& function) : function(function) {
    function.getFunctionParameters(xsize, ysize, values);
}

int ReferenceDiscrete2DFunction::getNumArguments() const {
    return 2;
}

double ReferenceDiscrete2DFunction::evaluate(const double* arguments) const {
    int i = (int) round(arguments[0]);
    int j = (int) round(arguments[1]);
    if (i < 0 || i >= xsize || j < 0 || j >= ysize)
        throw OpenMMException("ReferenceDiscrete2DFunction: argument out of range");
    return values[i+j*xsize];
}

double ReferenceDiscrete2DFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    return 0.0;
}

CustomFunction* ReferenceDiscrete2DFunction::clone() const {
    return new ReferenceDiscrete2DFunction(function);
}

ReferenceDiscrete3DFunction::ReferenceDiscrete3DFunction(const Discrete3DFunction& function) : function(function) {
    function.getFunctionParameters(xsize, ysize, zsize, values);
}

int ReferenceDiscrete3DFunction::getNumArguments() const {
    return 3;
}

double ReferenceDiscrete3DFunction::evaluate(const double* arguments) const {
    int i = (int) round(arguments[0]);
    int j = (int) round(arguments[1]);
    int k = (int) round(arguments[2]);
    if (i < 0 || i >= xsize || j < 0 || j >= ysize || k < 0 || k >= zsize)
        throw OpenMMException("ReferenceDiscrete3DFunction: argument out of range");
    return values[i+(j+k*ysize)*xsize];
}

double ReferenceDiscrete3DFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    return 0.0;
}

CustomFunction* ReferenceDiscrete3DFunction::clone() const {
    return new ReferenceDiscrete3DFunction(function);
}

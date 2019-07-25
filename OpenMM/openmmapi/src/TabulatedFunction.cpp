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

Continuous1DFunction::Continuous1DFunction(const vector<double>& values, double min, double max) {
    if (max <= min)
        throw OpenMMException("Continuous1DFunction: max <= min for a tabulated function.");
    if (values.size() < 2)
        throw OpenMMException("Continuous1DFunction: a tabulated function must have at least two points");
    this->values = values;
    this->min = min;
    this->max = max;
}

void Continuous1DFunction::getFunctionParameters(vector<double>& values, double& min, double& max) const {
    values = this->values;
    min = this->min;
    max = this->max;
}

void Continuous1DFunction::setFunctionParameters(const vector<double>& values, double min, double max) {
    if (max <= min)
        throw OpenMMException("Continuous1DFunction: max <= min for a tabulated function.");
    if (values.size() < 2)
        throw OpenMMException("Continuous1DFunction: a tabulated function must have at least two points");
    this->values = values;
    this->min = min;
    this->max = max;
}

Continuous1DFunction* Continuous1DFunction::Copy() const {
    vector<double> new_vec(values.size());
    for (size_t i = 0; i < values.size(); i++)
        new_vec[i] = values[i];
    return new Continuous1DFunction(new_vec, min, max);
}

Continuous2DFunction::Continuous2DFunction(int xsize, int ysize, const vector<double>& values, double xmin, double xmax, double ymin, double ymax) {
    if (xsize < 2 || ysize < 2)
        throw OpenMMException("Continuous2DFunction: must have at least two points along each axis");
    if (values.size() != xsize*ysize)
        throw OpenMMException("Continuous2DFunction: incorrect number of values");
    if (xmax <= xmin)
        throw OpenMMException("Continuous2DFunction: xmax <= xmin for a tabulated function.");
    if (ymax <= ymin)
        throw OpenMMException("Continuous2DFunction: ymax <= ymin for a tabulated function.");
    this->values = values;
    this->xsize = xsize;
    this->ysize = ysize;
    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
}

void Continuous2DFunction::getFunctionParameters(int& xsize, int& ysize, vector<double>& values, double& xmin, double& xmax, double& ymin, double& ymax) const {
    values = this->values;
    xsize = this->xsize;
    ysize = this->ysize;
    xmin = this->xmin;
    xmax = this->xmax;
    ymin = this->ymin;
    ymax = this->ymax;
}

void Continuous2DFunction::setFunctionParameters(int xsize, int ysize, const vector<double>& values, double xmin, double xmax, double ymin, double ymax) {
    if (xsize < 2 || ysize < 2)
        throw OpenMMException("Continuous2DFunction: must have at least two points along each axis");
    if (values.size() != xsize*ysize)
        throw OpenMMException("Continuous2DFunction: incorrect number of values");
    if (xmax <= xmin)
        throw OpenMMException("Continuous2DFunction: xmax <= xmin for a tabulated function.");
    if (ymax <= ymin)
        throw OpenMMException("Continuous2DFunction: ymax <= ymin for a tabulated function.");
    this->values = values;
    this->xsize = xsize;
    this->ysize = ysize;
    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
}

Continuous2DFunction* Continuous2DFunction::Copy() const {
    vector<double> new_vec(values.size());
    for (size_t i = 0; i < values.size(); i++)
        new_vec[i] = values[i];
    return new Continuous2DFunction(xsize, ysize, new_vec, xmin, xmax, ymin, ymax);
}

Continuous3DFunction::Continuous3DFunction(int xsize, int ysize, int zsize, const vector<double>& values, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
    if (xsize < 2 || ysize < 2 || zsize < 2)
        throw OpenMMException("Continuous3DFunction: must have at least two points along each axis");
    if (values.size() != xsize*ysize*zsize)
        throw OpenMMException("Continuous3DFunction: incorrect number of values");
    if (xmax <= xmin)
        throw OpenMMException("Continuous3DFunction: xmax <= xmin for a tabulated function.");
    if (ymax <= ymin)
        throw OpenMMException("Continuous3DFunction: ymax <= ymin for a tabulated function.");
    if (zmax <= zmin)
        throw OpenMMException("Continuous3DFunction: zmax <= zmin for a tabulated function.");
    this->values = values;
    this->xsize = xsize;
    this->ysize = ysize;
    this->zsize = zsize;
    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
    this->zmin = zmin;
    this->zmax = zmax;
}

void Continuous3DFunction::getFunctionParameters(int& xsize, int& ysize, int& zsize, vector<double>& values, double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) const {
    values = this->values;
    xsize = this->xsize;
    ysize = this->ysize;
    zsize = this->zsize;
    xmin = this->xmin;
    xmax = this->xmax;
    ymin = this->ymin;
    ymax = this->ymax;
    zmin = this->zmin;
    zmax = this->zmax;
}

void Continuous3DFunction::setFunctionParameters(int xsize, int ysize, int zsize, const vector<double>& values, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
    if (xsize < 2 || ysize < 2 || zsize < 2)
        throw OpenMMException("Continuous3DFunction: must have at least two points along each axis");
    if (values.size() != xsize*ysize*zsize)
        throw OpenMMException("Continuous3DFunction: incorrect number of values");
    if (xmax <= xmin)
        throw OpenMMException("Continuous3DFunction: xmax <= xmin for a tabulated function.");
    if (ymax <= ymin)
        throw OpenMMException("Continuous3DFunction: ymax <= ymin for a tabulated function.");
    if (zmax <= zmin)
        throw OpenMMException("Continuous3DFunction: zmax <= zmin for a tabulated function.");
    this->values = values;
    this->xsize = xsize;
    this->ysize = ysize;
    this->zsize = zsize;
    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
    this->zmin = zmin;
    this->zmax = zmax;
}

Continuous3DFunction* Continuous3DFunction::Copy() const {
    vector<double> new_vec(values.size());
    for (size_t i = 0; i < values.size(); i++)
        new_vec[i] = values[i];
    return new Continuous3DFunction(xsize, ysize, zsize, new_vec, xmin, xmax, ymin, ymax, zmin, zmax);
}


Discrete1DFunction::Discrete1DFunction(const vector<double>& values) {
    this->values = values;
}

void Discrete1DFunction::getFunctionParameters(vector<double>& values) const {
    values = this->values;
}

void Discrete1DFunction::setFunctionParameters(const vector<double>& values) {
    this->values = values;
}

Discrete1DFunction* Discrete1DFunction::Copy() const {
    vector<double> new_vec(values.size());
    for (size_t i = 0; i < values.size(); i++)
        new_vec[i] = values[i];
    return new Discrete1DFunction(new_vec);
}

Discrete2DFunction::Discrete2DFunction(int xsize, int ysize, const vector<double>& values) {
    if (values.size() != xsize*ysize)
        throw OpenMMException("Discrete2DFunction: incorrect number of values");
    this->xsize = xsize;
    this->ysize = ysize;
    this->values = values;
}

void Discrete2DFunction::getFunctionParameters(int& xsize, int& ysize, vector<double>& values) const {
    xsize = this->xsize;
    ysize = this->ysize;
    values = this->values;
}

void Discrete2DFunction::setFunctionParameters(int xsize, int ysize, const vector<double>& values) {
    if (values.size() != xsize*ysize)
        throw OpenMMException("Discrete2DFunction: incorrect number of values");
    this->xsize = xsize;
    this->ysize = ysize;
    this->values = values;
}

Discrete2DFunction* Discrete2DFunction::Copy() const {
    vector<double> new_vec(values.size());
    for (size_t i = 0; i < values.size(); i++)
        new_vec[i] = values[i];
    return new Discrete2DFunction(xsize, ysize, new_vec);
}

Discrete3DFunction::Discrete3DFunction(int xsize, int ysize, int zsize, const vector<double>& values) {
    if (values.size() != xsize*ysize*zsize)
        throw OpenMMException("Discrete3DFunction: incorrect number of values");
    this->xsize = xsize;
    this->ysize = ysize;
    this->zsize = zsize;
    this->values = values;
}

void Discrete3DFunction::getFunctionParameters(int& xsize, int& ysize, int& zsize, vector<double>& values) const {
    xsize = this->xsize;
    ysize = this->ysize;
    zsize = this->zsize;
    values = this->values;
}

void Discrete3DFunction::setFunctionParameters(int xsize, int ysize, int zsize, const vector<double>& values) {
    if (values.size() != xsize*ysize*zsize)
        throw OpenMMException("Discrete3DFunction: incorrect number of values");
    this->xsize = xsize;
    this->ysize = ysize;
    this->zsize = zsize;
    this->values = values;
}

Discrete3DFunction* Discrete3DFunction::Copy() const {
    vector<double> new_vec(values.size());
    for (size_t i = 0; i < values.size(); i++)
        new_vec[i] = values[i];
    return new Discrete3DFunction(xsize, ysize, zsize, new_vec);
}

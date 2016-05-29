#ifndef OPENMM_REFERENCETABULATEDFUNCTION_H_
#define OPENMM_REFERENCETABULATEDFUNCTION_H_

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

#include "openmm/TabulatedFunction.h"
#include "openmm/internal/windowsExport.h"
#include "lepton/CustomFunction.h"
#include <vector>

namespace OpenMM {

/**
 * Given a TabulatedFunction, wrap it in an appropriate subclass of Lepton::CustomFunction.
 */
extern "C" OPENMM_EXPORT Lepton::CustomFunction* createReferenceTabulatedFunction(const TabulatedFunction& function);

/**
 * This class adapts a Continuous1DFunction into a Lepton::CustomFunction.
 */
class OPENMM_EXPORT ReferenceContinuous1DFunction : public Lepton::CustomFunction {
public:
    ReferenceContinuous1DFunction(const Continuous1DFunction& function);
    int getNumArguments() const;
    double evaluate(const double* arguments) const;
    double evaluateDerivative(const double* arguments, const int* derivOrder) const;
    CustomFunction* clone() const;
private:
    ReferenceContinuous1DFunction(const ReferenceContinuous1DFunction& other);
    const Continuous1DFunction& function;
    double min, max;
    std::vector<double> x, values, derivs;
};

/**
 * This class adapts a Continuous2DFunction into a Lepton::CustomFunction.
 */
class OPENMM_EXPORT ReferenceContinuous2DFunction : public Lepton::CustomFunction {
public:
    ReferenceContinuous2DFunction(const Continuous2DFunction& function);
    int getNumArguments() const;
    double evaluate(const double* arguments) const;
    double evaluateDerivative(const double* arguments, const int* derivOrder) const;
    CustomFunction* clone() const;
private:
    ReferenceContinuous2DFunction(const ReferenceContinuous2DFunction& other);
    const Continuous2DFunction& function;
    int xsize, ysize;
    double xmin, xmax, ymin, ymax;
    std::vector<double> x, y, values;
    std::vector<std::vector<double> > c;
};

/**
 * This class adapts a Continuous3DFunction into a Lepton::CustomFunction.
 */
class OPENMM_EXPORT ReferenceContinuous3DFunction : public Lepton::CustomFunction {
public:
    ReferenceContinuous3DFunction(const Continuous3DFunction& function);
    int getNumArguments() const;
    double evaluate(const double* arguments) const;
    double evaluateDerivative(const double* arguments, const int* derivOrder) const;
    CustomFunction* clone() const;
private:
    ReferenceContinuous3DFunction(const ReferenceContinuous3DFunction& other);
    const Continuous3DFunction& function;
    int xsize, ysize, zsize;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    std::vector<double> x, y, z, values;
    std::vector<std::vector<double> > c;
};

/**
 * This class adapts a Discrete1DFunction into a Lepton::CustomFunction.
 */
class OPENMM_EXPORT ReferenceDiscrete1DFunction : public Lepton::CustomFunction {
public:
    ReferenceDiscrete1DFunction(const Discrete1DFunction& function);
    int getNumArguments() const;
    double evaluate(const double* arguments) const;
    double evaluateDerivative(const double* arguments, const int* derivOrder) const;
    CustomFunction* clone() const;
private:
    const Discrete1DFunction& function;
    std::vector<double> values;
};

/**
 * This class adapts a Discrete2DFunction into a Lepton::CustomFunction.
 */
class OPENMM_EXPORT ReferenceDiscrete2DFunction : public Lepton::CustomFunction {
public:
    ReferenceDiscrete2DFunction(const Discrete2DFunction& function);
    int getNumArguments() const;
    double evaluate(const double* arguments) const;
    double evaluateDerivative(const double* arguments, const int* derivOrder) const;
    CustomFunction* clone() const;
private:
    const Discrete2DFunction& function;
    int xsize, ysize;
    std::vector<double> values;
};

/**
 * This class adapts a Discrete3DFunction into a Lepton::CustomFunction.
 */
class OPENMM_EXPORT ReferenceDiscrete3DFunction : public Lepton::CustomFunction {
public:
    ReferenceDiscrete3DFunction(const Discrete3DFunction& function);
    int getNumArguments() const;
    double evaluate(const double* arguments) const;
    double evaluateDerivative(const double* arguments, const int* derivOrder) const;
    CustomFunction* clone() const;
private:
    const Discrete3DFunction& function;
    int xsize, ysize, zsize;
    std::vector<double> values;
};

} // namespace OpenMM

#endif /*OPENMM_REFERENCETABULATEDFUNCTION_H_*/

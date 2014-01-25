#ifndef OPENMM_TABULATEDFUNCTION_H_
#define OPENMM_TABULATEDFUNCTION_H_

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

#include "internal/windowsExport.h"
#include <vector>

namespace OpenMM {

/**
 * A TabulatedFunction uses a set of tabulated values to define a mathematical function.
 * It can be used by various custom forces.
 * 
 * TabulatedFunction is an abstract class with concrete subclasses for more specific
 * types of functions.  There are subclasses for:
 * 
 * <ul>
 * <li>1, 2, and 3 dimensional functions.  The dimensionality of a function means
 * the number of input arguments it takes.</li>
 * <li>Continuous and discrete functions.  A continuous function is interpolated by
 * fitting a natural cubic spline to the tabulated values.  A discrete function is
 * only defined for integer values of its arguments (that is, at the tabulated points),
 * and does not try to interpolate between them.  Discrete function can be evaluated
 * more quickly than continuous ones.</li>
 * </ul>
 */

class OPENMM_EXPORT TabulatedFunction {
public:
    virtual ~TabulatedFunction() {
    }
};

/**
 * This is a TabulatedFunction that computes a continuous one dimensional function.
 */
class OPENMM_EXPORT Continuous1DFunction : public TabulatedFunction {
public:
    /**
     * Create a Continuous1DFunction f(x) based on a set of tabulated values.
     * 
     * @param values         the tabulated values of the function f(x) at uniformly spaced values of x between min
     *                       and max.  A natural cubic spline is used to interpolated between the tabulated values.
     *                       The function is assumed to be zero for x &lt; min or x &gt; max.
     * @param min            the value of x corresponding to the first element of values
     * @param max            the value of x corresponding to the last element of values
     */
    Continuous1DFunction(const std::vector<double>& values, double min, double max);
    /**
     * Get the parameters for the tabulated function.
     *
     * @param values         the tabulated values of the function f(x) at uniformly spaced values of x between min
     *                       and max.  A natural cubic spline is used to interpolated between the tabulated values.
     *                       The function is assumed to be zero for x &lt; min or x &gt; max.
     * @param min            the value of x corresponding to the first element of values
     * @param max            the value of x corresponding to the last element of values
     */
    void getFunctionParameters(std::vector<double>& values, double& min, double& max) const;
    /**
     * Set the parameters for the tabulated function.
     *
     * @param values         the tabulated values of the function f(x) at uniformly spaced values of x between min
     *                       and max.  A natural cubic spline is used to interpolated between the tabulated values.
     *                       The function is assumed to be zero for x &lt; min or x &gt; max.
     * @param min            the value of x corresponding to the first element of values
     * @param max            the value of x corresponding to the last element of values
     */
    void setFunctionParameters(const std::vector<double>& values, double min, double max);
private:
    std::vector<double> values;
    double min, max;
};

/**
 * This is a TabulatedFunction that computes a discrete one dimensional function f(x).
 * To evaluate it, x is rounded to the nearest integer and the table element with that
 * index is returned.  If the index is outside the range [0, size), the result is undefined.
 */
class OPENMM_EXPORT Discrete1DFunction : public TabulatedFunction {
public:
    /**
     * Create a Discrete1DFunction f(x) based on a set of tabulated values.
     * 
     * @param values         the tabulated values of the function f(x)
     */
    Discrete1DFunction(const std::vector<double>& values);
    /**
     * Get the parameters for the tabulated function.
     *
     * @param values         the tabulated values of the function f(x)
     */
    void getFunctionParameters(std::vector<double>& values) const;
    /**
     * Set the parameters for the tabulated function.
     *
     * @param values         the tabulated values of the function f(x)
     */
    void setFunctionParameters(const std::vector<double>& values);
private:
    std::vector<double> values;
};

/**
 * This is a TabulatedFunction that computes a discrete two dimensional function f(x,y).
 * To evaluate it, x and y are each rounded to the nearest integer and the table element with those
 * indices is returned.  If either index is outside the range [0, size), the result is undefined.
 */
class OPENMM_EXPORT Discrete2DFunction : public TabulatedFunction {
public:
    /**
     * Create a Discrete2DFunction f(x,y) based on a set of tabulated values.
     * 
     * @param xsize     the number of table elements along the x direction
     * @param ysize     the number of table elements along the y direction
     * @param values    the tabulated values of the function f(x,y), ordered so that
     *                  values[i+xsize*j] = f(i,j).  This must be of length xsize*ysize.
     */
    Discrete2DFunction(int xsize, int ysize, const std::vector<double>& values);
    /**
     * Get the parameters for the tabulated function.
     *
     * @param xsize     the number of table elements along the x direction
     * @param ysize     the number of table elements along the y direction
     * @param values    the tabulated values of the function f(x,y), ordered so that
     *                  values[i+xsize*j] = f(i,j).  This must be of length xsize*ysize.
     */
    void getFunctionParameters(int& xsize, int& ysize, std::vector<double>& values) const;
    /**
     * Set the parameters for the tabulated function.
     *
     * @param xsize     the number of table elements along the x direction
     * @param ysize     the number of table elements along the y direction
     * @param values    the tabulated values of the function f(x,y), ordered so that
     *                  values[i+xsize*j] = f(i,j).  This must be of length xsize*ysize.
     */
    void setFunctionParameters(int xsize, int ysize, const std::vector<double>& values);
private:
    int xsize, ysize;
    std::vector<double> values;
};

/**
 * This is a TabulatedFunction that computes a discrete three dimensional function f(x,y,z).
 * To evaluate it, x, y, and z are each rounded to the nearest integer and the table element with those
 * indices is returned.  If any index is outside the range [0, size), the result is undefined.
 */
class OPENMM_EXPORT Discrete3DFunction : public TabulatedFunction {
public:
    /**
     * Create a Discrete3DFunction f(x,y,z) based on a set of tabulated values.
     * 
     * @param xsize     the number of table elements along the x direction
     * @param ysize     the number of table elements along the y direction
     * @param zsize     the number of table elements along the z direction
     * @param values    the tabulated values of the function f(x,y,z), ordered so that
     *                  values[i+xsize*j+xsize*ysize*k] = f(i,j,k).  This must be of length xsize*ysize*zsize.
     */
    Discrete3DFunction(int xsize, int ysize, int zsize, const std::vector<double>& values);
    /**
     * Get the parameters for the tabulated function.
     *
     * @param xsize     the number of table elements along the x direction
     * @param ysize     the number of table elements along the y direction
     * @param zsize     the number of table elements along the z direction
     * @param values    the tabulated values of the function f(x,y,z), ordered so that
     *                  values[i+xsize*j+xsize*ysize*k] = f(i,j,k).  This must be of length xsize*ysize*zsize.
     */
    void getFunctionParameters(int& xsize, int& ysize, int& zsize, std::vector<double>& values) const;
    /**
     * Set the parameters for the tabulated function.
     *
     * @param xsize     the number of table elements along the x direction
     * @param ysize     the number of table elements along the y direction
     * @param zsize     the number of table elements along the z direction
     * @param values    the tabulated values of the function f(x,y,z), ordered so that
     *                  values[i+xsize*j+xsize*ysize*k] = f(i,j,k).  This must be of length xsize*ysize*zsize.
     */
    void setFunctionParameters(int xsize, int ysize, int zsize, const std::vector<double>& values);
private:
    int xsize, ysize, zsize;
    std::vector<double> values;
};

} // namespace OpenMM

#endif /*OPENMM_TABULATEDFUNCTION_H_*/

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
    /**
     * @deprecated This will be removed in a future release.
     */
    virtual TabulatedFunction* Copy() const = 0;
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
     *                       and max.  A natural cubic spline is used to interpolate between the tabulated values.
     *                       The function is assumed to be zero for x &lt; min or x &gt; max.
     * @param min            the value of x corresponding to the first element of values
     * @param max            the value of x corresponding to the last element of values
     */
    Continuous1DFunction(const std::vector<double>& values, double min, double max);
    /**
     * Get the parameters for the tabulated function.
     *
     * @param[out] values         the tabulated values of the function f(x) at uniformly spaced values of x between min
     *                            and max.  A natural cubic spline is used to interpolate between the tabulated values.
     *                            The function is assumed to be zero for x &lt; min or x &gt; max.
     * @param[out] min            the value of x corresponding to the first element of values
     * @param[out] max            the value of x corresponding to the last element of values
     */
    void getFunctionParameters(std::vector<double>& values, double& min, double& max) const;
    /**
     * Set the parameters for the tabulated function.
     *
     * @param values         the tabulated values of the function f(x) at uniformly spaced values of x between min
     *                       and max.  A natural cubic spline is used to interpolate between the tabulated values.
     *                       The function is assumed to be zero for x &lt; min or x &gt; max.
     * @param min            the value of x corresponding to the first element of values
     * @param max            the value of x corresponding to the last element of values
     */
    void setFunctionParameters(const std::vector<double>& values, double min, double max);
    /**
     * Create a deep copy of the tabulated function.
     * 
     * @deprecated This will be removed in a future release.
     */
    Continuous1DFunction* Copy() const;
private:
    std::vector<double> values;
    double min, max;
};

/**
 * This is a TabulatedFunction that computes a continuous two dimensional function.
 */
class OPENMM_EXPORT Continuous2DFunction : public TabulatedFunction {
public:
    /**
     * Create a Continuous2DFunction f(x,y) based on a set of tabulated values.
     *
     * @param values     the tabulated values of the function f(x,y) at xsize uniformly spaced values of x between xmin
     *                   and xmax, and ysize values of y between ymin and ymax.  A natural cubic spline is used to interpolate between the tabulated values.
     *                   The function is assumed to be zero when x or y is outside its specified range.  The values should be ordered so that
     *                   values[i+xsize*j] = f(x_i,y_j), where x_i is the i'th uniformly spaced value of x.  This must be of length xsize*ysize.
     * @param xsize      the number of table elements along the x direction
     * @param ysize      the number of table elements along the y direction
     * @param xmin       the value of x corresponding to the first element of values
     * @param xmax       the value of x corresponding to the last element of values
     * @param ymin       the value of y corresponding to the first element of values
     * @param ymax       the value of y corresponding to the last element of values
     */
    Continuous2DFunction(int xsize, int ysize, const std::vector<double>& values, double xmin, double xmax, double ymin, double ymax);
    /**
     * Get the parameters for the tabulated function.
     *
     * @param[out] values     the tabulated values of the function f(x,y) at xsize uniformly spaced values of x between xmin
     *                        and xmax, and ysize values of y between ymin and ymax.  A natural cubic spline is used to interpolate between the tabulated values.
     *                        The function is assumed to be zero when x or y is outside its specified range.  The values should be ordered so that
     *                        values[i+xsize*j] = f(x_i,y_j), where x_i is the i'th uniformly spaced value of x.  This must be of length xsize*ysize.
     * @param[out] xsize      the number of table elements along the x direction
     * @param[out] ysize      the number of table elements along the y direction
     * @param[out] xmin       the value of x corresponding to the first element of values
     * @param[out] xmax       the value of x corresponding to the last element of values
     * @param[out] ymin       the value of y corresponding to the first element of values
     * @param[out] ymax       the value of y corresponding to the last element of values
     */
    void getFunctionParameters(int& xsize, int& ysize, std::vector<double>& values, double& xmin, double& xmax, double& ymin, double& ymax) const;
    /**
     * Set the parameters for the tabulated function.
     *
     * @param values     the tabulated values of the function f(x,y) at xsize uniformly spaced values of x between xmin
     *                   and xmax, and ysize values of y between ymin and ymax.  A natural cubic spline is used to interpolate between the tabulated values.
     *                   The function is assumed to be zero when x or y is outside its specified range.  The values should be ordered so that
     *                   values[i+xsize*j] = f(x_i,y_j), where x_i is the i'th uniformly spaced value of x.  This must be of length xsize*ysize.
     * @param xsize      the number of table elements along the x direction
     * @param ysize      the number of table elements along the y direction
     * @param xmin       the value of x corresponding to the first element of values
     * @param xmax       the value of x corresponding to the last element of values
     * @param ymin       the value of y corresponding to the first element of values
     * @param ymax       the value of y corresponding to the last element of values
     */
    void setFunctionParameters(int xsize, int ysize, const std::vector<double>& values, double xmin, double xmax, double ymin, double ymax);
    /**
     * Create a deep copy of the tabulated function
     * 
     * @deprecated This will be removed in a future release.
     */
    Continuous2DFunction* Copy() const;
private:
    std::vector<double> values;
    int xsize, ysize;
    double xmin, xmax, ymin, ymax;
};

/**
 * This is a TabulatedFunction that computes a continuous three dimensional function.
 */
class OPENMM_EXPORT Continuous3DFunction : public TabulatedFunction {
public:
    /**
     * Create a Continuous3DFunction f(x,y,z) based on a set of tabulated values.
     *
     * @param values     the tabulated values of the function f(x,y,z) at xsize uniformly spaced values of x between xmin
     *                   and xmax, ysize values of y between ymin and ymax, and zsize values of z between zmin and zmax.
     *                   A natural cubic spline is used to interpolate between the tabulated values.  The function is
     *                   assumed to be zero when x, y, or z is outside its specified range.  The values should be ordered so
     *                   that values[i+xsize*j+xsize*ysize*k] = f(x_i,y_j,z_k), where x_i is the i'th uniformly spaced value of x.
     *                   This must be of length xsize*ysize*zsize.
     * @param xsize      the number of table elements along the x direction
     * @param ysize      the number of table elements along the y direction
     * @param ysize      the number of table elements along the z direction
     * @param xmin       the value of x corresponding to the first element of values
     * @param xmax       the value of x corresponding to the last element of values
     * @param ymin       the value of y corresponding to the first element of values
     * @param ymax       the value of y corresponding to the last element of values
     * @param zmin       the value of z corresponding to the first element of values
     * @param zmax       the value of z corresponding to the last element of values
     */
    Continuous3DFunction(int xsize, int ysize, int zsize, const std::vector<double>& values, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
    /**
     * Get the parameters for the tabulated function.
     *
     * @param[out] values     the tabulated values of the function f(x,y,z) at xsize uniformly spaced values of x between xmin
     *                        and xmax, ysize values of y between ymin and ymax, and zsize values of z between zmin and zmax.
     *                        A natural cubic spline is used to interpolate between the tabulated values.  The function is
     *                        assumed to be zero when x, y, or z is outside its specified range.  The values should be ordered so
     *                        that values[i+xsize*j+xsize*ysize*k] = f(x_i,y_j,z_k), where x_i is the i'th uniformly spaced value of x.
     *                        This must be of length xsize*ysize*zsize.
     * @param[out] xsize      the number of table elements along the x direction
     * @param[out] ysize      the number of table elements along the y direction
     * @param[out] zsize      the number of table elements along the z direction
     * @param[out] xmin       the value of x corresponding to the first element of values
     * @param[out] xmax       the value of x corresponding to the last element of values
     * @param[out] ymin       the value of y corresponding to the first element of values
     * @param[out] ymax       the value of y corresponding to the last element of values
     * @param[out] zmin       the value of z corresponding to the first element of values
     * @param[out] zmax       the value of z corresponding to the last element of values
     */
    void getFunctionParameters(int& xsize, int& ysize, int& zsize, std::vector<double>& values, double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) const;
    /**
     * Set the parameters for the tabulated function.
     *
     * @param values     the tabulated values of the function f(x,y,z) at xsize uniformly spaced values of x between xmin
     *                   and xmax, ysize values of y between ymin and ymax, and zsize values of z between zmin and zmax.
     *                   A natural cubic spline is used to interpolate between the tabulated values.  The function is
     *                   assumed to be zero when x, y, or z is outside its specified range.  The values should be ordered so
     *                   that values[i+xsize*j+xsize*ysize*k] = f(x_i,y_j,z_k), where x_i is the i'th uniformly spaced value of x.
     *                   This must be of length xsize*ysize*zsize.
     * @param xsize      the number of table elements along the x direction
     * @param ysize      the number of table elements along the y direction
     * @param zsize      the number of table elements along the z direction
     * @param xmin       the value of x corresponding to the first element of values
     * @param xmax       the value of x corresponding to the last element of values
     * @param ymin       the value of y corresponding to the first element of values
     * @param ymax       the value of y corresponding to the last element of values
     * @param zmin       the value of z corresponding to the first element of values
     * @param zmax       the value of z corresponding to the last element of values
     */
    void setFunctionParameters(int xsize, int ysize, int zsize, const std::vector<double>& values, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
    /**
     * Create a deep copy of the tabulated function
     * 
     * @deprecated This will be removed in a future release.
     */
    Continuous3DFunction* Copy() const;
private:
    std::vector<double> values;
    int xsize, ysize, zsize;
    double xmin, xmax, ymin, ymax, zmin, zmax;
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
     * @param[out] values    the tabulated values of the function f(x)
     */
    void getFunctionParameters(std::vector<double>& values) const;
    /**
     * Set the parameters for the tabulated function.
     *
     * @param values         the tabulated values of the function f(x)
     */
    void setFunctionParameters(const std::vector<double>& values);
    /**
     * Create a deep copy of the tabulated function
     * 
     * @deprecated This will be removed in a future release.
     */
    Discrete1DFunction* Copy() const;
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
     * @param[out] xsize     the number of table elements along the x direction
     * @param[out] ysize     the number of table elements along the y direction
     * @param[out] values    the tabulated values of the function f(x,y), ordered so that
     *                       values[i+xsize*j] = f(i,j).  This must be of length xsize*ysize.
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
    /**
     * Create a deep copy of the tabulated function
     * 
     * @deprecated This will be removed in a future release.
     */
    Discrete2DFunction* Copy() const;
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
     * @param[out] xsize     the number of table elements along the x direction
     * @param[out] ysize     the number of table elements along the y direction
     * @param[out] zsize     the number of table elements along the z direction
     * @param[out] values    the tabulated values of the function f(x,y,z), ordered so that
     *                       values[i+xsize*j+xsize*ysize*k] = f(i,j,k).  This must be of length xsize*ysize*zsize.
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
    /**
     * Create a deep copy of the tabulated function
     * 
     * @deprecated This will be removed in a future release.
     */
    Discrete3DFunction* Copy() const;
private:
    int xsize, ysize, zsize;
    std::vector<double> values;
};

} // namespace OpenMM

#endif /*OPENMM_TABULATEDFUNCTION_H_*/

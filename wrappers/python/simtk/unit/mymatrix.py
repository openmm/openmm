"""
Pure python inversion of small matrices, to avoid requiring numpy or similar in SimTK.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Christopher M. Bruns
Contributors: Peter Eastman

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import print_function, division, absolute_import

import sys

def eye(size):
    """
    Returns identity matrix.

    >>> print(eye(3))
    [[1, 0, 0]
     [0, 1, 0]
     [0, 0, 1]]
    """
    result = []
    for row in range(size):
        r = []
        for col in range(size):
            if row == col:
                r.append(1)
            else:
                r.append(0)
        result.append(r)
    return MyMatrix(result)

def zeros(m, n=None):
    """
    Returns matrix of zeroes

    >>> print(zeros(3))
    [[0, 0, 0]
     [0, 0, 0]
     [0, 0, 0]]
    """
    if n is None:
        n = m
    result = []
    for row in range(m):
        r = []
        for col in range(n):
            r.append(0)
        result.append(r)
    return MyMatrix(result)

class MyVector(object):
    """
    Parent class of MyMatrix and type of Matrix Row.
    """
    def __init__(self, collection):
        if isinstance(collection, MyVector):
            self.data = collection.data
        else:
            self.data = collection

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return self.__class__.__name__ + "(" + repr(self.data) + ")"

    def __getitem__(self, key):
        return self.data[key]

    def __contains__(self, item):
        return item in self.data

    def __delitem__(self, key):
        del self.data[key]

    def __iter__(self):
        for item in self.data:
            yield item

    def __len__(self):
        return len(self.data)

    def __setitem__(self, key, value):
        self.data[key] = value

    def __rmul__(self, lhs):
        try:
            len(lhs)
            # left side is not scalar, delegate mul to that class
            return NotImplemented
        except TypeError:
            new_vec = []
            for element in self:
                new_vec.append(lhs * element)
            return self.__class__(new_vec)

class MyMatrix(MyVector):
    """
    Pure python linear algebra matrix for internal matrix inversion in UnitSystem.

    >>> m = MyMatrix([[1,0,],[0,1,]])
    >>> print(m)
    [[1, 0]
     [0, 1]]
    >>> print(~m)
    [[1.0, 0.0]
     [0.0, 1.0]]
    >>> print(eye(5))
    [[1, 0, 0, 0, 0]
     [0, 1, 0, 0, 0]
     [0, 0, 1, 0, 0]
     [0, 0, 0, 1, 0]
     [0, 0, 0, 0, 1]]
    >>> m = eye(5)
    >>> m[1][1]
    1
    >>> m[1:4]
    MyMatrixTranspose([[0, 0, 0],[1, 0, 0],[0, 1, 0],[0, 0, 1],[0, 0, 0]])
    >>> print(m[1:4])
    [[0, 0, 0]
     [1, 0, 0]
     [0, 1, 0]
     [0, 0, 1]
     [0, 0, 0]]
    >>> print(m[1:4][0:2])
    [[0, 1]
     [0, 0]
     [0, 0]]
    >>> m[1:4][0:2] = [[9,8],[7,6],[5,4]]
    >>> print(m)
    [[1, 0, 0, 0, 0]
     [9, 8, 0, 0, 0]
     [7, 6, 1, 0, 0]
     [5, 4, 0, 1, 0]
     [0, 0, 0, 0, 1]]
    """
    def numRows(self):
        return len(self.data)

    def numCols(self):
        if len(self.data) == 0:
            return 0
        else:
            return len(self.data[0])

    def __len__(self):
        return self.numRows()

    def __str__(self):
        result = ""
        start_char = "["
        for m in range(self.numRows()):
            result += start_char
            result += str(self[m])
            if m < self.numRows() - 1:
                result += "\n"
            start_char = " "
        result += "]"
        return result

    def __repr__(self):
        return 'MyMatrix(' + MyVector.__repr__(self) + ')'

    def is_square(self):
        return self.numRows() == self.numCols()

    def __iter__(self):
        for item in self.data:
            yield MyVector(item)

    def __getitem__(self, m):
        if isinstance(m, slice):
            return MyMatrixTranspose(self.data[m])
        else:
            return MyVector(self.data[m])

    def __setitem__(self, key, rhs):
        if isinstance(key, slice):
            self.data[key] = rhs
        else:
            assert len(rhs) == self.numCols()
            self.data[key] = MyVector(rhs)

    def __mul__(self, rhs):
        """
        Matrix multiplication.

        >>> a = MyMatrix([[1,2],[3,4]])
        >>> b = MyMatrix([[5,6],[7,8]])
        >>> print(a)
        [[1, 2]
         [3, 4]]
        >>> print(b)
        [[5, 6]
         [7, 8]]
        >>> print(a*b)
        [[19, 22]
         [43, 50]]

        """
        m = self.numRows()
        n = len(rhs[0])
        r = len(rhs)
        if self.numCols() != r:
            raise ArithmeticError("Matrix multplication size mismatch (%d vs %d)" % (self.numCols(), r))
        result = zeros(m, n)
        for i in range(m):
            for j in range(n):
                for k in range(r):
                    result[i][j] += self[i][k]*rhs[k][j]
        return result

    def __add__(self, rhs):
        """
        Matrix addition.

        >>> print(MyMatrix([[1, 2],[3, 4]]) + MyMatrix([[5, 6],[7, 8]]))
        [[6, 8]
         [10, 12]]
        """
        m = self.numRows()
        n = self.numCols()
        assert len(rhs) == m
        assert len(rhs[0]) == n
        result = zeros(m,n)
        for i in range(m):
            for j in range(n):
                result[i][j] = self[i][j] + rhs[i][j]
        return result

    def __sub__(self, rhs):
        """
        Matrix subtraction.

        >>> print(MyMatrix([[1, 2],[3, 4]]) - MyMatrix([[5, 6],[7, 8]]))
        [[-4, -4]
         [-4, -4]]
        """
        m = self.numRows()
        n = self.numCols()
        assert len(rhs) == m
        assert len(rhs[0]) == n
        result = zeros(m,n)
        for i in range(m):
            for j in range(n):
                result[i][j] = self[i][j] - rhs[i][j]
        return result

    def __pos__(self):
        return self

    def __neg__(self):
        m = self.numRows()
        n = self.numCols()
        result = zeros(m, n)
        for i in range(m):
            for j in range(n):
                result[i][j] = -self[i][j]
        return result

    def __invert__(self):
        """
        >>> m = MyMatrix([[1,1],[0,1]])
        >>> print(m)
        [[1, 1]
         [0, 1]]
        >>> print(~m)
        [[1.0, -1.0]
         [0.0, 1.0]]
        >>> print(m*~m)
        [[1.0, 0.0]
         [0.0, 1.0]]
        >>> print(~m*m)
        [[1.0, 0.0]
         [0.0, 1.0]]
        >>> m = MyMatrix([[1,0,0],[0,0,1],[0,-1,0]])
        >>> print(m)
        [[1, 0, 0]
         [0, 0, 1]
         [0, -1, 0]]
        >>> print(~m)
        [[1.0, 0.0, 0.0]
         [0.0, 0.0, -1.0]
         [0.0, 1.0, 0.0]]
        >>> print(m*~m)
        [[1.0, 0.0, 0.0]
         [0.0, 1.0, 0.0]
         [0.0, 0.0, 1.0]]
        >>> print(~m*m)
        [[1.0, 0.0, 0.0]
         [0.0, 1.0, 0.0]
         [0.0, 0.0, 1.0]]
        """
        assert self.is_square()
        if self.numRows() == 0:
            return self
        elif self.numRows() == 1:
            val = self[0][0]
            val = 1.0/val
            return MyMatrix([[val]])
        elif self.numRows() == 2: # 2x2 is analytic
            # http://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_2.C3.972_matrices
            a = self[0][0]
            b = self[0][1]
            c = self[1][0]
            d = self[1][1]
            determinant = a*d - b*c
            if determinant == 0:
                raise ArithmeticError("Cannot invert 2x2 matrix with zero determinant")
            else:
                return 1.0/(a*d - b*c) * MyMatrix([[d, -b],[-c, a]])
        else:
            # Gauss Jordan elimination from numerical recipes
            n = self.numRows()
            m1 = self.numCols()
            assert n == m1
            # Copy initial matrix into result matrix
            a = zeros(n, n)
            for i in range (0,n):
                for j in range (0,n):
                    a[i][j] = self[i][j]
            # These arrays are used for bookkeeping on the pivoting
            indxc = [0] * n
            indxr = [0] * n
            ipiv = [0] * n
            for i in range (0,n):
                big = 0.0
                for j in range (0,n):
                    if ipiv[j] != 1:
                        for k in range (0,n):
                            if ipiv[k] == 0:
                                if abs(a[j][k]) >= big:
                                    big = abs(a[j][k])
                                    irow = j
                                    icol = k
                ipiv[icol] += 1
                # We now have the pivot element, so we interchange rows...
                if irow != icol:
                    for l in range(n):
                        temp = a[irow][l]
                        a[irow][l] = a[icol][l]
                        a[icol][l] = temp
                indxr[i] = irow
                indxc[i] = icol
                if a[icol][icol] == 0:
                    raise ArithmeticError("Cannot invert singular matrix")
                pivinv = 1.0/a[icol][icol]
                a[icol][icol] = 1.0
                for l in range(n):
                    a[icol][l] *= pivinv
                for ll in range(n): # next we reduce the rows
                    if ll == icol:
                        continue # except the pivot one, of course
                    dum = a[ll][icol]
                    a[ll][icol] = 0.0
                    for l in range(n):
                        a[ll][l] -= a[icol][l]*dum
            # Unscramble the permuted columns
            for l in range(n-1, -1, -1):
                if indxr[l] == indxc[l]:
                    continue
                for k in range(n):
                    temp = a[k][indxr[l]]
                    a[k][indxr[l]] = a[k][indxc[l]]
                    a[k][indxc[l]] = temp
            return a

    def transpose(self):
        return MyMatrixTranspose(self.data)


class MyMatrixTranspose(MyMatrix):

    def transpose(self):
        return MyMatrix(self.data)

    def numRows(self):
        if len(self.data) == 0:
            return 0
        else:
            return len(self.data[0])

    def numCols(self):
        return len(self.data)

    def __getitem__(self, key):
        result = []
        for row in self.data:
            result.append(row[key])
        if isinstance(key, slice):
            return MyMatrix(result)
        else:
            return MyVector(result)

    def __setitem__(self, key, rhs):
        for n in range(len(self.data)):
            self.data[n][key] = rhs[n]

    def __str__(self):
        if len(self.data) == 0:
            return "[[]]"
        start_char = "["
        result = ""
        for m in range(len(self.data[0])):
            result += start_char
            result += "["
            sep_char = ""
            for n in range(len(self.data)):
                result += sep_char
                result += str(self.data[n][m])
                sep_char = ", "
            result += "]"
            if m < len(self.data[0]) - 1:
                result += "\n"
            start_char = " "
        result += "]"
        return result

    def __repr__(self):
        if len(self.data) == 0:
            return "MyMatrixTranspose([[]])"
        start_char = "["
        result = 'MyMatrixTranspose('
        for m in range(len(self.data[0])):
            result += start_char
            result += "["
            sep_char = ""
            for n in range(len(self.data)):
                result += sep_char
                result += repr(self.data[n][m])
                sep_char = ", "
            result += "]"
            if m < len(self.data[0]) - 1:
                result += ","
            start_char = ""
        result += '])'
        return result


# run module directly for testing
if __name__=='__main__':

    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])

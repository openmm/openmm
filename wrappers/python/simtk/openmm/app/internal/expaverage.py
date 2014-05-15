"""
expaverage.py: Exponentially Weighted Moving Average

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2013 Stanford University and the Authors.
Authors: Peter Eastman
Contributors: Robert McGibbon

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
__author__ = "Robert McGibbon"
__version__ = "1.0"

class ExponentiallyWeightedMovingAverage(object):
    """Compute an online exponentially weighted moving average
    
    """
    def __init__(self, alpha):
        """
        Parameters:
         - alpha (float) Smoothing factor in [0, 1]. Larger values of alpha
           reduce the level of smoothing. In the limiting case of alpha=1,
           the state is always equal to the last value of the input series.
        """
        self.alpha = alpha
        self.alphac = 1-alpha
        self.state = None

    def push(self, x):
        """
        Parameters:
        - x (float) Data point to add to the averager
        """
        if self.state is None:
            self.state = x
        else:
            self.state = (self.alpha * x) + (self.alphac * self.state)

    def getState(self):
        """Get the current exponentially-smoothed moving average
        """
        return self.state

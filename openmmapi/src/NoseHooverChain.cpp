/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
 * Authors: Andreas Krämer and Andrew C. Simmonett                            *
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

#include "openmm/NoseHooverChain.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

NoseHooverChain::NoseHooverChain(double temperature, double relativeTemperature, double collisionFrequency,
                                 double relativeCollisionFrequency,
                                 int numDOFs_, int chainLength_, int numMTS_,
                                 int numYoshidaSuzuki, int chainID_,
                                 const std::vector<int>& thermostatedAtoms, const std::vector<std::pair<int, int> > &thermostatedPairs):
        temp(temperature), relativeTemp(relativeTemperature), freq(collisionFrequency),
        relativeFreq(relativeCollisionFrequency), numDOFs(numDOFs_),
        chainLength(chainLength_), numMTS(numMTS_), numYS(numYoshidaSuzuki),
        chainID(chainID_), thermostatedAtoms(thermostatedAtoms), thermostatedPairs(thermostatedPairs)
{}


std::vector<double> NoseHooverChain::getYoshidaSuzukiWeights() const {
    switch (numYS) {
        case 1:
            return {1};
        case 3:
            return {0.828981543588751, -0.657963087177502, 0.828981543588751};
        case 5:
            return {0.2967324292201065, 0.2967324292201065, -0.186929716880426, 0.2967324292201065,
                    0.2967324292201065};
        default:
            throw OpenMMException("The number of Yoshida-Suzuki weights must be 1,3, or 5.");
    }
}

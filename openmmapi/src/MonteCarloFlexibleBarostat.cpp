/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Sander Vandenhaute                                 *
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

#include "openmm/MonteCarloFlexibleBarostat.h"
#include "openmm/internal/MonteCarloFlexibleBarostatImpl.h"

using namespace OpenMM;

MonteCarloFlexibleBarostat::MonteCarloFlexibleBarostat(double defaultPressure, double defaultTemperature,
                                                       int frequency, bool scaleMoleculesAsRigid) : scaleMoleculesAsRigid(scaleMoleculesAsRigid) {
    setDefaultPressure(defaultPressure);
    setDefaultTemperature(defaultTemperature);
    setFrequency(frequency);
    setRandomNumberSeed(0);
}

void MonteCarloFlexibleBarostat::setDefaultPressure(double pressure) {
    defaultPressure = pressure;
}

void MonteCarloFlexibleBarostat::setFrequency(int freq) {
    if (freq < 0)
        throw OpenMMException("Frequency cannot be negative");
    frequency = freq;
}

void MonteCarloFlexibleBarostat::setDefaultTemperature(double temp) {
    if (temp < 0)
        throw OpenMMException("Temperature cannot be negative");
    defaultTemperature = temp;
}

void MonteCarloFlexibleBarostat::computeCurrentPressure(Context& context, std::vector<double>& pressure) const {
    dynamic_cast<MonteCarloFlexibleBarostatImpl&>(getImplInContext(context)).computeCurrentPressure(getContextImpl(context), pressure);
}

ForceImpl* MonteCarloFlexibleBarostat::createImpl() const {
    return new MonteCarloFlexibleBarostatImpl(*this);
}

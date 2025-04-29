/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2025 Stanford University and the Authors.      *
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

#include "openmm/MonteCarloBarostat.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/MonteCarloBarostatImpl.h"
#include "openmm/MonteCarloAnisotropicBarostat.h"

using namespace OpenMM;

MonteCarloBarostat::MonteCarloBarostat(double defaultPressure, double defaultTemperature, int frequency) {
    setDefaultPressure(defaultPressure);
    setDefaultTemperature(defaultTemperature);
    setFrequency(frequency);
    setRandomNumberSeed(0);
}

void MonteCarloBarostat::setDefaultPressure(double pressure) {
    defaultPressure = pressure;
}

void MonteCarloBarostat::setFrequency(int freq) {
    if (freq < 0)
        throw OpenMMException("Frequency cannot be negative");
    frequency = freq;
}

void MonteCarloBarostat::setDefaultTemperature(double temp) {
    if (temp < 0)
        throw OpenMMException("Temperature cannot be negative");
    defaultTemperature = temp;
}

double MonteCarloBarostat::computeCurrentPressure(Context& context) const {
    return dynamic_cast<MonteCarloBarostatImpl&>(getImplInContext(context)).computeCurrentPressure(getContextImpl(context));
}

ForceImpl* MonteCarloBarostat::createImpl() const {
    return new MonteCarloBarostatImpl(*this);
}

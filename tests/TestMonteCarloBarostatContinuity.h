/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2026 Stanford University and the Authors.           *
 * Authors: Evan Pretti                                                       *
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

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/State.h"
#include <vector>

using namespace OpenMM;
using namespace std;

void checkContinuity(State& state, vector<Vec3>& last, vector<Vec3>& current, double threshold) {
    // Copy the previously current positions before we overwrite them.

    last.assign(current.begin(), current.end());

    // Get the current box vectors.

    Vec3 a, b, c;
    state.getPeriodicBoxVectors(a, b, c);

    // Get the current positions in reduced coordinates accounting for only the
    // box matrix diagonal, corresponding with how MonteCarloFlexibleBarostat
    // works: when the box vectors change, coordinates are rescaled with the
    // diagonal elements but not sheared with the off-diagonal elements.

    const vector<Vec3>& positions = state.getPositions();
    current.resize(positions.size());
    for(int i = 0; i < positions.size(); i++) {
        Vec3 pos = positions[i];
        current[i] = Vec3(pos[0] / a[0], pos[1] / b[1], pos[2] / c[2]);
    }

    // Make sure particles have not moved too far since the last step.  Use
    // last.size() since on the first step last will be empty and current will
    // contain the initial coordinates.

    for (int j = 0; j < last.size(); j++) {
        Vec3 delta = current[j] - last[j];
        ASSERT_USUALLY_TRUE(sqrt(delta.dot(delta)) < threshold);
    }
}

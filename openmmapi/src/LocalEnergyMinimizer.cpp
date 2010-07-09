/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

#include "openmm/LocalEnergyMinimizer.h"
#include "openmm/OpenMMException.h"
#include "lbfgs.h"
#include "openmm/Platform.h"
#include <cmath>
#include <sstream>
#include <vector>

using namespace OpenMM;
using namespace std;

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    Context* context = reinterpret_cast<Context*>(instance);
    int numParticles = context->getSystem().getNumParticles();
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++)
        positions[i] = Vec3(x[3*i], x[3*i+1], x[3*i+2]);
    context->setPositions(positions);
    State state = context->getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    for (int i = 0; i < numParticles; i++) {
        g[3*i] = -forces[i][0];
        g[3*i+1] = -forces[i][1];
        g[3*i+2] = -forces[i][2];
    }
    return state.getPotentialEnergy();
}

void LocalEnergyMinimizer::minimize(Context& context, double tolerance, int maxIterations) {
    int numParticles = context.getSystem().getNumParticles();
    lbfgsfloatval_t *x = lbfgs_malloc(numParticles*3);
    if (x == NULL)
        throw OpenMMException("LocalEnergyMinimizer: Failed to allocate memory");

    // Record the initial positions and determine a normalization constant for scaling the tolerance.
    
    const vector<Vec3>& positions = context.getState(State::Positions).getPositions();
    double norm = 0.0;
    for (int i = 0; i < numParticles; i++) {
        x[3*i] = positions[i][0];
        x[3*i+1] = positions[i][1];
        x[3*i+2] = positions[i][2];
        norm += positions[i].dot(positions[i]);
    }
    norm /= numParticles;
    norm = (norm < 1 ? 1 : sqrt(norm));

    // Perform the minimization.

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    if (!context.getPlatform().supportsDoublePrecision())
        param.xtol = 1e-7;
    param.epsilon = tolerance/norm;
    param.max_iterations = maxIterations;
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
    lbfgsfloatval_t fx;
    lbfgs(numParticles*3, x, &fx, evaluate, NULL, &context, &param);
    lbfgs_free(x);
}


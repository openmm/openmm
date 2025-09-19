#ifndef OPENMM_DRUDEHELPERS_H_
#define OPENMM_DRUDEHELPERS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
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

#include "openmm/DrudeForce.h"
#include "openmm/System.h"
#include "openmm/internal/ContextImpl.h"

namespace OpenMM {

/**
 * Find the DrudeForce contained in the System.
 */
const DrudeForce* getDrudeForce(const ContextImpl& context);

/**
 * Identify normal particles (not part of a pair) and Drude particle pairs.
 */
void findParticlesAndPairs(const ContextImpl& context, std::vector<int>& normalParticles, std::vector<std::pair<int, int> >& pairParticles);

std::vector<Vec3> assignDrudeVelocities(const ContextImpl& context, double temperature, double drudeTemperature, int randomSeed);

/**
 * Computes the instantaneous temperatures of the system and the internal Drude motion and returns a pair (T_system, T_drude)
 */
std::pair<double, double> computeTemperaturesFromVelocities(const ContextImpl& context, const std::vector<Vec3>& velocities);

}

#endif /*OPENMM_DRUDEHELPERS_H_*/

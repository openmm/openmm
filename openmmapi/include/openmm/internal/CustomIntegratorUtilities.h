#ifndef OPENMM_CUSTOMINTEGRATORUTILITIES_H_
#define OPENMM_CUSTOMINTEGRATORUTILITIES_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2015-2016 Stanford University and the Authors.      *
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

#include "openmm/CustomIntegrator.h"
#include "openmm/internal/ContextImpl.h"
#include "lepton/CustomFunction.h"
#include "lepton/ParsedExpression.h"
#include <map>
#include <string>
#include <vector>

namespace OpenMM {

class System;

/**
 * This class defines a set of utility functions that are useful in implementing CustomIntegrator.
 */

class OPENMM_EXPORT CustomIntegratorUtilities {
public:
    enum Comparison {
        EQUAL = 0, LESS_THAN = 1, GREATER_THAN = 2, NOT_EQUAL = 3, LESS_THAN_OR_EQUAL = 4, GREATER_THAN_OR_EQUAL = 5
    };
    /**
     * Parse the expression for the condition of an "if" or "while" block, and split it into
     * a left hand side, right hand side, and comparison operator.
     */
    static void parseCondition(const std::string& expression, std::string& lhs, std::string& rhs, Comparison& comparison);
    /**
     * Analyze the sequence of steps in a CustomIntegrator.  For each step:
     *
     * 1. Parse all expressions involved in the step, and identify the comparison operator (for conditional steps).
     * 2. For each conditional block, identify what step marks the end of the block.
     * 3. Determine whether the step causes previously computed forces and energies to become invalid.
     * 4. Determine whether the step itself needs forces and/or energies.
     * 5. Decide whether forces and energies should both be computed at that step (because
     *    it is more efficient to compute both at once, even if one won't be needed
     *    until a later step).
     * 6. Identify what force group each step needs forces and/or energies for.
     */
    static void analyzeComputations(const ContextImpl& context, const CustomIntegrator& integrator, std::vector<std::vector<Lepton::ParsedExpression> >& expressions,
            std::vector<Comparison>& comparisons, std::vector<int>& blockEnd, std::vector<bool>& invalidatesForces, std::vector<bool>& needsForces,
            std::vector<bool>& needsEnergy, std::vector<bool>& computeBoth, std::vector<int>& forceGroup, const std::map<std::string, Lepton::CustomFunction*>& functions);
    /**
     * Determine whether an expression involves a particular variable.
     */
    static bool usesVariable(const Lepton::ParsedExpression& expression, const std::string& variable);
private:
    static bool usesVariable(const Lepton::ExpressionTreeNode& node, const std::string& variable);
    static void enumeratePaths(int firstStep, std::vector<int> steps, std::vector<int> jumps, const std::vector<int>& blockEnd,
            const std::vector<CustomIntegrator::ComputationType>& stepType, const std::vector<bool>& needsForces, const std::vector<bool>& needsEnergy,
            const std::vector<bool>& invalidatesForces, const std::vector<int>& forceGroup, std::vector<bool>& computeBoth);
    static void analyzeForceComputationsForPath(std::vector<int>& steps, const std::vector<bool>& needsForces, const std::vector<bool>& needsEnergy,
            const std::vector<bool>& invalidatesForces, const std::vector<int>& forceGroup, std::vector<bool>& computeBoth);
    static void validateDerivatives(const Lepton::ExpressionTreeNode& node, const std::vector<std::string>& derivNames);
};

} // namespace OpenMM

#endif /*OPENMM_CUSTOMINTEGRATORUTILITIES_H_*/

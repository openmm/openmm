#ifndef LEPTON_EXPRESSION_TREE_NODE_H_
#define LEPTON_EXPRESSION_TREE_NODE_H_

/* -------------------------------------------------------------------------- *
 *                                   Lepton                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the Lepton expression parser originating from              *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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

#include "windowsIncludes.h"
#include <string>
#include <vector>

namespace Lepton {

class Operation;

class LEPTON_EXPORT ExpressionTreeNode {
public:
    ExpressionTreeNode(Operation* operation, const std::vector<ExpressionTreeNode>& children);
    ExpressionTreeNode(Operation* operation, const ExpressionTreeNode& child1, const ExpressionTreeNode& child2);
    ExpressionTreeNode(Operation* operation, const ExpressionTreeNode& child);
    ExpressionTreeNode(Operation* operation);
    ExpressionTreeNode(const ExpressionTreeNode& node);
    ExpressionTreeNode();
    ~ExpressionTreeNode();
    ExpressionTreeNode& operator=(const ExpressionTreeNode& node);
    const Operation& getOperation() const;
    const std::vector<ExpressionTreeNode>& getChildren() const;
private:
    Operation* operation;
    std::vector<ExpressionTreeNode> children;
};

} // namespace Lepton

#endif /*LEPTON_EXPRESSION_TREE_NODE_H_*/

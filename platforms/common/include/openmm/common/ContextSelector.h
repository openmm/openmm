#ifndef OPENMM_CONTEXTSELECTOR_H_
#define OPENMM_CONTEXTSELECTOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2021 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "ComputeContext.h"

namespace OpenMM {

/**
 * This class provides a safe and easy way to select a ComputeContext as current
 * for a block of code.  The constructor calls pushAsCurrent() on the context.
 * When it goes out of scope, the destructor calls popAsCurrent() on it.  Simply
 * create a local variable of this class, and the context will be current for
 * the remainder of the code block in which it is declared.  
 */

class OPENMM_EXPORT_COMMON ContextSelector {
public:
    ContextSelector(ComputeContext& context) : context(context) {
        context.pushAsCurrent();
    }
    ~ContextSelector() {
        context.popAsCurrent();
    }
private:
    ComputeContext& context;
};

} // namespace OpenMM

#endif /*OPENMM_CONTEXTSELECTOR_H_*/

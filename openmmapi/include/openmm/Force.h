#ifndef OPENMM_FORCE_H_
#define OPENMM_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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

#include "internal/windowsExport.h"

namespace OpenMM {

class Context;
class ContextImpl;
class ForceImpl;

/**
 * Force objects apply forces to the particles in a System, or alter their behavior in other
 * ways.  This is an abstract class.  Subclasses define particular forces.
 * 
 * More specifically, a Force object can do any or all of the following:
 * 
 * <ul>
 * <li>Add a contribution to the force on each particle</li>
 * <li>Add a contribution to the potential energy of the System</li>
 * <li>Modify the positions and velocities of particles at the start of each time step</li>
 * <li>Define parameters which are stored in the Context and can be modified by the user</li>
 * <li>Change the values of parameters defined by other Force objects at the start of each time step</li>
 * </ul>
 * 
 * Forces may be organized into "force groups".  This is used for multiple time step integration,
 * and allows subsets of the Forces in a System to be evaluated at different times.  By default,
 * all Forces are in group 0.  Call setForceGroup() to change this.  Some Force subclasses may
 * provide additional methods to further split their computations into multiple groups.  Be aware
 * that particular Platforms may place restrictions on the use of force groups, such as requiring
 * all nonbonded forces to be in the same group.
 */

class OPENMM_EXPORT Force {
public:
    Force() : forceGroup(0) {
    }
    virtual ~Force() {
    }
    /**
     * Get the force group this Force belongs to.
     */
    int getForceGroup() const;
    /**
     * Set the force group this Force belongs to.
     * 
     * @param group    the group index.  Legal values are between 0 and 31 (inclusive).
     */
    void setForceGroup(int group);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions. This method should be overridden for all Force subclasses, or
     * a OpenMM::OpenMMException will be thrown
     *
     * @return true if Force uses periodic boundaries or false if it does not
     */
    virtual bool usesPeriodicBoundaryConditions() const;
protected:
    friend class ContextImpl;
    /**
     * When a Context is created, it invokes this method on each Force in the System.
     * It should create a new ForceImpl object which can be used by the context for calculating forces.
     * The ForceImpl will be deleted automatically when the Context is deleted.
     */
    virtual ForceImpl* createImpl() const = 0;
    /**
     * Get the ForceImpl corresponding to this Force in a Context.
     */
    ForceImpl& getImplInContext(Context& context);
    /**
     * Get a const reference to the ForceImpl corresponding to this Force in a Context.
     */
    const ForceImpl& getImplInContext(const Context& context) const;
    /**
     * Get the ContextImpl corresponding to a Context.
     */
    ContextImpl& getContextImpl(Context& context);
private:
    int forceGroup;
};

} // namespace OpenMM

#endif /*OPENMM_FORCE_H_*/

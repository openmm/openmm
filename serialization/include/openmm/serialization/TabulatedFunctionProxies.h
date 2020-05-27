#ifndef OPENMM_TABULATEDFUNCTION_PROXIES_H_
#define OPENMM_TABULATEDFUNCTION_PROXIES_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

#include "openmm/internal/windowsExport.h"
#include "openmm/serialization/SerializationProxy.h"

namespace OpenMM {

/**
 * This is a proxy for serializing Continuous1DFunction objects.
 */
    
class OPENMM_EXPORT Continuous1DFunctionProxy : public SerializationProxy {
public:
    Continuous1DFunctionProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

/**
 * This is a proxy for serializing Continuous2DFunction objects.
 */
    
class OPENMM_EXPORT Continuous2DFunctionProxy : public SerializationProxy {
public:
    Continuous2DFunctionProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

/**
 * This is a proxy for serializing Continuous3DFunction objects.
 */
    
class OPENMM_EXPORT Continuous3DFunctionProxy : public SerializationProxy {
public:
    Continuous3DFunctionProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

/**
 * This is a proxy for serializing Discrete1DFunction objects.
 */
    
class OPENMM_EXPORT Discrete1DFunctionProxy : public SerializationProxy {
public:
    Discrete1DFunctionProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

/**
 * This is a proxy for serializing Discrete2DFunction objects.
 */
    
class OPENMM_EXPORT Discrete2DFunctionProxy : public SerializationProxy {
public:
    Discrete2DFunctionProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

/**
 * This is a proxy for serializing Discrete3DFunction objects.
 */
    
class OPENMM_EXPORT Discrete3DFunctionProxy : public SerializationProxy {
public:
    Discrete3DFunctionProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

} // namespace OpenMM

#endif /*OPENMM_TABULATEDFUNCTION_PROXIES_H_*/

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
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

#include "openmm/common/ComputeArray.h"
#include "openmm/common/ComputeContext.h"

using namespace OpenMM;

ComputeArray::ComputeArray() : impl(NULL) {
}

ComputeArray::~ComputeArray() {
    if (impl != NULL)
        delete impl;
}

ArrayInterface& ComputeArray::getArray() {
    if (impl == NULL)
        throw OpenMMException("ComputeArray has not been initialized");
    return *impl;
}

void ComputeArray::initialize(ComputeContext& context, int size, int elementSize, const std::string& name) {
    if (impl != NULL)
        throw OpenMMException("The array "+getName()+" has already been initialized");
    impl = context.createArray();
    impl->initialize(context, size, elementSize, name);
}

void ComputeArray::resize(int size) {
    if (impl == NULL)
        throw OpenMMException("ComputeArray has not been initialized");
    impl->resize(size);
}

bool ComputeArray::isInitialized() const {
    return (impl != NULL);
}

int ComputeArray::getSize() const {
    if (impl == NULL)
        throw OpenMMException("ComputeArray has not been initialized");
    return impl->getSize();
}

int ComputeArray::getElementSize() const {
    if (impl == NULL)
        throw OpenMMException("ComputeArray has not been initialized");
    return impl->getElementSize();
}

const std::string& ComputeArray::getName() const {
    if (impl == NULL)
        throw OpenMMException("ComputeArray has not been initialized");
    return impl->getName();
}

ComputeContext& ComputeArray::getContext() {
    if (impl == NULL)
        throw OpenMMException("ComputeArray has not been initialized");
    return impl->getContext();
}

void ComputeArray::upload(const void* data, bool blocking) {
    if (impl == NULL)
        throw OpenMMException("ComputeArray has not been initialized");
    impl->upload(data, blocking);
}

void ComputeArray::download(void* data, bool blocking) const {
    if (impl == NULL)
        throw OpenMMException("ComputeArray has not been initialized");
    impl->download(data, blocking);
}

void ComputeArray::copyTo(ArrayInterface& dest) const {
    if (impl == NULL)
        throw OpenMMException("ComputeArray has not been initialized");
    impl->copyTo(dest);
}
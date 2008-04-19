/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "Platform.h"
#include "PlatformException.h"
#include "Kernel.h"
#include "Stream.h"
#include "KernelFactory.h"
#include "StreamFactory.h"
#include <set>

using namespace OpenMM;
using namespace std;

std::vector<Platform*> Platform::platforms;

Platform::~Platform() {
    set<KernelFactory*> uniqueKernelFactories;
    set<StreamFactory*> uniqueStreamFactories;
    for (map<string, KernelFactory*>::const_iterator iter = kernelFactories.begin(); iter != kernelFactories.end(); ++iter)
        uniqueKernelFactories.insert(iter->second);
    for (map<string, StreamFactory*>::const_iterator iter = streamFactories.begin(); iter != streamFactories.end(); ++iter)
        uniqueStreamFactories.insert(iter->second);
    for (set<KernelFactory*>::const_iterator iter = uniqueKernelFactories.begin(); iter != uniqueKernelFactories.end(); ++iter)
        delete *iter;
    for (set<StreamFactory*>::const_iterator iter = uniqueStreamFactories.begin(); iter != uniqueStreamFactories.end(); ++iter)
        delete *iter;
}

void Platform::registerKernelFactory(std::string name, KernelFactory* factory) {
    kernelFactories[name] = factory;
}

void Platform::registerStreamFactory(std::string name, StreamFactory* factory) {
    streamFactories[name] = factory;
}

bool Platform::supportsKernels(std::vector<std::string> kernelNames) const {
    for (int i = 0; i < (int) kernelNames.size(); ++i)
        if (kernelFactories.find(kernelNames[i]) == kernelFactories.end())
            return false;
    return true;
}

Kernel Platform::createKernel(std::string name) const {
    if (kernelFactories.find(name) == kernelFactories.end())
        throw PlatformException("Called createKernel() on a Platform which does not support the requested kernel");
    return Kernel(kernelFactories.find(name)->second->createKernelImpl(name, *this));
}

Stream Platform::createStream(std::string name, int size, Stream::DataType type) const {
    if (streamFactories.find(name) == streamFactories.end())
        return Stream(getDefaultStreamFactory().createStreamImpl(name, size, type, *this));
    return Stream(streamFactories.find(name)->second->createStreamImpl(name, size, type, *this));
}

void Platform::registerPlatform(Platform* platform) {
    platforms.push_back(platform);
}

int Platform::getNumPlatforms() {
    return platforms.size();
}

Platform& Platform::getPlatform(int index) {
    return *platforms[index];
}

Platform& Platform::findPlatform(std::vector<std::string> kernelNames) {
    Platform* best = 0;
    double speed = 0.0;
    for (int i = 0; i < (int) platforms.size(); ++i) {
        if (platforms[i]->supportsKernels(kernelNames) && platforms[i]->getSpeed() > speed) {
            best = platforms[i];
            speed = best->getSpeed();
        }
    }
    if (best == 0)
        throw PlatformException("No Platform supports all the requested kernels");
    return *best;
}

#ifndef OPENMM_PLATFORM_H_
#define OPENMM_PLATFORM_H_

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

#include "Stream.h"
#include <map>
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

class Kernel;
class KernelFactory;
class OpenMMContextImpl;
class StreamFactory;

/**
 * A Platform defines an implementation of all the kernels and streams needed to perform some calculation.
 * More precisely, a Platform object acts as a registry for a set of KernelFactory and StreamFactory
 * objects which together implement the kernels and streams.  The Platform class, in turn, provides a
 * static registry of all available Platform objects.
 * 
 * To get a Platform object, call
 * 
 * <pre>
 * Platform& platform Platform::findPlatform(kernelNames);
 * </pre>
 * 
 * passing in the names of all kernels that will be required for the calculation you plan to perform.  It
 * will return the fastest available Platform which provides implementations of all the specified kernels.
 * You can then call createKernel() and createStream() to construct particular kernels and streams as needed.
 */

class OPENMM_EXPORT Platform {
public:
    virtual ~Platform();
    /**
     * Get the name of this platform.  This should be a unique identifier which can be used to recognized it.
     */
    virtual std::string getName() const = 0;
    /**
     * Get an estimate of how fast this Platform class is.  This need not be precise.  It only is expected to
     * return an order or magnitude estimate of the relative performance of different Platform classes.  An
     * unoptimized reference implementation should return 1.0, and all other Platforms should return a larger
     * value that is an estimate of how many times faster they are than the reference implementation.
     */
    virtual double getSpeed() const = 0;
    /**
     * Get whether this Platform supports double precision arithmetic.  If this returns false, the platform
     * is permitted to implement double precision streams internally as single precision.
     */
    virtual bool supportsDoublePrecision() const = 0;
    /**
     * Get the default StreamFactory for this Platform.  It will be used to create Streams whenever a
     * different StreamFactory has not been registered for the requested stream name.
     */
    virtual const StreamFactory& getDefaultStreamFactory() const = 0;
    /**
     * This is called whenever a new OpenMMContext is created.  It gives the Platform a chance to initialize
     * the context and store platform-specific data in it.
     */
    virtual void contextCreated(OpenMMContextImpl& context) const;
    /**
     * This is called whenever an OpenMMContext is deleted.  It gives the Platform a chance to clean up
     * any platform-specific data that was stored in it.
     */
    virtual void contextDestroyed(OpenMMContextImpl& context) const;
    /**
     * Register a KernelFactory which should be used to create Kernels with a particular name.
     * The Platform takes over ownership of the factory, and will delete it when the Platform itself
     * is deleted.
     * 
     * @param name     the kernel name for which the factory should be used
     * @param factory  the factory to use for creating Kernels with the specified name
     */
    void registerKernelFactory(std::string name, KernelFactory* factory);
    /**
     * Register a StreamFactory which should be used to create Streams with a particular name.
     * The Platform takes over ownership of the factory, and will delete it when the Platform itself
     * is deleted.
     * 
     * @param name     the stream name for which the factory should be used
     * @param factory  the factory to use for creating Streams with the specified name
     */
    void registerStreamFactory(std::string name, StreamFactory* factory);
    /**
     * Determine whether this Platforms provides implementations of a set of kernels.
     * 
     * @param kernelNames the names of the kernels of interests
     * @return true if this Platform provides implementations of all the kernels in the list,
     * false if there are any which it does not support
     */
    bool supportsKernels(std::vector<std::string> kernelNames) const ;
    /**
     * Create a Kernel object.  If you call this method multiple times for different contexts with the same name,
     * the returned Kernels are independent and do not interact with each other.  This means
     * that it is possible to have multiple simulations in progress at one time without them
     * interfering.
     * 
     * If no KernelFactory has been registered for the specified name, this will throw an exception.
     * 
     * @param name the name of the Kernel to get
     * @param context the context for which to create a Kernel
     * @return a newly created Kernel object
     */
    Kernel createKernel(std::string name, OpenMMContextImpl& context) const;
    /**
     * Create a Stream object.  If you call this method multiple times for different contexts with the same name,
     * the returned Streams are independent and do not interact with each other.  This means
     * that it is possible to have multiple simulations in progress at one time without them
     * interfering.
     * 
     * If a StreamFactory has been registered for the specified name, it will be used to create
     * the Stream.  Otherwise, the default StreamFactory will be used.
     * 
     * @param name the name of the Stream to get
     * @param context the context for which to create a Stream
     * @return a newly created Stream object
     */
    Stream createStream(std::string name, int size, Stream::DataType type, OpenMMContextImpl& context) const;
    /**
     * Register a new Platform.
     */
    static void registerPlatform(Platform* platform);
    /**
     * Get the number of Platforms that have been registered.
     */
    static int getNumPlatforms();
    /**
     * Get a registered Platform by index.
     */
    static Platform& getPlatform(int index);
    /**
     * Find a Platform which can be used to perform a calculation.
     * 
     * @param kernelNames the names of all kernels which will be needed for the calculation
     * @return the fastest registered Platform which supports all of the requested kernels.  If no
     * Platform exists which supports all of them, this will throw an exception.
     */
    static Platform& findPlatform(std::vector<std::string> kernelNames);
private:

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif

    std::map<std::string, KernelFactory*> kernelFactories;
    std::map<std::string, StreamFactory*> streamFactories;
    static std::vector<Platform*> platforms;

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

} // namespace OpenMM

#endif /*OPENMM_PLATFORM_H_*/

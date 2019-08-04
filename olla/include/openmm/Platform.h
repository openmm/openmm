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
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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

#include <map>
#include <string>
#include <vector>
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

class Context;
class ContextImpl;
class Kernel;
class KernelFactory;

/**
 * A Platform defines an implementation of all the kernels needed to perform some calculation.
 * More precisely, a Platform object acts as a registry for a set of KernelFactory
 * objects which together implement the kernels.  The Platform class, in turn, provides a
 * static registry of all available Platform objects.
 *
 * To get a Platform object, call
 *
 * <pre>
 * Platform& platform = Platform::findPlatform(kernelNames);
 * </pre>
 *
 * passing in the names of all kernels that will be required for the calculation you plan to perform.  It
 * will return the fastest available Platform which provides implementations of all the specified kernels.
 * You can then call createKernel() to construct particular kernels as needed.
 */

class OPENMM_EXPORT Platform {
public:
    virtual ~Platform();
    /**
     * Get the name of this platform.  This should be a unique identifier which can be used to recognized it.
     */
    virtual const std::string& getName() const = 0;
    /**
     * Get an estimate of how fast this Platform class is.  This need not be precise.  It only is expected to
     * return an order or magnitude estimate of the relative performance of different Platform classes.  An
     * unoptimized reference implementation should return 1.0, and all other Platforms should return a larger
     * value that is an estimate of how many times faster they are than the reference implementation.
     */
    virtual double getSpeed() const = 0;
    /**
     * Get whether this Platform supports double precision arithmetic.  If this returns false, the platform
     * is permitted to represent double precision values internally as single precision.
     *
     * @deprecated This method is not well defined, and is too simplistic to describe the actual behavior of
     * some Platforms, such as ones that offer multiple precision modes.  It will be removed in a future release.
     */
    virtual bool supportsDoublePrecision() const = 0;
    /**
     * Get the names of all Platform-specific properties this Platform supports.
     */
    const std::vector<std::string>& getPropertyNames() const;
    /**
     * Get the value of a Platform-specific property for a Context.
     *
     * @param context     the Context for which to get the property
     * @param property    the name of the property to get
     * @return the value of the property
     */
    virtual const std::string& getPropertyValue(const Context& context, const std::string& property) const;
    /**
     * Set the value of a Platform-specific property for a Context.
     *
     * @param context     the Context for which to set the property
     * @param property    the name of the property to set
     * @param value       the value to set for the property
     */
    virtual void setPropertyValue(Context& context, const std::string& property, const std::string& value) const;
    /**
     * Get the default value of a Platform-specific property.  This is the value that will be used for
     * newly created Contexts.
     *
     * @param property    the name of the property to get
     * @return the default value of the property
     */
    const std::string& getPropertyDefaultValue(const std::string& property) const;
    /**
     * Set the default value of a Platform-specific property.  This is the value that will be used for
     * newly created Contexts.
     *
     * @param property    the name of the property to set
     * @param value       the value to set for the property
     */
    void setPropertyDefaultValue(const std::string& property, const std::string& value);
    /**
     * This is called whenever a new Context is created.  It gives the Platform a chance to initialize
     * the context and store platform-specific data in it.
     *
     * @param context    the newly created context
     * @param properties a set of values for platform-specific properties.  Keys are the property names.
     */
    virtual void contextCreated(ContextImpl& context, const std::map<std::string, std::string>& properties) const;
    /**
     * This is called whenever a new Context is created using ContextImpl::createLinkedContext().  It gives the
     * Platform a chance to initialize the context and store platform-specific data in it.
     *
     * @param context          the newly created context
     * @param originalContext  the original context it is linked to
     */
    virtual void linkedContextCreated(ContextImpl& context, ContextImpl& originalContext) const;
    /**
     * This is called whenever a Context is deleted.  It gives the Platform a chance to clean up
     * any platform-specific data that was stored in it.
     */
    virtual void contextDestroyed(ContextImpl& context) const;
    /**
     * Register a KernelFactory which should be used to create Kernels with a particular name.
     * The Platform takes over ownership of the factory, and will delete it when the Platform itself
     * is deleted.
     *
     * @param name     the kernel name for which the factory should be used
     * @param factory  the factory to use for creating Kernels with the specified name
     */
    void registerKernelFactory(const std::string& name, KernelFactory* factory);
    /**
     * Determine whether this Platforms provides implementations of a set of kernels.
     *
     * @param kernelNames the names of the kernels of interests
     * @return true if this Platform provides implementations of all the kernels in the list,
     * false if there are any which it does not support
     */
    bool supportsKernels(const std::vector<std::string>& kernelNames) const ;
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
    Kernel createKernel(const std::string& name, ContextImpl& context) const;
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
     * Get any failures caused during the last call to loadPluginsFromDirectory
     */
    static std::vector<std::string> getPluginLoadFailures();
    /**
     * Get the registered Platform with a particular name.  If no Platform with that name has been
     * registered, this throws an exception.
     */
    static Platform& getPlatformByName(const std::string& name);
    /**
     * Find a Platform which can be used to perform a calculation.
     *
     * @param kernelNames the names of all kernels which will be needed for the calculation
     * @return the fastest registered Platform which supports all of the requested kernels.  If no
     * Platform exists which supports all of them, this will throw an exception.
     */
    static Platform& findPlatform(const std::vector<std::string>& kernelNames);
    /**
     * Load a dynamic library (DLL) which contains an OpenMM plugin.  Typically, each Platform
     * is distributed as a separate dynamic library.  This method can then be called at runtime
     * to load each available library.  Each library should contain an initializer function to
     * register any Platforms and KernelFactories that it contains.
     *
     * If the file does not exist or cannot be loaded, an exception is thrown.
     *
     * @param file   the path to the dynamic library file.  This is interpreted using the operating
     *               system's rules for loading libraries.  Typically it may be either an absolute path
     *               or relative to a set of standard locations.
     */
    static void loadPluginLibrary(const std::string& file);
    /**
     * Load multiple dynamic libraries (DLLs) which contain OpenMM plugins from one or more directories.
     * Multiple fully-qualified paths can be joined together with ':' on unix-like systems
     * (or ';' on windows-like systems); each will be searched for plugins, in-order. For example,
     * '/foo/plugins:/bar/plugins' will search both `/foo/plugins` and `/bar/plugins`. If an
     * identically-named plugin is encountered twice it will be loaded at both points; be careful!!!
     *
     * This method loops over every file contained in the specified directories and calls loadPluginLibrary()
     * for each one.  If an error occurs while trying to load a particular file, that file is simply
     * ignored. You can retrieve a list of all such errors by calling getPluginLoadFailures().
     *
     * @param directory    a ':' (unix) or ';' (windows) deliminated list of paths containing libraries to load
     * @return the names of all files which were successfully loaded as libraries
     */
    static std::vector<std::string> loadPluginsFromDirectory(const std::string& directory);
    /**
     * Get the default directory from which to load plugins.  If the environment variable
     * OPENMM_PLUGIN_DIR is set, this returns its value.  Otherwise, it returns a platform
     * specific default location. 
     *
     * @return the path to the default plugin directory
     */
    static const std::string& getDefaultPluginsDirectory();
    /**
     * Get a string containing the version number of the OpenMM library.
     */
    static const std::string& getOpenMMVersion();
protected:
    /**
     * Get the ContextImpl for a Context.
     */
    ContextImpl& getContextImpl(Context& context) const;
    /**
     * Get the ContextImpl for a Context.
     */
    const ContextImpl& getContextImpl(const Context& context) const;
    std::vector<std::string> platformProperties;
    std::map<std::string, std::string> deprecatedPropertyReplacements;
private:
    friend class ContextImpl;
    std::map<std::string, KernelFactory*> kernelFactories;
    std::map<std::string, std::string> defaultProperties;
    static std::vector<Platform*>& getPlatforms();
    static std::vector<std::string> pluginLoadFailures;
};


} // namespace OpenMM

#endif /*OPENMM_PLATFORM_H_*/

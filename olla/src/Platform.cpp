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

#include "openmm/Platform.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/Kernel.h"
#include "openmm/KernelFactory.h"
#ifdef WIN32
#include <windows.h>
#include <sstream>
#else
#include <dlfcn.h>
#include <dirent.h>
#include <cstdlib>
#endif
#include <set>

#include "ReferencePlatform.h"

using namespace OpenMM;
using namespace std;

static int registerPlatforms() {

    // Register the Platforms built into the main library.  This should eventually be moved elsewhere.
    
    ReferencePlatform* platform = new ReferencePlatform();
    Platform::registerPlatform(platform);
    return 0;
}

static int platformInitializer = registerPlatforms();

Platform::~Platform() {
    set<KernelFactory*> uniqueKernelFactories;
    for (map<string, KernelFactory*>::const_iterator iter = kernelFactories.begin(); iter != kernelFactories.end(); ++iter)
        uniqueKernelFactories.insert(iter->second);
    for (set<KernelFactory*>::const_iterator iter = uniqueKernelFactories.begin(); iter != uniqueKernelFactories.end(); ++iter)
        delete *iter;
}

const vector<string>& Platform::getPropertyNames() const {
    return platformProperties;
}

const string& Platform::getPropertyValue(const Context& context, const string& property) const {
    throw OpenMMException("getPropertyValue: Illegal property name");
}

void Platform::setPropertyValue(Context& context, const string& property, const string& value) const {
    throw OpenMMException("setPropertyValue: Illegal property name");
}

const string& Platform::getPropertyDefaultValue(const string& property) const {
    map<string, string>::const_iterator value = defaultProperties.find(property);
    if (value == defaultProperties.end())
        throw OpenMMException("getPropertyDefaultValue: Illegal property name");
    return value->second;
}

void Platform::setPropertyDefaultValue(const string& property, const string& value) {
    defaultProperties[property] = value;
}

void Platform::contextCreated(ContextImpl& context, const map<string, string>& properties) const {
}

void Platform::contextDestroyed(ContextImpl& context) const {
}

void Platform::registerKernelFactory(const string& name, KernelFactory* factory) {
    kernelFactories[name] = factory;
}

bool Platform::supportsKernels(const vector<string>& kernelNames) const {
    for (int i = 0; i < (int) kernelNames.size(); ++i)
        if (kernelFactories.find(kernelNames[i]) == kernelFactories.end())
            return false;
    return true;
}

Kernel Platform::createKernel(const string& name, ContextImpl& context) const {
    if (kernelFactories.find(name) == kernelFactories.end())
        throw OpenMMException("Called createKernel() on a Platform which does not support the requested kernel");
    return Kernel(kernelFactories.find(name)->second->createKernelImpl(name, *this, context));
}
vector<Platform*>& Platform::getPlatforms() {
    static vector<Platform*> platforms;
    return platforms;
}

void Platform::registerPlatform(Platform* platform) {
    getPlatforms().push_back(platform);
}

int Platform::getNumPlatforms() {
    return getPlatforms().size();
}

Platform& Platform::getPlatform(int index) {
    return *getPlatforms()[index];
}

Platform& Platform::getPlatformByName(const string& name) {
    for (int i = 0; i < getNumPlatforms(); i++)
        if (getPlatform(i).getName() == name)
            return getPlatform(i);
    throw OpenMMException("There is no registered Platform called \""+name+"\"");
}

Platform& Platform::findPlatform(const vector<string>& kernelNames) {
    Platform* best = 0;
    vector<Platform*>& platforms = getPlatforms();
    double speed = 0.0;
    for (int i = 0; i < (int) platforms.size(); ++i) {
        if (platforms[i]->supportsKernels(kernelNames) && platforms[i]->getSpeed() > speed) {
            best = platforms[i];
            speed = best->getSpeed();
        }
    }
    if (best == 0)
        throw OpenMMException("No Platform supports all the requested kernels");
    return *best;
}

#ifdef WIN32
static HMODULE loadOneLibrary(const string& file) {
    // Tell Windows not to bother the user with ugly error boxes.
    const UINT oldErrorMode = SetErrorMode(SEM_FAILCRITICALERRORS);
    HMODULE handle = LoadLibrary(file.c_str());
    SetErrorMode(oldErrorMode); // Restore previous error mode.
    if (handle == NULL) {
        stringstream message;
        message << "Error loading library " << file << ": " << GetLastError();
        throw OpenMMException(message.str());
    }
    return handle;
}

static void initializePlugins(vector<HMODULE>& plugins) {
    for (int i = 0; i < (int) plugins.size(); i++) {
        void (*init)();
        *(void **)(&init) = GetProcAddress(plugins[i], "registerPlatforms");
        if (init != NULL)
            (*init)();
    }
    for (int i = 0; i < (int) plugins.size(); i++) {
        void (*init)();
        *(void **)(&init) = GetProcAddress(plugins[i], "registerKernelFactories");
        if (init != NULL)
            (*init)();
    }
}
#else
static void* loadOneLibrary(const string& file) {
    void *handle = dlopen(file.c_str(), RTLD_LAZY | RTLD_GLOBAL);
    if (handle == NULL)
        throw OpenMMException("Error loading library "+file+": "+dlerror());
    return handle;
}

static void initializePlugins(vector<void*>& plugins) {
    for (int i = 0; i < (int) plugins.size(); i++) {
        void (*init)();
        *(void **)(&init) = dlsym(plugins[i], "registerPlatforms");
        if (init != NULL)
            (*init)();
    }
    for (int i = 0; i < (int) plugins.size(); i++) {
        void (*init)();
        *(void **)(&init) = dlsym(plugins[i], "registerKernelFactories");
        if (init != NULL)
            (*init)();
    }
}
#endif

void Platform::loadPluginLibrary(const string& file) {
#ifdef WIN32
    vector<HMODULE> plugins;
#else
    vector<void*> plugins;
#endif
    plugins.push_back(loadOneLibrary(file));
    initializePlugins(plugins);
}

vector<string> Platform::loadPluginsFromDirectory(const string& directory) {
    vector<string> files;
    char dirSeparator;
#ifdef WIN32
    dirSeparator = '\\';
    WIN32_FIND_DATA fileInfo;
    string filePattern(directory + dirSeparator + "*.dll");
    HANDLE findHandle = FindFirstFile(filePattern.c_str(), &fileInfo);
    if (findHandle != INVALID_HANDLE_VALUE) {
        do {
            if (fileInfo.cFileName[0] != '.')
                files.push_back(string(fileInfo.cFileName));
        } while (FindNextFile(findHandle, &fileInfo));
        FindClose(findHandle);
    }
    vector<HMODULE> plugins;
#else
    dirSeparator = '/';
    DIR* dir;
    struct dirent *entry;
    dir = opendir(directory.c_str());
    if (dir != NULL) {
        while ((entry = readdir(dir)) != NULL) {
            if (entry->d_name[0] != '.')
                files.push_back(string(entry->d_name));
        }
        closedir(dir);
    }
    vector<void*> plugins;
#endif
    vector<string> loadedLibraries;
    for (unsigned int i = 0; i < files.size(); ++i) {
        try {
            plugins.push_back(loadOneLibrary(directory+dirSeparator+files[i]));
            loadedLibraries.push_back(files[i]);
        } catch (OpenMMException& ex) {
            // Just ignore it.
        }
    }
    initializePlugins(plugins);
    return loadedLibraries;
}

const string& Platform::getDefaultPluginsDirectory() {
    char* dir = getenv("OPENMM_PLUGIN_DIR");
    static string directory;
#ifdef _MSC_VER
    if (dir != NULL)
        directory = string(dir);
    else {
        dir = getenv("PROGRAMFILES");
        if (dir == NULL)
            directory = "C:\\\\Program Files\\OpenMM\\lib\\plugins";
        else
            directory = string(dir)+"\\OpenMM\\lib\\plugins";
    }
#else
    if (dir == NULL)
        directory = "/usr/local/openmm/lib/plugins";
    else
        directory = string(dir);
#endif
    return directory;
}

// Some bizarre preprocessor magic required to convert a macro to a string...
#define STRING1(x) #x
#define STRING(x) STRING1(x)

const string& Platform::getOpenMMVersion() {
    static const string version = STRING(OPENMM_MAJOR_VERSION) "." STRING(OPENMM_MINOR_VERSION);
    return version;
}

ContextImpl& Platform::getContextImpl(Context& context) const {
    return *context.impl;
}

const ContextImpl& Platform::getContextImpl(const Context& context) const {
    return *context.impl;
}


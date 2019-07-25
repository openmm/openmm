#ifndef OPENMM_PLUGININITIALIZER_H_
#define OPENMM_PLUGININITIALIZER_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2011 Stanford University and the Authors.      *
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

/**
 * This file contains declarations for initialization functions that must be defined
 * by plugins, and that are invoked after the plugins have been loaded.  There are
 * two such functions: one for registering new Platforms, and one for adding new
 * KernelFactories to existing Platforms.
 */

/**
 * This function registers new Platforms that are defined by the plugin.  It will be
 * invoked after the plugin is loaded, and should register the new Platforms by
 * calling Platform::registerPlatform().
 */
extern "C" void registerPlatforms();

/**
 * This function registers new KernelFactories for existing Platforms.
 * It will be invoked after the plugin is loaded, and should register
 * the new factories by calling registerKernelFactory() on the appropriate Platform objects.
 * It is not invoked until after registerPlatforms() has been called on every plugin,
 * thus avoiding initialization order problems when one plugin adds a KernelFactory
 * to a Platform defined by another plugin.
 */
extern "C" void registerKernelFactories();

#endif /*OPENMM_PLUGININITIALIZER_H_*/

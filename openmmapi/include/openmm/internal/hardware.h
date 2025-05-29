#ifndef OPENMM_HARDWARE_H_
#define OPENMM_HARDWARE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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
 * This file defines a collection of functions for querying the specific hardware being used.
 */

#ifdef __APPLE__
   #include <sys/sysctl.h>
   #include <dlfcn.h>
#else
   #ifdef WIN32
      #define NOMINMAX
      #include <windows.h>
      #include <intrin.h>
   #else
      #ifdef __ANDROID__
        #include <cpu-features.h>
      #else
        #include <unistd.h>
      #endif
   #endif
#endif

/**
 * Get the number of CPU cores available.
 */
static int getNumProcessors() {
#ifdef __APPLE__
    int ncpu;
    size_t len = 4;
    if (sysctlbyname("hw.logicalcpu", &ncpu, &len, NULL, 0) == 0)
       return ncpu;
    else
       return 1;
#else
#ifdef WIN32
    SYSTEM_INFO siSysInfo;
    int ncpu;
    GetSystemInfo(&siSysInfo);
    ncpu = siSysInfo.dwNumberOfProcessors;
    if (ncpu < 1)
        ncpu = 1;
    return ncpu;
#else
    #ifdef __ANDROID__
        return android_getCpuCount();
    #else
      long nProcessorsOnline = sysconf(_SC_NPROCESSORS_ONLN);
      if (nProcessorsOnline == -1)
          return 1;
      else
          return (int) nProcessorsOnline;
    #endif
#endif
#endif
}

/**
 * Get a description of the CPU's capabilities.
 */
#ifdef WIN32
#define cpuid __cpuid
#else
#if !defined(__ANDROID__) && !defined(__PNACL__) && !defined(__PPC__) \
    && !defined(__ARM__) && !defined(__ARM64__) && !defined(__LOONGARCH64__)
    static void cpuid(int cpuInfo[4], int infoType) {
    #ifdef __LP64__
        __asm__ __volatile__ (
            "cpuid":
            "=a" (cpuInfo[0]),
            "=b" (cpuInfo[1]),
            "=c" (cpuInfo[2]),
            "=d" (cpuInfo[3]) :
            "a" (infoType)
        );
    #else
        __asm__ __volatile__ (
            "pushl %%ebx\n"
            "cpuid\n"
            "movl %%ebx, %1\n"
            "popl %%ebx\n" :
            "=a" (cpuInfo[0]),
            "=r" (cpuInfo[1]),
            "=c" (cpuInfo[2]),
            "=d" (cpuInfo[3]) :
            "a" (infoType)
        );
    #endif
    }
#else
    static void cpuid(int cpuInfo[4], int infoType) {
        cpuInfo[0] = cpuInfo[1] = cpuInfo[2] = 0;
    }
#endif
#endif

/**
 * Get whether this is an x86 CPU that supports AVX.
 */
static bool isAvxSupported() {
    int cpuInfo[4];
    cpuid(cpuInfo, 0);
    if (cpuInfo[0] >= 1) {
        cpuid(cpuInfo, 1);
        return ((cpuInfo[2] & ((int) 1 << 28)) != 0);
    }
    return false;
}

/**
 * Get the maximum supported size for vectors in multiples of four bytes.  This
 * is the number of int or float values that can be contained in a vector.
 */
static int getVectorWidth() {
    if (isAvxSupported())
        return 8;
    return 4;
}
#endif // OPENMM_HARDWARE_H_

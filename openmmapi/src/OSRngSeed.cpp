/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
 * Authors: Robert T. McGibbon                                                *
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

#if defined(_WIN32) || defined(__CYGWIN__)
#include <windows.h>
static HCRYPTPROV hCryptProv = 0;
#pragma comment(lib, "advapi32.lib")
#else
#include <fcntl.h>
#include <unistd.h>
#endif
#include "openmm/OpenMMException.h"
#include "openmm/internal/OSRngSeed.h"

using OpenMM::OpenMMException;

int osrngseed(void) {
    int value;
#if defined(_WIN32) || defined(__CYGWIN__)
    if (!::CryptAcquireContextW(&hCryptProv, 0, 0, PROV_RSA_FULL, CRYPT_VERIFYCONTEXT | CRYPT_SILENT)) {
        throw OpenMMException("Failed to initialize Windows random API (CryptoGen)");
    }
    if (!CryptGenRandom(hCryptProv, sizeof(int), (BYTE*) &value)) {
        ::CryptReleaseContext(hCryptProv, 0);
        throw OpenMMException("Failed to get random numbers");
    }
    if (!::CryptReleaseContext(hCryptProv, 0)) {
        throw OpenMMException("Failed to release Windows random API context");
    }
#else
    int m_fd = open("/dev/urandom", O_RDONLY);
    if (m_fd == -1) {
        throw OpenMMException("Failed to open /dev/urandom");
    }
    if (read(m_fd, &value, sizeof(int)) != sizeof(int)) {
        throw OpenMMException("Failed to read bytes from /dev/urandom");
    }
    close(m_fd);
#endif
    return value;
}
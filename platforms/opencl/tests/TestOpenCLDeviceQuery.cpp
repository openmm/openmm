/**
 * This file is adapted from vexcl's (https://github.com/ddemidov/vexcl)
 * example "devlist.cpp", which is
 *
 * Copyright (c) 2012-2014 Denis Demidov <dennis.demidov@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <set>
#include <algorithm>
#include <vector>

#include "OpenCLContext.h"

using namespace std;

#define SHOW_DEVPROP(name)                                      \
    cout << "    " << left << setw(32) << #name << " = "        \
    << d.getInfo< name >() << endl


int main() {
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    cout << "OpenCL devices:" << endl << endl;

    for (int j = 0; j < platforms.size(); j++) {
        vector<cl::Device> devices;
        try {
            platforms[j].getDevices(CL_DEVICE_TYPE_ALL, &devices);
        }
        catch (...) {
            // There are no devices available for this platform.
            continue;
        }

        for (int i = 0; i < devices.size(); i++) {
            cl::Device d = devices[i];
            cout << "OpenCLPlatformIndex " << j << ", OpenCLDeviceIndex " << i << ": \"" << d.getInfo<CL_DEVICE_NAME>()
                 << "\"" << endl << "    " << left << setw(32) << "CL_PLATFORM_NAME" << " = "
                 << cl::Platform(d.getInfo<CL_DEVICE_PLATFORM>()).getInfo<CL_PLATFORM_NAME>()
                 << endl
                 << "    " << left << setw(32) << "CL_PLATFORM_VENDOR" << " = "
                 << platforms[j].getInfo<CL_PLATFORM_VENDOR>()
                 << endl;

            SHOW_DEVPROP(CL_DEVICE_VENDOR);
            SHOW_DEVPROP(CL_DEVICE_VERSION);
            cout << "    " << left << setw(32) << "CL_DEVICE_TYPE" << " = ";
            if (d.getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU) {
                cout << "CL_DEVICE_TYPE_CPU" << endl;
            } else if (d.getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_GPU) {
                cout << "CL_DEVICE_TYPE_GPU" << endl;
            } else if (d.getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_ACCELERATOR) {
                cout << "CL_DEVICE_TYPE_ACCELERATOR" << endl;
            } else {
                cout << "Unknown" << endl;
            }


            SHOW_DEVPROP(CL_DEVICE_MAX_COMPUTE_UNITS);
            cout << "    " << left << setw(32) << "CL_DEVICE_MAX_WORK_ITEM_SIZES" << " = [";
            for (int k = 0; k < d.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>().size(); k++) {
                cout << d.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()[k];
                if (k < d.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>().size() - 1)
                    cout << ", ";
            }
            cout << "]" << endl;

            SHOW_DEVPROP(CL_DEVICE_HOST_UNIFIED_MEMORY);
            SHOW_DEVPROP(CL_DEVICE_GLOBAL_MEM_SIZE);
            SHOW_DEVPROP(CL_DEVICE_LOCAL_MEM_SIZE);
            SHOW_DEVPROP(CL_DEVICE_MAX_MEM_ALLOC_SIZE);
            SHOW_DEVPROP(CL_DEVICE_ADDRESS_BITS);
            SHOW_DEVPROP(CL_DEVICE_MAX_CLOCK_FREQUENCY);

            int processingElementsPerComputeUnit;
            if (d.getInfo<CL_DEVICE_TYPE>() != CL_DEVICE_TYPE_GPU) {
                processingElementsPerComputeUnit = 1;
            } else if  (d.getInfo<CL_DEVICE_EXTENSIONS>().find("cl_nv_device_attribute_query") != string::npos) {
                cl_uint computeCapabilityMajor;
#ifdef CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV
                clGetDeviceInfo(d(), CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, sizeof(cl_uint), &computeCapabilityMajor, NULL);
                processingElementsPerComputeUnit = (computeCapabilityMajor < 2 ? 8 : 32);
#endif
            } else if (d.getInfo<CL_DEVICE_EXTENSIONS>().find("cl_amd_device_attribute_query") != string::npos) {
#ifdef CL_DEVICE_SIMD_PER_COMPUTE_UNIT_AMD
                try {
                    processingElementsPerComputeUnit = d.getInfo<CL_DEVICE_SIMD_PER_COMPUTE_UNIT_AMD>() *
                        d.getInfo<CL_DEVICE_SIMD_WIDTH_AMD>() *
                        d.getInfo<CL_DEVICE_SIMD_INSTRUCTION_WIDTH_AMD>();
                } catch (cl::Error err) {}
#endif
            }
            cout << "    processingElementsPerComputeUnit" << " = " << processingElementsPerComputeUnit << endl;
            int speed = devices[i].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>()*processingElementsPerComputeUnit*d.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>();
            cout << "    estimatedSpeed                   = " << speed << endl;


            cout << "    " << left << setw(32) << "CL_DEVICE_EXTENSIONS" << " = ";
            istringstream iss(d.getInfo<CL_DEVICE_EXTENSIONS>());
            set<string> extensions;

            extensions.insert(istream_iterator<string>(iss), istream_iterator<string>());
            size_t w = 40;
            for (set<string>::iterator s = extensions.begin(); s != extensions.end(); ++s) {
                w += s->length() + 1;
                if (w > 80) {
                    cout << endl << setw(w = 8) << "";
                    w += s->length() + 1;
                }
                cout << *s << " ";
            }
            cout << endl << endl;
        }
    }
}

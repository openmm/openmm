/**
 * This file is adapted from vexcl's (https://github.com/ddemidov/vexcl)
 * example "devlist.cpp", which is
 *
 *  Copyright (c) 2012-2014 Denis Demidov <dennis.demidov@gmail.com>
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

#define SHOW_DEVPROP(name)                               \
    cout << "    " << left << setw(32) << #name << " = " \
         << d.getInfo< name >() << endl


int main() {
  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);
  cout << "OpenCL devices:" << endl << endl;

  for (int j = 0; j < platforms.size(); j++) {
    vector<cl::Device> devices;
    platforms[j].getDevices(CL_DEVICE_TYPE_ALL, &devices);
    
    for (int i = 0; i < devices.size(); i++) {
      cl::Device d = devices[i];
      cout << "OpenCL Platform " << j << ", Device " << i << ": \"" << d.getInfo<CL_DEVICE_NAME>()
           << "\"" << endl << "    " << left << setw(32) << "CL_PLATFORM_NAME" << " = "
           << cl::Platform(d.getInfo<CL_DEVICE_PLATFORM>()).getInfo<CL_PLATFORM_NAME>()
           << endl
           << "    " << left << setw(32) << "CL_PLATFORM_VENDOR" << " = "
           << platforms[j].getInfo<CL_PLATFORM_VENDOR>()
           << endl;

      SHOW_DEVPROP(CL_DEVICE_VENDOR);
      SHOW_DEVPROP(CL_DEVICE_VERSION);
      SHOW_DEVPROP(CL_DEVICE_MAX_COMPUTE_UNITS);
      SHOW_DEVPROP(CL_DEVICE_HOST_UNIFIED_MEMORY);
      SHOW_DEVPROP(CL_DEVICE_GLOBAL_MEM_SIZE);
      SHOW_DEVPROP(CL_DEVICE_LOCAL_MEM_SIZE);
      SHOW_DEVPROP(CL_DEVICE_MAX_MEM_ALLOC_SIZE);
      SHOW_DEVPROP(CL_DEVICE_ADDRESS_BITS);
      SHOW_DEVPROP(CL_DEVICE_MAX_CLOCK_FREQUENCY);
      cout << "    " << left << setw(32) << "CL_DEVICE_EXTENSIONS" << " = ";
      {
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
      }
      cout << endl << endl;
    }
  }
}

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
           << "    " << left << setw(32) << "CL_PLATFORM_NAME" << " = "
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

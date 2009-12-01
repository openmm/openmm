#include "OpenMM.h"
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>

using namespace OpenMM;
using namespace std;

void mysetenv(const char* name, const char* value) {
#ifdef _MSC_VER
    char buffer[1000];
    sprintf(buffer, "%s=%s", name, value);
    putenv(buffer);
#else
    setenv(name, value, 1);
#endif
}

void testLoadNMLPlugin() 
{
    string pluginDir = Platform::getDefaultPluginsDirectory();
    cout << "Default plugins directory = " << pluginDir << endl;
    Platform::loadPluginsFromDirectory(pluginDir);
    vector<string> kernelName;
    kernelName.push_back("IntegrateNMLStepKernel");
    // Was NormalModeLangevin plugin loaded?
    Platform& platform = Platform::findPlatform(kernelName); // throws if no platform with kernel
}

int main(int argc, const char* argv[]) 
{
    try 
    {
        if (argc > 1) {
            const char* plugin_dir = argv[1];
            // 0 => don't set if variable exists
            mysetenv("OPENMM_PLUGIN_DIR", plugin_dir);
            cout << plugin_dir << endl;
            cout << getenv("OPENMM_PLUGIN_DIR") << endl;
        }
        testLoadNMLPlugin();
        cout << "tests passed" << endl;
        return 0;
    } 
    catch (...) 
    {
        cout << "FAILED" << endl;
        return 1;
    }
}

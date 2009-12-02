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

    // Create a context, just to initialize all plugins
    System system;
    VerletIntegrator integrator(0.01);
    Context context(system, integrator);

    // Was NormalModeLangevin plugin loaded?
    vector<string> kernelName;
    kernelName.push_back("IntegrateNMLStep");
    cout << "Searching for kernel IntegrateNMLStep" << endl;
    Platform& platform = Platform::findPlatform(kernelName); // throws if no platform with kernel
}

int main(int argc, const char* argv[]) 
{
    try 
    {
        // Set OPENMM_PLUGIN_DIR from first command line argument
        if (argc > 1) {
            const char* plugin_dir = argv[1];
            mysetenv("OPENMM_PLUGIN_DIR", plugin_dir);
            cout << plugin_dir << endl;
            cout << getenv("OPENMM_PLUGIN_DIR") << endl;
        }
        testLoadNMLPlugin();
        cout << "tests passed" << endl;
        return 0;
    } 
    catch (const std::exception& exc) 
    {
        cout << "FAILED: " << exc.what() << endl;
        return 1;
    }
}


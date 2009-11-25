#include "OpenMM.h"
#include <vector>
#include <string>
#include <iostream>

using namespace OpenMM;
using namespace std;

void testLoadNMLPlugin() 
{
    Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
    vector<string> kernelName;
    kernelName.push_back("IntegrateNMLStepKernel");
    // Was NormalModeLangevin plugin loaded?
    Platform& platform = Platform::findPlatform(kernelName); // throws if no platform with kernel
}

int main() 
{
    try 
    {
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

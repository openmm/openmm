#include "OpenMM.h"
#include <vector>
#include <string>

using namespace OpenMM;
using namespace std;

void testLoadNMLPlugin() 
{
    Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
    vector<string> kernelName;
    kernelName.push_back("IntegrateNMLStepKernel");
    // Was NormalModeLangevin plugin loaded?
    const Platform& platform = Platform::findPlatform(kernelName); // throws if no platform with kernel
}

int main() 
{
    testLoadNMLPlugin();
    return 0;
}

#include "gpu_usage.h"
#ifdef __APPLE__

#include <CoreFoundation/CoreFoundation.h>
#include <Cocoa/Cocoa.h>
#include <IOKit/IOKitLib.h>
#include <IOKit/graphics/IOGraphicsLib.h>

float getGPUCoreUsage(int gpu_number) {

    static ssize_t gpuUsage=0;
    int i = 0;
    CFMutableDictionaryRef dict = IOServiceMatching(kIOAcceleratorClassName);
    io_iterator_t it;

    if (IOServiceGetMatchingServices(kIOMasterPortDefault,dict,
                                     &it) == kIOReturnSuccess) {
        io_registry_entry_t entry;
        while ((entry = IOIteratorNext(it))) {
            CFMutableDictionaryRef serviceDictionary;
            if (IORegistryEntryCreateCFProperties(entry,
                                                  &serviceDictionary,
                                                  kCFAllocatorDefault,
                                                  kNilOptions) != kIOReturnSuccess) {
                IOObjectRelease(entry);
                IOObjectRelease(it);
                return 0;
            }

            CFMutableDictionaryRef properties = (CFMutableDictionaryRef) CFDictionaryGetValue( serviceDictionary, CFSTR("PerformanceStatistics") );
            if (properties) {
                const void* gpuCoreUtilization = CFDictionaryGetValue(properties, CFSTR("GPU Core Utilization"));
                if (gpuCoreUtilization) {
                    CFNumberGetValue( (CFNumberRef) gpuCoreUtilization, kCFNumberSInt64Type, &gpuUsage);
                }
            }

            CFRelease(serviceDictionary);
            IOObjectRelease(entry);
            if (i == gpu_number) {
               IOObjectRelease(it);
               return gpuUsage/10000000;
            }
            i++;
        }
        IOObjectRelease(it);
    }
    return 0;
}
#endif


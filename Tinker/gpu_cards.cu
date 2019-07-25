#include "gpu_usage.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

#ifdef __linux__
#include <nvml.h>
#endif

#define GIGABYTE 1073741824.0f

#define cudaSuccess_SAFE_CALL(function)    \
do {                                       \
   if (function != cudaSuccess) {          \
      bestDevice = -1;                     \
      return -1;                           \
   }                                       \
} while (0)

#define CUDA_SUCCESS_SAFE_CALL(function)   \
do {                                       \
   if (function != CUDA_SUCCESS) {         \
      bestDevice = -1;                     \
      return -1;                           \
   }                                       \
} while (0)

static int bestDevice = -1;

int findBestCUDACard() {
   int nDevices;
   const int verbose = 1;
   
   // We only want to check for the best device once. Running
   // analyze first before the simulation causes misreporting
   // of resources and leads to the wrong gpu being selected

   if (bestDevice != -1)
      return bestDevice;

   cudaSuccess_SAFE_CALL(cudaGetDeviceCount(&nDevices));

   float device_loads[nDevices];
   float device_gflops[nDevices];
   int device_used_mem[nDevices];


   printf ("\n Number of CUDA Dards Detected :    %d\n", nDevices);

   // Loop through all the CUDA devices and pick the best one

   for (int i = 0; i < nDevices; ++i) {
      cudaDeviceProp prop;
      CUdevice device;
      CUcontext context;
      size_t freeMemory;
      size_t totalMemory;
      int usedMemory;
      double GFLOPS;
      int coresPerMP = 0, MPCount, IPC = 1;
      int major, minor;

      cudaSuccess_SAFE_CALL (cudaGetDeviceProperties(&prop, i));
      CUDA_SUCCESS_SAFE_CALL (cuDeviceGet(&device, i));
      CUDA_SUCCESS_SAFE_CALL (cuCtxCreate(&context, 0, device));
      CUDA_SUCCESS_SAFE_CALL (cuMemGetInfo(&freeMemory, &totalMemory));

      usedMemory = totalMemory - freeMemory;
      major = prop.major;
      minor = prop.minor;
      MPCount = prop.multiProcessorCount;

      // This will need to be updated when new cards are released
      // Data from https://en.wikipedia.org/wiki/CUDA from the
      // Architecture specifications chart. IPC is equal to the
      // "Number of instructions issued at once by scheduler"

      switch (major) {
         case 1:
            coresPerMP = 8;
         break;
         case 2:
            if (minor == 0)
               coresPerMP = 32;
            else if (minor == 1) {
               coresPerMP = 48;
               IPC = 2;
            }
         break;
         case 3:
            coresPerMP = 192;
            IPC = 2;
         break;
         case 5:
            coresPerMP = 128;
            IPC = 2;
         break;
         case 6:
            IPC = 2;
            if (minor == 0)
               coresPerMP = 64;
            else
               coresPerMP = 128;
         break;
      }


      // GFLOPS = CUDA Cores * Clockspeed * Instructions per Clock

      GFLOPS = MPCount * coresPerMP * (double)prop.clockRate / (1000.0 * 1000.0) * IPC;

      // If the data for coresPerMP was not present, then this
      // is a newer GPU and we will assume it is fast
 
      if (coresPerMP == 0) {
         GFLOPS = INT_MAX;
      }
      float gpu_load = 0.0f;
#ifdef __APPLE__
      gpu_load = getGPUCoreUsage(i);
#endif
#ifdef __linux__
      if (NVML_SUCCESS != nvmlInit())
          printf("Failure to initialize NVML\n");
      else {
          nvmlUtilization_t gpuUtil;
          nvmlDevice_t nvmlDevice;
          unsigned int nvmlClock;
          if (NVML_SUCCESS != nvmlDeviceGetHandleByIndex(nDevices-i-1, &nvmlDevice))
              printf("Failure to find NVML device\n");
          else {
              nvmlDeviceGetUtilizationRates(nvmlDevice, &gpuUtil);
              gpu_load = (float)gpuUtil.gpu;
              nvmlDeviceGetClockInfo(nvmlDevice, NVML_CLOCK_SM, &nvmlClock);
          }
      }
      //nvmlDeviceGetUtilizationRates(nvmlDevice, &gpuUtil);
#endif

#ifdef __WIN32
      // TODO; Windows version not available at this time
#endif

      device_loads[i] = gpu_load;
      device_gflops[i] = GFLOPS;
      device_used_mem[i] = usedMemory;

      // Total used memory is used to guess if the device
      // is already running a job

      if (verbose) {
         printf("\n Device Number :  %d\n", i);
         printf("\tDevice Name              %s\n", prop.name);
         if (coresPerMP != 0) {
            printf("\tCUDA Cores               %d\n", MPCount * coresPerMP);
            printf("\tTFLOPS                   %.2f\n", GFLOPS/1000);
         }
         printf("\tClockspeed (GHz)         %.3f\n", ((float)prop.clockRate) / 1000000.0f);
         printf("\tTotal Memory (GB)        %.2f\n", round(10.0f * ((float)totalMemory/GIGABYTE) / 10.0f));
         printf("\tFree Memory (GB)         %.2f\n", round(1000.0f * (float)freeMemory/GIGABYTE) / 1024.0f);
         printf("\tGPU load                 %.2f%%\n", gpu_load);
      }
   }

   if (nDevices == 1) {
      bestDevice = 0;
      return bestDevice;
   }

   // Figure out which GPU we want to use

   float lowestLoad = 101;
   bool nonzero_load = false;

   // Pick the device that has the lowest load

   for (int i = 0; i < nDevices; ++i) {
      if (device_loads[i] > 0)
         nonzero_load = true;
      if (device_loads[i] < lowestLoad) {
         lowestLoad = device_loads[i];
         bestDevice = i;
      }
   }
   if (nonzero_load)
      return bestDevice;
   bestDevice = 0;

   // If none of the devices are in use, find the fastest device

   float highestGFLOPS = 0;
   bool highest_performance[nDevices];
   for (int i = 0; i < nDevices; ++i) {
      highest_performance[i] = false;
      if (device_gflops[i] >= highestGFLOPS) {
         highest_performance[i] = true;
         if (device_gflops[i] == highestGFLOPS) {
            for (int j = 0; j < i; ++j)
               highest_performance[j] = false;
         }
      }
   }
   
   // From the list of fastest devices, find the one with
   // the least amount of used memory

   int lowestMem = device_used_mem[0];
   bestDevice = 0;
   for (int i = 1; i < nDevices; ++i) {
      if (highest_performance[i] && device_used_mem[i] < lowestMem) {
         lowestMem = device_used_mem[i];
         bestDevice = i;
      }
   }
   return bestDevice;
}

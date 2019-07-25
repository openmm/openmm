/* 
 * file has_opencl_gpu.c 
 * created May 6, 2010
 * by Christopher Bruns
 * Returns zero if an OpenCL-capable GPU is found.
 * Returns one otherwise.
 */

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <stdio.h>

/* 
 * check_devices() looks for a GPU among all OpenCL devices 
 * in a particular OpenCL platform.
 * Returns zero if a GPU is found.  Returns one otherwise.
 */
int check_devices(cl_platform_id platform_id) 
{
    cl_int err;
    cl_device_id devices[10];
    cl_uint num_devices;
    size_t d;
    char dname[500];
    size_t namesize;
    
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 10, devices, &num_devices);
    if (err != CL_SUCCESS) {
        printf("clGetDeviceIDs error: %d", err);
        switch(err) {
            case CL_INVALID_PLATFORM :
                printf(" INVALID_PLATFORM; platform is not a valid platform.\n");
                break;
            case CL_INVALID_DEVICE_TYPE :
                printf(" CL_INVALID_DEVICE_TYPE; device_type is not a valid value.\n");
                break;
            case CL_INVALID_VALUE :
                printf(" CL_INVALID_VALUE; num_entries is equal to zero and device_type is not NULL or if both num_devices and device_type  are NULL.\n");
                break;
            case CL_DEVICE_NOT_FOUND :
                printf(" CL_DEVICE_NOT_FOUND; no OpenCL devices that matched device_type were found.\n");
                break;
            default :
                printf(" unrecognized error code '%d'\n", err);
                break;
        }
        return 1;
    }
    printf("%d devices found\n", num_devices);
    if (num_devices < 1)
        return 1;
    for (d = 0; d < num_devices; ++d) {
        clGetDeviceInfo(devices[d], CL_DEVICE_NAME, 500, dname, &namesize);
        printf("Device #%d name = %s\n", d, dname);
    }
    return 0;
}

int main(int argc, char** argv) 
{
    size_t p;
    cl_int err;
    cl_platform_id platforms[10];
    cl_uint num_platforms;
    int status;

    /* Investigate all valid platform IDs */
    err = clGetPlatformIDs(10, platforms, &num_platforms);
    if (num_platforms < 1) {
        printf("No OpenCL platforms found.\n");
        return 1;
    }
    for (p = 0; p < num_platforms; ++p) {
        printf("Checking platform ID %d\n", p);
        status = check_devices(platforms[p]);
        if (status == 0)
            return status; // found GPU
    }
    return 1; // did NOT find GPU
}

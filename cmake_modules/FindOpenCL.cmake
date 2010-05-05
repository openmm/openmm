
### OPENCL_INCLUDE_DIR ###
# Try OPENCL_DIR variable before looking elsewhere
find_path(OPENCL_INCLUDE_DIR 
    NAMES OpenCL/cl.h CL/cl.h 
    PATHS "$ENV{OPENCL_DIR}"
    PATH_SUFFIXES "include"
    NO_DEFAULT_PATH
)
# As a last resort, look in default include areas and elsewhere
find_path(OPENCL_INCLUDE_DIR 
    NAMES OpenCL/cl.h CL/cl.h
    PATHS
        "$ENV{CUDA_INC_PATH}"
        "C:/CUDA"
    PATH_SUFFIXES "include"
)

### OPENCL_LIBRARY ###
if("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
    set(path_suffixes "lib/x86")
    if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
      set(path_suffixes "lib/x86_64")
    endif("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
endif("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
# First look in OPENCL_DIR variable location
find_library(OPENCL_LIBRARY
    NAMES OpenCL
    PATHS $ENV{OPENCL_DIR} ${OPENCL_LIB_SEARCH_PATH} ""
    PATH_SUFFIXES ${path_suffixes} "lib"
    NO_DEFAULT_PATH
)
# If above fails, look in system default and other locations
find_library(OPENCL_LIBRARY
    NAMES OpenCL
    PATHS
        "$ENV{CUDA_LIB_PATH}"
        "C:/CUDA"
    PATH_SUFFIXES ${path_suffixes} "lib"
)

find_package_handle_standard_args(OPENCL DEFAULT_MSG OPENCL_LIBRARY OPENCL_INCLUDE_DIR)

if(OPENCL_FOUND)
  set(OPENCL_LIBRARIES ${OPENCL_LIBRARY})
  mark_as_advanced(CLEAR OPENCL_INCLUDE_DIR)
  mark_as_advanced(CLEAR OPENCL_LIBRARY)
else(OPENCL_FOUND)
  set(OPENCL_LIBRARIES)
  mark_as_advanced(OPENCL_INCLUDE_DIR)
  mark_as_advanced(OPENCL_LIBRARY)
endif(OPENCL_FOUND)

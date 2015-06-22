
### OPENCL_INCLUDE_DIR ###
# Try OPENCL_DIR variable before looking elsewhere
find_path(OPENCL_INCLUDE_DIR
    NAMES OpenCL/opencl.h CL/opencl.h
    PATHS $ENV{OPENCL_DIR}
    PATH_SUFFIXES "include"
    NO_DEFAULT_PATH
)
# Next look in environment variables set by OpenCL SDK installations
find_path(OPENCL_INCLUDE_DIR
    NAMES OpenCL/opencl.h CL/opencl.h
    PATHS
        $ENV{CUDA_PATH}
        $ENV{AMDAPPSDKROOT}
    PATH_SUFFIXES "include"
    NO_DEFAULT_PATH
)
# On Macs, look inside the platform SDK
if(DEFINED CMAKE_OSX_SYSROOT)
    find_path(OPENCL_INCLUDE_DIR
        NAMES opencl.h opencl.h
        PATHS
            "${CMAKE_OSX_SYSROOT}/System/Library/Frameworks/OpenCL.framework/Headers"
        NO_DEFAULT_PATH
    )
endif(DEFINED CMAKE_OSX_SYSROOT)
# As a last resort, look in default system areas followed by other possible locations
find_path(OPENCL_INCLUDE_DIR
    NAMES OpenCL/opencl.h CL/opencl.h
    PATHS
        "C:/CUDA"
        "/usr/local/cuda"
        "/usr/local/streamsdk"
        "/usr"
        "${CUDA_TOOLKIT_ROOT_DIR}"
    PATH_SUFFIXES "include"
)

### OPENCL_LIBRARY ###
if("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
    if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
      set(path_suffixes "lib/x86_64")
    else("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
      set(path_suffixes "lib/x86")
    endif("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
elseif(MSVC)
    if(CMAKE_CL_64)
      set(path_suffixes "lib/x64" "lib/x86_64")
    else(CMAKE_CL_64)
      set(path_suffixes "lib/Win32" "lib/x86")
    endif(CMAKE_CL_64)
else(MSVC)
    set(path_suffixes "lib")
endif("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
# Try OPENCL_DIR variable before looking elsewhere
find_library(OPENCL_LIBRARY
    NAMES OpenCL
    PATHS
      $ENV{OPENCL_DIR}
      ${OPENCL_LIB_SEARCH_PATH}
    PATH_SUFFIXES ${path_suffixes}
    NO_DEFAULT_PATH
)
# Next look in environment variables set by OpenCL SDK installations
find_library(OPENCL_LIBRARY
    NAMES OpenCL
    PATHS
      $ENV{CUDA_PATH}
      $ENV{AMDAPPSDKROOT}
    PATH_SUFFIXES ${path_suffixes}
    NO_DEFAULT_PATH
)
# As a last resort, look in default system areas followed by other possible locations
find_library(OPENCL_LIBRARY
    NAMES OpenCL
    PATHS
        "C:/CUDA"
        "/usr/local/cuda"
        "/usr/local/streamsdk"
        "/usr"
        "${CUDA_TOOLKIT_ROOT_DIR}"
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

# Given the root directory of a CUDA toolkit, find the version number and paths
# to headers and libraries.

function(parse_cuda_toolkit toolkit_root cuda_version include_dir libcuda nvrtc cufft)
    if(EXISTS "${toolkit_root}/version.json")
        file(READ ${toolkit_root}/version.json version_json)
        string(JSON full_cuda_version GET ${version_json} cuda version)
    else()
        file(READ ${toolkit_root}/version.txt full_cuda_version)
    endif()
    string(REGEX MATCH "[0-9]*\\.[0-9]*" short_version ${full_cuda_version})
    set(${cuda_version} "${short_version}" PARENT_SCOPE)
    set(${include_dir} "${toolkit_root}/include" PARENT_SCOPE)
    set(${libcuda} "${toolkit_root}/lib64/stubs/libcuda.so" PARENT_SCOPE)
    set(${nvrtc} "${toolkit_root}/lib64/libnvrtc.so" PARENT_SCOPE)
    set(${cufft} "${toolkit_root}/lib64/libcufft.so" PARENT_SCOPE)
endfunction()


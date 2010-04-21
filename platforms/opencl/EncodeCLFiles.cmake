FILE(GLOB OPENCL_KERNELS ${CL_SOURCE_DIR}/kernels/*.cl)
SET(CL_FILE_DECLARATIONS)
SET(CL_FILE_DEFINITIONS)
CONFIGURE_FILE(${CL_SOURCE_DIR}/OpenCLKernelSources.cpp.in ${CL_KERNELS_CPP})
FOREACH(file ${OPENCL_KERNELS})
    # Load the file contents and process it.
    FILE(STRINGS ${file} file_content NEWLINE_CONSUME)
    STRING(REPLACE \\\n \\\\\n file_content "${file_content}")
    STRING(REPLACE \" \\\" file_content "${file_content}")
    STRING(REPLACE \n \\n\"\n\" file_content "${file_content}")

    # Determine a name for the variable that will contain this file's contents
    FILE(RELATIVE_PATH filename ${CL_SOURCE_DIR}/kernels ${file})
    STRING(LENGTH ${filename} filename_length)
    MATH(EXPR filename_length ${filename_length}-3)
    STRING(SUBSTRING ${filename} 0 ${filename_length} variable_name)

    # Record the variable declaration and definition.
    SET(CL_FILE_DECLARATIONS ${CL_FILE_DECLARATIONS}static\ const\ std::string\ ${variable_name};\n)
    FILE(APPEND ${CL_KERNELS_CPP} const\ string\ OpenCLKernelSources::${variable_name}\ =\ \"${file_content}\"\;\n)
ENDFOREACH(file)
CONFIGURE_FILE(${CL_SOURCE_DIR}/OpenCLKernelSources.h.in ${CL_KERNELS_H})

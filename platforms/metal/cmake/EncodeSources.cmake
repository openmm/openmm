# Metal Platform code: copyright (c) 2026 Chun-Chi Hung.
# Distributed under the GNU Lesser General Public License, version 3 or later.
# Embed unchanged source files. Runtime adaptation belongs to the Metal backend.

if(SOURCE_FILE)
    set(kernels "${SOURCE_FILE}")
else()
    file(GLOB kernels "${SOURCE_DIR}/*.${EXTENSION}")
endif()
file(MAKE_DIRECTORY "${OUTPUT_DIR}")
set(header "// Generated from kernel sources; do not edit.\n#pragma once\n#include <string>\nnamespace OpenMM {\nclass ${CLASS} {\npublic:\n")
set(source "// Generated from kernel sources; do not edit.\n#include \"${CLASS}.h\"\n")
foreach(kernel IN LISTS kernels)
    get_filename_component(name "${kernel}" NAME_WE)
    file(READ "${kernel}" content)
    string(APPEND header "    static const std::string ${name};\n")
    string(APPEND source "const std::string OpenMM::${CLASS}::${name} = R\"OPENMM_METAL(${content})OPENMM_METAL\";\n")
endforeach()
string(APPEND header "};\n}\n")
file(WRITE "${OUTPUT_DIR}/${CLASS}.h" "${header}")
file(WRITE "${OUTPUT_DIR}/${CLASS}.cpp" "${source}")

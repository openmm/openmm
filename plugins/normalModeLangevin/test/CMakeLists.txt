
add_executable(TestNormalModeLangevin TestNormalModeLangevin.cpp)
target_link_libraries(TestNormalModeLangevin OpenMM)
add_test(TestNormalModeLangevin ${EXECUTABLE_OUTPUT_PATH}/TestNormalModeLangevin)
set_tests_properties(TestNormalModeLangevin PROPERTIES ENVIRONMENT
    "OPENMM_PLUGIN_DIR=\"${LIBRARY_OUTPUT_DIRECTORY}\""
)

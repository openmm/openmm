
add_executable(TestNormalModeLangevin TestNormalModeLangevin.cpp)
target_link_libraries(TestNormalModeLangevin ${SHARED_TARGET})

# pass OPENMM_PLUGIN_DIR as command line argument
add_test(TestNormalModeLangevin 
    "${EXECUTABLE_OUTPUT_PATH}/TestNormalModeLangevin"
    "${test_plugin_dir}")

# does not work - nice try - maybe in cmake 2.8
set_tests_properties(TestNormalModeLangevin  PROPERTIES 
    ENVIRONMENT "OPENMM_PLUGIN_DIR=${test_plugin_dir}")

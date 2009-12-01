
# Copy plugin to test_plugin_dir
set(test_plugin_dir "${CMAKE_BINARY_DIR}/test_plugin_dir")
file(MAKE_DIRECTORY "${test_plugin_dir}")
get_target_property(old_loc NormalModeLangevin LOCATION)
get_filename_component(new_loc ${old_loc} NAME)
set(new_loc "${test_plugin_dir}/${new_loc}")
add_custom_command(
    DEPENDS ${old_loc} NormalModeLangevin
    OUTPUT ${new_loc}
    COMMAND ${CMAKE_COMMAND} -E copy ${old_loc} ${new_loc})
add_custom_target(CopyTestNmlPlugin ALL
    DEPENDS "${new_loc}"
    COMMENT "Copying normal mode langevin plugin for testing")

add_executable(TestNormalModeLangevin TestNormalModeLangevin.cpp)
target_link_libraries(TestNormalModeLangevin ${SHARED_TARGET})

# pass OPENMM_PLUGIN_DIR as command line argument
add_test(TestNormalModeLangevin 
    "${EXECUTABLE_OUTPUT_PATH}/TestNormalModeLangevin"
    "${test_plugin_dir}")

# does not work - nice try - maybe in cmake 2.8
set_tests_properties(TestNormalModeLangevin  PROPERTIES 
    ENVIRONMENT "OPENMM_PLUGIN_DIR=${test_plugin_dir}")



add_executable(TestNormalModeLangevin TestNormalModeLangevin.cpp)
target_link_libraries(TestNormalModeLangevin ${SHARED_TARGET})

add_test(TestNormalModeLangevin "${EXECUTABLE_OUTPUT_PATH}/TestNormalModeLangevin")


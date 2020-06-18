add_executable(vcdimtests vcdimtests.cpp)
target_link_libraries(vcdimtests gtest_main)
add_test(NAME vcdim COMMAND vcdimtests)

add_executable(projectiontests projectiontests.cpp)
target_link_libraries(projectiontests gtest_main)
add_test(NAME projection COMMAND projectiontests)

add_executable(projectiontests-inv projectiontests.cpp)
target_link_libraries(projectiontests-inv gtest_main)
target_compile_definitions(projectiontests-inv PUBLIC "WHEELWRIGHT_INVARIANTS")
add_test(NAME projection-inv COMMAND projectiontests-inv)

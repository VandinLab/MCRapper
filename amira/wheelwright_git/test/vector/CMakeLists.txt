add_executable(vectortests vectortests.cpp)
target_link_libraries(vectortests gtest_main)
add_test(NAME vector COMMAND vectortests)

add_executable(vectortests-inv vectortests.cpp)
target_link_libraries(vectortests-inv gtest_main)
target_compile_definitions(vectortests-inv PUBLIC "WHEELWRIGHT_INVARIANTS")
add_test(NAME vector-inv COMMAND vectortests-inv)

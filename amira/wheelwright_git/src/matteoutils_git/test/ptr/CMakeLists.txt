add_executable(ptrtests ptrtests.cpp)
target_link_libraries(ptrtests PRIVATE gtest_main)
add_test(NAME ptr COMMAND ptrtests)

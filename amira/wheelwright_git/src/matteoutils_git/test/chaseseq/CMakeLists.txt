add_executable(chaseseqtests chaseseqtests.cpp)
target_link_libraries(chaseseqtests PRIVATE gtest_main)
add_test(NAME chaseseq COMMAND chaseseqtests)

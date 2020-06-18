cmake_minimum_required(VERSION 3.8.0) # required because of compiler meta-features such as cxx_std_11
project(wheelwrighttests VERSION 0.01 LANGUAGES CXX)

find_package(Boost REQUIRED)

## Beginning of googletest setup, from
# https://github.com/google/googletest/tree/master/googletest#incorporating-into-an-existing-cmake-project
if (NOT TARGET gtest)
	# Download and unpack googletest at configure time
	configure_file(./CMakeLists.txt.in ${CMAKE_BINARY_DIR}/googletest-download/CMakeLists.txt)
	execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
		RESULT_VARIABLE result
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
	if(result)
		message(FATAL_ERROR "CMake step for googletest failed: ${result}")
	endif()
	execute_process(COMMAND ${CMAKE_COMMAND} --build .
		RESULT_VARIABLE result
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)
	if(result)
		message(FATAL_ERROR "Build step for googletest failed: ${result}")
	endif()
	# Prevent overriding the parent project's compiler/linker settings on
	# Windows
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
	# Add googletest directly to our build. This defines the gtest and
	# gtest_main targets.
	add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
		${CMAKE_BINARY_DIR}/googletest-build EXCLUDE_FROM_ALL)
endif()
## End of googletest setup

if(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
	if(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 4.0)
		set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3 -Og")
	else()
		set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3 -O0")
	endif()
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
	set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3 -fvar-tracking -fvar-tracking-assignments")
else()
	set(CMAKE_CXX_FLAGS_DEBUG "-g -fvar-tracking -fvar-tracking-assignments")
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL FreeBSD AND ${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
	if (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 5.0 AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 6.0)
		set(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath=/usr/local/lib/gcc5")
		set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-rpath=/usr/local/lib/gcc5")
	elseif (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 6.0 AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 7.0)
		set(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath=/usr/local/lib/gcc6")
		set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-rpath=/usr/local/lib/gcc6")
	elseif (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 7.0)
		set(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath=/usr/local/lib/gcc7")
		set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-rpath=/usr/local/lib/gcc7")
	else ()
		MESSAGE("You may want to set CMAKE_SHARED_LINKER_FLAGS to '-Wl,-rpath=/usr/local/lib/gccX', where 'X' is your g++ major version number.")
	endif()
endif()

set(WARNS_COMPILE_FLAGS -Wextra -Wall -Werror -pedantic-errors -Wformat=2
	-Wpointer-arith -Wcast-qual)
add_compile_options(${WARNS_COMPILE_FLAGS})
link_libraries(wheelwright)

add_subdirectory(bounds)
add_subdirectory(projection)
add_subdirectory(vector)

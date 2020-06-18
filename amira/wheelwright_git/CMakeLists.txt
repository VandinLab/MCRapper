cmake_minimum_required(VERSION 3.9.0) # required because of doxygen_add_docs

project(Wheelwright
	VERSION 0.01
	DESCRIPTION "A C++17 header-only library providing algorithms and data structures related to Rademacher averages"
	LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # must come after project()

find_package(Boost REQUIRED)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
	set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT FALSE)
endif()

add_subdirectory(src/matteoutils_git)

add_library(wheelwright INTERFACE)
target_include_directories(wheelwright INTERFACE ${Boost_INCLUDE_DIRS})
target_include_directories(wheelwright INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
	$<INSTALL_INTERFACE:include/wheelwright>
)
target_link_libraries(wheelwright INTERFACE matteoutils)
target_compile_features(wheelwright INTERFACE cxx_std_17)

if(NOT WHEELWRIGHT_INVARIANTS)
	if(CMAKE_BUILD_TYPE)
		if(${CMAKE_BUILD_TYPE} STREQUAL Debug)
			set(WHEELWRIGHT_INVARIANTS_DEF OFF)
		endif()
	else()
		set(WHEELWRIGHT_INVARIANTS_DEF ON)
	endif()
else()
	set(WHEELWRIGHT_INVARIANTS_DEF ${WHEELWRIGHT_INVARIANTS})
endif()
option(WHEELWRIGHT_INVARIANTS
	"Wheelwright: enable invariants checks" ${WHEELWRIGHT_INVARIANTS_DEF})
target_compile_definitions(wheelwright INTERFACE
	$<$<BOOL:${WHEELWRIGHT_INVARIANTS}>:WHEELWRIGHT_INVARIANTS>)

file(GLOB WHEELWRIGHT_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h")

option(WHEELWRIGHT_USE_OPENMP
	"Enable parallelism with OpenMP in Wheelwright"
	On)
if(WHEELWRIGHT_USE_OPENMP)
	find_package(OpenMP)
	if(OpenMP_CXX_FOUND)
		target_compile_options(wheelwright INTERFACE ${OpenMP_CXX_FLAGS})
		target_link_libraries(wheelwright INTERFACE ${OpenMP_CXX_FLAGS})
		target_link_libraries(wheelwright INTERFACE ${OpenMP_CXX_LIBRARIES})
	else()
		message(WARNING
			"Wheelwright: OpenMP parallelism requested, but no OpenMP support found. Falling back to no OpenMP parallelism."
			)
	endif()
endif()

option(WHEELWRIGHT_USE_NLOPT
	"Enable emp. Rade. Avg. bounds methods using NLOpt in Wheelwright."
	On)
if(WHEELWRIGHT_USE_NLOPT)
	find_package(NLopt)
	if(NOT NLopt_LIBRARY_DIRS)
		find_package(PkgConfig REQUIRED)
		pkg_search_module(NLopt REQUIRED nlopt)
	endif()
	target_compile_definitions(wheelwright INTERFACE WHEELWRIGHT_NLOPT)
	target_include_directories(wheelwright INTERFACE ${NLopt_INCLUDE_DIRS})
	target_link_libraries(wheelwright INTERFACE ${NLopt_LIBRARIES})
endif()

option(WHEELWRIGHT_BUILD_DOCS
	"Create and install the HTML API docs for Wheelwright (requires Doxygen)"
	Off)
if(WHEELWRIGHT_BUILD_DOCS)
	find_package(Doxygen REQUIRED OPTIONAL_COMPONENTS dot mscgen dia)
	set(DOXIGEN_EXCLUDE_SYMBOLS "*::impl::*")
	set(DOXYGEN_BUILTIN_STL_SUPPORT "YES")
	set(DOXYGEN_JAVADOC_AUTOBRIEF "YES")
	set(DOXYGEN_OPTIMIZE_FOR_C "YES")
	set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/docs/wheelwright")
	doxygen_add_docs(docs
		${WHEELWRIGHT_SRCS}
		${CMAKE_CURRENT_SOURCE_DIR}/README.md
		COMMENT "Generating Wheelwright docs with Doxygen"
	)
endif()

option(WHEELWRIGHT_BUILD_TESTS "Compile the tests for Wheelwright" Off)
if(WHEELWRIGHT_BUILD_TESTS)
	enable_testing()
	add_subdirectory(test)
endif()

install(FILES ${WHEELWRIGHT_SRCS} DESTINATION include/wheelwright)

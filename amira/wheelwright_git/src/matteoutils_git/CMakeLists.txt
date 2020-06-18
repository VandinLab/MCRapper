cmake_minimum_required(VERSION 3.9.0) # required because of doxygen_add_docs

project(MatteoUtils
		VERSION 0.01
		DESCRIPTION "A C++17 header-only library providing some utility functions"
		LANGUAGES CXX
		)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # must come after project()

find_package(Boost REQUIRED)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR} CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
	set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT FALSE)
endif()

add_library(matteoutils INTERFACE)
target_include_directories(matteoutils INTERFACE ${Boost_INCLUDE_DIRS})
target_include_directories(matteoutils INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
	$<INSTALL_INTERFACE:include/matteoutils>
)
file(GLOB MATTEOUTILS_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h")

option(MATTEOUTILS_USE_OPENMP
	"Enable parallelism with OpenMP in MatteoUtils"
	On)
if (MATTEOUTILS_USE_OPENMP)
	find_package(OpenMP)
	if (OpenMP_CXX_FOUND)
		target_compile_options(matteoutils INTERFACE ${OpenMP_CXX_FLAGS})
		target_link_libraries(matteoutils INTERFACE ${OpenMP_CXX_FLAGS})
		target_link_libraries(matteoutils INTERFACE ${OpenMP_CXX_LIBRARIES})
	else()
		message(WARNING
			"MatteoUtils: OpenMP parallelism requested, but no OpenMP support found. Falling back to no OpenMP parallelism."
			)
	endif()
endif()

option(MATTEOUTILS_ENABLE_INVARIANTS
	"Enable invariants checking in MatteoUtils"
	Off)
if (MATTEOUTILS_ENABLE_INVARIANTS)
	target_compile_definitions(matteoutils INTERFACE "MATTEOUTILS_INVARIANTS")
endif()

option(MATTEOUTILS_BUILD_DOCS
	"Create and install the HTML API docs for MatteoUtils (requires Doxygen)"
	Off)
if (MATTEOUTILS_BUILD_DOCS)
	find_package(Doxygen REQUIRED OPTIONAL_COMPONENTS dot mscgen dia)
	set(DOXIGEN_EXCLUDE_SYMBOLS "MatteoUtils::impl::*")
	set(DOXYGEN_BUILTIN_STL_SUPPORT "YES")
	set(DOXYGEN_JAVADOC_AUTOBRIEF "YES")
	set(DOXYGEN_OPTIMIZE_FOR_C "YES")
	set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/docs/matteoutils")
	doxygen_add_docs(docs
		${MATTEOUTILS_SRCS}
		${CMAKE_CURRENT_SOURCE_DIR}/README.md
		COMMENT "Generating MatteoUtils docs with Doxygen"
	)
endif()

option(BUILD_MATTEOUTILS_TESTS "Compile the tests for MatteoUtils" Off)
if (BUILD_MATTEOUTILS_TESTS)
	enable_testing()
	add_subdirectory(test)
endif()

install(FILES ${MATTEOUTILS_SRCS} DESTINATION include/matteoutils)

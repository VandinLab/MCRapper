# - Try to find NLopt
# Once done this will define
#  NLopt_FOUND - System has NLopt
#  NLopt_INCLUDE_DIRS - The NLopt include directories
#  NLopt_LIBRARIES - The libraries needed to use NLopt
#  NLopt_DEFINITIONS - Compiler switches required for using NLopt

find_package(PkgConfig)
pkg_check_modules(PC_NLopt QUIET nlopt)
set(NLopt_DEFINITIONS ${PC_NLopt_CFLAGS_OTHER})

find_path(NLopt_INCLUDE_DIR nlopt.hpp
          HINTS ${PC_NLopt_INCLUDEDIR} ${PC_NLopt_INCLUDE_DIRS}
          PATH_SUFFIXES NLopt)

find_library(NLopt_LIBRARY NAMES nlopt NLopt
             HINTS ${PC_NLopt_LIBDIR} ${PC_NLopt_LIBRARY_DIRS})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NLopt_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(NLopt  DEFAULT_MSG
                                  NLopt_LIBRARY NLopt_INCLUDE_DIR)

mark_as_advanced(NLopt_INCLUDE_DIR NLopt_LIBRARY)

set(NLopt_LIBRARIES ${NLopt_LIBRARY})
set(NLopt_INCLUDE_DIRS ${NLopt_INCLUDE_DIR})

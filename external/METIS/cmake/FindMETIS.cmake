###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2014 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
# - Find METIS include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(METIS
#               [REQUIRED]             # Fail with error if metis is not found
#              )
#
# This module finds headers and metis library.
# Results are reported in variables:
#  METIS_FOUND           - True if headers and requested libraries were found
#  METIS_INCLUDE_DIRS    - metis include directories
#  METIS_LIBRARY_DIRS    - Link directories for metis libraries
#  METIS_LIBRARIES       - metis component libraries to be linked
#  METIX_INTSIZE         - Number of octets occupied by a idx_t (IDXTYPEWIDTH)
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DMETIS_DIR=path/to/metis):
#  METIS_DIR             - Where to find the base directory of metis
#  METIS_INCDIR          - Where to find the header files
#  METIS_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: METIS_DIR, METIS_INCDIR, METIS_LIBDIR

#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013      Florent Pruvost
#
# This software is a computer program whose purpose is to process
# Matrices Over Runtime Systems @ Exascale (MORSE). More information
# can be found on the following website: http://www.inria.fr/en/teams/morse.
#
# This software is governed by the CeCILL-C license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL-C
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL-C license and that you accept its terms.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

if (NOT METIS_FOUND)
  set(METIS_DIR "" CACHE PATH "Installation directory of METIS library")
  if (NOT METIS_FIND_QUIETLY)
    message(STATUS "A cache variable, namely METIS_DIR, has been set to specify the install directory of METIS")
  endif()
endif()

# Looking for include
# -------------------

# Add system include paths to search include
# ------------------------------------------
unset(_inc_env)
set(ENV_METIS_DIR "$ENV{METIS_DIR}")
set(ENV_METIS_INCDIR "$ENV{METIS_INCDIR}")
if(ENV_METIS_INCDIR)
  list(APPEND _inc_env "${ENV_METIS_INCDIR}")
elseif(ENV_METIS_DIR)
  list(APPEND _inc_env "${ENV_METIS_DIR}")
  list(APPEND _inc_env "${ENV_METIS_DIR}/include")
  list(APPEND _inc_env "${ENV_METIS_DIR}/include/metis")
else()
  if(WIN32)
    string(REPLACE ":" ";" _inc_env "$ENV{INCLUDE}")
  else()
    string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{C_INCLUDE_PATH}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{CPATH}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{INCLUDE_PATH}")
    list(APPEND _inc_env "${_path_env}")
  endif()
endif()
list(APPEND _inc_env "${CMAKE_PLATFORM_IMPLICIT_INCLUDE_DIRECTORIES}")
list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
list(REMOVE_DUPLICATES _inc_env)


# Try to find the metis header in the given paths
# -------------------------------------------------
# call cmake macro to find the header path
if(METIS_INCDIR)
  set(METIS_metis.h_DIRS "METIS_metis.h_DIRS-NOTFOUND")
  find_path(METIS_metis.h_DIRS
    NAMES metis.h
    HINTS ${METIS_INCDIR})
else()
  if(METIS_DIR)
    set(METIS_metis.h_DIRS "METIS_metis.h_DIRS-NOTFOUND")
    find_path(METIS_metis.h_DIRS
      NAMES metis.h
      HINTS ${METIS_DIR}
      PATH_SUFFIXES "include" "include/metis")
  else()
    set(METIS_metis.h_DIRS "METIS_metis.h_DIRS-NOTFOUND")
    find_path(METIS_metis.h_DIRS
      NAMES metis.h
      HINTS ${_inc_env})
  endif()
endif()
mark_as_advanced(METIS_metis.h_DIRS)


# If found, add path to cmake variable
# ------------------------------------
if (METIS_metis.h_DIRS)
  set(METIS_INCLUDE_DIRS "${METIS_metis.h_DIRS}")
else ()
  set(METIS_INCLUDE_DIRS "METIS_INCLUDE_DIRS-NOTFOUND")
  if(NOT METIS_FIND_QUIETLY)
    message(STATUS "Looking for metis -- metis.h not found")
  endif()
endif()


# Looking for lib
# ---------------

# Add system library paths to search lib
# --------------------------------------
unset(_lib_env)
set(ENV_METIS_LIBDIR "$ENV{METIS_LIBDIR}")
if(ENV_METIS_LIBDIR)
  list(APPEND _lib_env "${ENV_METIS_LIBDIR}")
elseif(ENV_METIS_DIR)
  list(APPEND _lib_env "${ENV_METIS_DIR}")
  list(APPEND _lib_env "${ENV_METIS_DIR}/lib")
else()
  if(WIN32)
    string(REPLACE ":" ";" _lib_env "$ENV{LIB}")
  else()
    if(APPLE)
      string(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
    else()
      string(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
    endif()
    list(APPEND _lib_env "${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
    list(APPEND _lib_env "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
  endif()
endif()
list(REMOVE_DUPLICATES _lib_env)

# Try to find the metis lib in the given paths
# ----------------------------------------------
# call cmake macro to find the lib path
if(METIS_LIBDIR)
  set(METIS_metis_LIBRARY "METIS_metis_LIBRARY-NOTFOUND")
  find_library(METIS_metis_LIBRARY
    NAMES metis
    HINTS ${METIS_LIBDIR})
else()
  if(METIS_DIR)
    set(METIS_metis_LIBRARY "METIS_metis_LIBRARY-NOTFOUND")
    find_library(METIS_metis_LIBRARY
      NAMES metis
      HINTS ${METIS_DIR}
      PATH_SUFFIXES lib lib32 lib64)
  else()
    set(METIS_metis_LIBRARY "METIS_metis_LIBRARY-NOTFOUND")
    find_library(METIS_metis_LIBRARY
      NAMES metis
      HINTS ${_lib_env})
  endif()
endif()
mark_as_advanced(METIS_metis_LIBRARY)


# If found, add path to cmake variable
# ------------------------------------
if (METIS_metis_LIBRARY)
  get_filename_component(metis_lib_path "${METIS_metis_LIBRARY}" PATH)
  # set cmake variables
  set(METIS_LIBRARIES    "${METIS_metis_LIBRARY}")
  set(METIS_LIBRARY_DIRS "${metis_lib_path}")
else ()
  set(METIS_LIBRARIES    "METIS_LIBRARIES-NOTFOUND")
  set(METIS_LIBRARY_DIRS "METIS_LIBRARY_DIRS-NOTFOUND")
  if(NOT METIS_FIND_QUIETLY)
    message(STATUS "Looking for metis -- lib metis not found")
  endif()
endif ()

# check a function to validate the find
if(METIS_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # METIS
  if (METIS_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS  "${METIS_INCLUDE_DIRS}")
  endif()
  if (METIS_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${METIS_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${METIS_LIBRARIES}")
  # m
  find_library(M_LIBRARY NAMES m)
  mark_as_advanced(M_LIBRARY)
  if(M_LIBRARY)
    list(APPEND REQUIRED_LIBS "-lm")
  endif()

  # set required libraries for link
  set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
  set(CMAKE_REQUIRED_LIBRARIES)
  foreach(lib_dir ${REQUIRED_LIBDIRS})
    list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
  endforeach()
  list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
  string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

  # test link
  unset(METIS_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(METIS_NodeND METIS_WORKS)
  mark_as_advanced(METIS_WORKS)

  if(NOT METIS_WORKS)
    if(NOT METIS_FIND_QUIETLY)
      message(STATUS "Looking for METIS : test of METIS_NodeND with METIS library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(METIS_LIBRARIES)

if (METIS_LIBRARIES)
  list(GET METIS_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(METIS_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of METIS library" FORCE)
  else()
    set(METIS_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of METIS library" FORCE)
  endif()
endif()
mark_as_advanced(METIS_DIR)
mark_as_advanced(METIS_DIR_FOUND)

# Check the size of METIS_Idx
# ---------------------------------
set(CMAKE_REQUIRED_INCLUDES ${METIS_INCLUDE_DIRS})

include(CheckCSourceRuns)
#stdio.h and stdint.h should be included by metis.h directly
set(METIS_C_TEST_METIS_Idx_4 "
#include <stdio.h>
#include <stdint.h>
#include <metis.h>
int main(int argc, char **argv) {
  if (sizeof(idx_t) == 4)
    return 0;
  else
    return 1;
}
")

set(METIS_C_TEST_METIS_Idx_8 "
#include <stdio.h>
#include <stdint.h>
#include <metis.h>
int main(int argc, char **argv) {
  if (sizeof(idx_t) == 8)
    return 0;
  else
    return 1;
}
")
unset(METIS_Idx_4 CACHE)
unset(METIS_Idx_8 CACHE)
check_c_source_runs("${METIS_C_TEST_METIS_Idx_4}" METIS_Idx_4)
check_c_source_runs("${METIS_C_TEST_METIS_Idx_8}" METIS_Idx_8)
if(NOT METIS_Idx_4)
  if(NOT METIS_Idx_8)
    set(METIS_INTSIZE -1)
  else()
    set(METIS_INTSIZE 8)
  endif()
else()
  set(METIS_INTSIZE 4)
endif()
set(CMAKE_REQUIRED_INCLUDES "")

# check that METIS has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS DEFAULT_MSG
  METIS_LIBRARIES
  METIS_WORKS)
#
# TODO: Add possibility to check for specific functions in the library
#

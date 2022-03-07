# - Find METIS include dirs and libraries
#   Use this module by invoking find_package with the form:
#  find_package(METIS
#               [REQUIRED] # optional, if forced fail with error if not found
#              )
#
# This module finds headers and metis library.
# Results are reported in variables:
#  METIS_FOUND               - True if headers and requested libs were found
#  METIS_INCLUDE_DIRS        - metis include directories
#  METIS_LIBRARIES           - metis libraries to be linked
#
# The User can give a specific path where to find header and libraries adding
# cmake options at configure (f.e. with -DMETIS_DIR=path/to/metis):
#
#  METIS_DIR             - Where to find the base directory of metis, that is
#                         the root folder containing header "include" and
#                         libraries "lib" or "lib64" folders
#
# Advanced vars METIS_INCDIR and METIS_LIB are exposed also for cross-check
# purposes of the search.
#
# METIS_DIR_FOUND is exposed as cached variable to notify the GUI/ccmake User of CGNS
# successfull search
# BEWARE: tested OSs are Linux and Windows under MSYS2/MINGW environment.
#
#----------------------------------------------------------------------------
#
#  mimic
#
#  Optimad Engineering S.r.l. ("COMPANY") CONFIDENTIAL
#  Copyright (c) 2015-2022 Optimad Engineering S.r.l., All Rights Reserved.
#
#  --------------------------------------------------------------------------
#
#  NOTICE:  All information contained herein is, and remains the property
#  of COMPANY. The intellectual and technical concepts contained herein are
#  proprietary to COMPANY and may be covered by Italian and Foreign Patents,
#  patents in process, and are protected by trade secret or copyright law.
#  Dissemination of this information or reproduction of this material is
#  strictly forbidden unless prior written permission is obtained from
#  COMPANY. Access to the source code contained herein is hereby forbidden
#  to anyone except current COMPANY employees, managers or contractors who
#  have executed Confidentiality and Non-disclosure agreements explicitly
#  covering such access.
#
#  The copyright notice above does not evidence any actual or intended
#  publication or disclosure of this source code, which includes information
#  that is confidential and/or proprietary, and is a trade secret, of
#  COMPANY. ANY REPRODUCTION, MODIFICATION, DISTRIBUTION, PUBLIC PERFORMANCE,
#  OR PUBLIC DISPLAY OF OR THROUGH USE  OF THIS  SOURCE CODE  WITHOUT THE
#  EXPRESS WRITTEN CONSENT OF COMPANY IS STRICTLY PROHIBITED, AND IN
#  VIOLATION OF APPLICABLE LAWS AND INTERNATIONAL TREATIES. THE RECEIPT OR
#  POSSESSION OF THIS SOURCE CODE AND/OR RELATED INFORMATION DOES NOT CONVEY
#  OR IMPLY ANY RIGHTS TO REPRODUCE, DISCLOSE OR DISTRIBUTE ITS CONTENTS, OR
#  TO MANUFACTURE, USE, OR SELL ANYTHING THAT IT  MAY DESCRIBE, IN WHOLE OR
#  IN PART.
#
#----------------------------------------------------------------------------

if (NOT METIS_DIR)
  set(METIS_DIR "" CACHE PATH "Force manually a path to a custom METIS installation dir. Leave it empty for auto-search.")
  if (NOT METIS_FIND_QUIETLY)
    message(STATUS "A METIS_DIR cache variable is exposed to specify manually the install directory of METIS")
  endif()
endif()

# Find include directories
# ------------------------------------------
unset(includePathsEnv)

 #looking into some standard env variables for include.
string(REPLACE ":" ";" temp "$ENV{INCLUDE}")
list(APPEND includePathsEnv "${temp}")
if(NOT MINGW)
    string(REPLACE ":" ";" temp "$ENV{CPATH}")
    list(APPEND includePathsEnv "${temp}")
    string(REPLACE ":" ";" temp "$ENV{INCLUDE_PATH}")
    list(APPEND includePathsEnv "${temp}")
    string(REPLACE ":" ";" temp "$ENV{C_INCLUDE_PATH}")
    list(APPEND includePathsEnv "${temp}")
endif()
# adding further paths to the pot
list(APPEND includePathsEnv "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
#remove duplicates
list(REMOVE_DUPLICATES includePathsEnv)

# find include paths. If available METIS_DIR as Cmake variable search into it
# otherwise scan the pot of env paths
unset(METIS_INCDIR CACHE)

if(METIS_DIR)
    find_path(METIS_INCDIR
              NAMES metis.h
              HINTS ${METIS_DIR}
              PATH_SUFFIXES "include" "include/metis")
else()
    find_path(METIS_INCDIR
              NAMES metis.h
              HINTS ${includePathsEnv})
endif()

mark_as_advanced(METIS_INCDIR)

#Update on unsuccessful search of headers in REQUIRED mode
if (NOT METIS_INCDIR)
    if(NOT METIS_FIND_QUIETLY)
        message(STATUS "Looking for METIS headers: metis.h not found")
    endif()
endif()

#Same run but lookin for libraries.
unset(libPathsEnv)
if(MINGW)
    string(REPLACE ":" ";" libPathsEnv "$ENV{LIB}")
#elseif(APPLE)
#    string(REPLACE ":" ";" libPathsEnv "$ENV{DYLD_LIBRARY_PATH}")
else()
    string(REPLACE ":" ";" libPathsEnv "$ENV{LD_LIBRARY_PATH}")
endif()
list(APPEND libPathsEnv "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
#remove diplicates
list(REMOVE_DUPLICATES libPathsEnv)

# find libraries. If available METIS_DIR as Cmake variable search into it
# otherwise scan the pot of env paths
unset(METIS_LIB CACHE)

list(APPEND METIS_NAMES "metis")
if(METIS_DIR)
    find_library(METIS_LIB
        NAMES ${METIS_NAMES}
        HINTS ${METIS_DIR}
        PATH_SUFFIXES lib lib32 lib64
        )
else()
    find_library(METIS_LIB
        NAMES ${METIS_NAMES}
        HINTS ${libPathsEnv}
    )
endif()

mark_as_advanced(METIS_LIB)

#Update on unsuccessful search of libraries in REQUIRED mode
if (NOT METIS_LIB)
    if(NOT METIS_FIND_QUIETLY)
        message(STATUS "Looking for METIS library: no suitable shared or static library found")
    endif()
endif()

#retrieve the package version
if (METIS_INCDIR )
    file(READ "${METIS_INCDIR}/metis.h" _metis_version_header)

    string(REGEX MATCH "define[ \t]+METIS_VER_MAJOR[ \t]+([0-9]+)" _metis_major_version_match "${_metis_version_header}")
    set(METIS_MAJOR_VERSION "${CMAKE_MATCH_1}")

    string(REGEX MATCH "define[ \t]+METIS_VER_MINOR[ \t]+([0-9]+)" _metis_minor_version_match "${_metis_version_header}")
    set(METIS_MINOR_VERSION "${CMAKE_MATCH_1}")

    string(REGEX MATCH "define[ \t]+METIS_VER_SUBMINOR[ \t]+([0-9]+)" _metis_subminor_version_match "${_metis_version_header}")
    set(METIS_SUBMINOR_VERSION "${CMAKE_MATCH_1}")

    if(NOT METIS_MAJOR_VERSION)
        set(METIS_VERSION METIS_VERSION-NOTFOUND)
    else()
        set(METIS_VERSION ${METIS_MAJOR_VERSION}.${METIS_MINOR_VERSION}.${METIS_SUBMINOR_VERSION})
    endif()
else ()
  set(METIS_VERSION METIS_VERSION-NOTFOUND)
endif ()

##expose the METIS_DIR_FOUND to notify the User
if (METIS_LIB)
    list(GET METIS_LIB 0 fentry)
    get_filename_component(fentrypath "${fentry}" PATH)
    if (${fentrypath} MATCHES "(/lib(32|64)?$)")
        string(REGEX REPLACE "(/lib(32|64)?$)" "" temp "${fentrypath}")
        set(METIS_DIR_FOUND "${temp}" CACHE PATH "Final METIS root directory found" FORCE)
    else()
        set(METIS_DIR_FOUND "${fentrypath}" CACHE PATH "Final METIS root directory found" FORCE)
    endif()
else()
    set(METIS_DIR_FOUND METIS_DIR-NOTFOUND CACHE PATH "Final METIS root directory found" FORCE)
endif()

#mark_as_advanced(METIS_DIR)
#mark_as_advanced(METIS_DIR_FOUND)

# handle the QUIETLY and REQUIRED arguments and set METIS_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS REQUIRED_VARS METIS_LIB METIS_INCDIR
                                   VERSION_VAR METIS_VERSION
                               FAIL_MESSAGE "METIS not found. Try providing manually a valid installation path into METIS_DIR")

#final step fill metis official INCLUDE_DIRS and LIBRARIES
if (METIS_FOUND)
  set(METIS_LIBRARIES "${METIS_LIB}")
  set(METIS_INCLUDE_DIRS "${METIS_INCDIR}")
  if (NOT TARGET METIS::METIS)
    add_library(METIS::METIS UNKNOWN IMPORTED)
    set_target_properties(METIS::METIS PROPERTIES
      IMPORTED_LOCATION "${METIS_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCDIR}")
  endif ()
endif ()

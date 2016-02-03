#---------------------------------------------------------------------------
#
#  bitpit
#
#  Copyright (C) 2015-2016 OPTIMAD engineering Srl
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of bitbit.
#
#  bitpit is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License v3 (LGPL)
#  as published by the Free Software Foundation.
#
#  bitpit is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
#
#---------------------------------------------------------------------------*/

# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Add a target to generate API documentation with Doxygen
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" OFF)

IF(BUILD_DOCUMENTATION)
  FIND_PACKAGE(Doxygen)
  IF(NOT DOXYGEN_FOUND)
    MESSAGE(FATAL_ERROR "Doxygen is needed to build the documentation.")
  ENDIF()

  # Documentation user options
  option(DOC_EXTRACT_PRIVATE "If set, all private members are included in the documentation" OFF)
  IF (DOC_EXTRACT_PRIVATE)
    SET(DOXY_EXTRACT_PRIVATE "YES")
  ELSE ()
    SET(DOXY_EXTRACT_PRIVATE "NO")
  ENDIF ()

  # Doxygen options
  #
  # DOXY_INPUT_DIRECTORIES will be substitued in the Doxygen input file,
  # therefore its elements are separated by spaces and not by semicolons.
  if (ENABLE_MPI)
    SET(DOXY_ENABLE_MPI 1)
  else ()
    SET(DOXY_ENABLE_MPI 0)
  endif ()

  SET(DOXYFILE_IN          ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
  SET(DOXYFILE_CSS         ${CMAKE_CURRENT_SOURCE_DIR}/bitpitStyleSheet.css)
  SET(DOXYFILE             ${PROJECT_BINARY_DIR}/Doxyfile)
  SET(DOXY_HTML_INDEX_FILE ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
  SET(DOXY_OUTPUT_ROOT     ${CMAKE_CURRENT_BINARY_DIR})
  SET(DOXY_PREDEFINED      "ENABLE_MPI=${DOXY_ENABLE_MPI}")

  SET(DOXY_INPUT_DIRECTORIES "${CMAKE_SOURCE_DIR}/src/common \\\n      ${CMAKE_SOURCE_DIR}/examples \\\n      ${CMAKE_SOURCE_DIR}/doc/extra")
  SET(DOXY_INPUT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/extra/bitpit.hpp;${CMAKE_CURRENT_SOURCE_DIR}/extra/bitpitStyleGuide.hpp;${COMMON_SOURCES};${COMMON_HEADERS};${EXAMPLES_SOURCES}")
  foreach (MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
    isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
    if (MODULE_ENABLED)
       string(TOUPPER ${MODULE_NAME} UPPERCASE_MODULE_NAME)

       set(DOXY_INPUT_DIRECTORIES "${DOXY_INPUT_DIRECTORIES} \\\n      ${CMAKE_SOURCE_DIR}/src/${MODULE_NAME}")
       list(APPEND DOXY_INPUT_FILES ${${UPPERCASE_MODULE_NAME}_SOURCES} ${${UPPERCASE_MODULE_NAME}_HEADERS})
    endif ()
  endforeach ()

  # Create Doxygen input file
  if (NOT "${DOC_EXTRACT_PRIVATE_INTERNAL}" STREQUAL "${DOC_EXTRACT_PRIVATE}")
    if (EXISTS ${DOXY_HTML_INDEX_FILE})
      FILE (REMOVE ${DOXY_HTML_INDEX_FILE})
    endif ()
    set(DOC_EXTRACT_PRIVATE_INTERNAL ${DOC_EXTRACT_PRIVATE} CACHE INTERNAL "This is the value of the last time DOC_EXTRACT_PRIVATE was set successfully." FORCE)
  endif ()

  CONFIGURE_FILE( ${DOXYFILE_IN} ${DOXYFILE} @ONLY )

  # Targets
  ADD_CUSTOM_COMMAND( OUTPUT ${DOXY_HTML_INDEX_FILE}
                      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE}
                      # The following should be ${DOXYFILE} only but it
                      # will break the dependency.
                      # The optimal solution would be creating a
                      # custom_command for ${DOXYFILE} generation
                      # but I still have to figure out how...
                      MAIN_DEPENDENCY ${DOXYFILE} ${DOXYFILE_IN} ${DOXYFILE_CSS}
                      DEPENDS ${DOXYFILE} ${DOXYFILE_IN} ${DOXYFILE_CSS} ${DOXY_INPUT_FILES}
                      COMMENT "Generating HTML documentation")

  ADD_CUSTOM_TARGET( doc ALL DEPENDS ${DOXY_HTML_INDEX_FILE} )
  add_custom_target(clean-doc COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  INSTALL( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION doc/${CMAKE_PROJECT_NAME}-${BITPIT_VERSION} )
ENDIF()

#---------------------------------------------------------------------------
#
#  bitpit
#
#  Copyright (C) 2015-2021 OPTIMAD engineering Srl
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of bitpit.
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

FIND_PACKAGE(Doxygen)
IF(NOT DOXYGEN_FOUND)
MESSAGE(FATAL_ERROR "Doxygen is needed to build the documentation.")
ENDIF()

# Documentation user options
option(BITPIT_DOC_EXTRACT_PRIVATE "If set, all private members are included in the documentation" OFF)
IF (BITPIT_DOC_EXTRACT_PRIVATE)
SET(DOXY_EXTRACT_PRIVATE "YES")
ELSE ()
SET(DOXY_EXTRACT_PRIVATE "NO")
ENDIF ()

# Doxygen options
#
# DOXY_INPUT_DIRECTORIES will be substitued in the Doxygen input file,
# therefore its elements are separated by spaces and not by semicolons.
if (BITPIT_ENABLE_MPI)
SET(DOXY_BITPIT_ENABLE_MPI 1)
else ()
SET(DOXY_BITPIT_ENABLE_MPI 0)
endif ()

SET(DOXYFILE_IN          ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
SET(DOXYFILE_CSS         ${CMAKE_CURRENT_SOURCE_DIR}/bitpitStyleSheet.css)
SET(DOXYFILE             ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
SET(DOXY_HTML_INDEX_FILE ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
SET(DOXY_OUTPUT_ROOT     ${CMAKE_CURRENT_BINARY_DIR})
SET(DOXY_PREDEFINED      "\"BITPIT_ENABLE_MPI=${DOXY_BITPIT_ENABLE_MPI}\" \"BITPIT_DEPRECATED(x)=x\"")

SET(DOXY_HTML_HEADER           "" CACHE INTERNAL "Header file for Doxygen documentation")
SET(DOXY_HTML_FOOTER           "" CACHE INTERNAL "Footer file for Doxygen documentation")
SET(DOXY_HTML_EXTRA_STYLESHEET "${CMAKE_SOURCE_DIR}/doc/bitpitStyleSheet.css" CACHE INTERNAL "Stylesheet file for Doxygen documentation")
SET(DOXY_DISABLE_INDEX         "YES" CACHE INTERNAL "Controls if the index should be disbaled")
SET(DOXY_GENERATE_TREEVIEW     "YES" CACHE INTERNAL "Controls if a tree-like index structure should be generated")

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
if (NOT "${BITPIT_DOC_EXTRACT_PRIVATE_INTERNAL}" STREQUAL "${BITPIT_DOC_EXTRACT_PRIVATE}")
if (EXISTS ${DOXY_HTML_INDEX_FILE})
    FILE (REMOVE ${DOXY_HTML_INDEX_FILE})
endif ()
set(BITPIT_DOC_EXTRACT_PRIVATE_INTERNAL ${BITPIT_DOC_EXTRACT_PRIVATE} CACHE INTERNAL "This is the value of the last time BITPIT_DOC_EXTRACT_PRIVATE was set successfully." FORCE)
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

INSTALL( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION ${CMAKE_INSTALL_DOCDIR}/${BITPIT_VERSION} )

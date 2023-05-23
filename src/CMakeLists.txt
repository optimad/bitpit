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
#----------------------------------------------------------------------------

#------------------------------------------------------------------------------------#
# Functions
#------------------------------------------------------------------------------------#

# Add a module to the project
#
# The module will be added to the project as a library built from the source
# files (*.cpp) available in the current directory.
#
# All configurable header files (*.hpp.in) and files with configurable template
# implementations (*.tpp.in) will be automatically configured and copied into
# the build directory.
#
# All header files (*.hpp) and files with template implementations (*.tpp) in
# the current directory and in the build directory will be installed inside
# the directory "include/${PROJECT_NAME}". Only files automatically generated
# by this function (e.g., configurable header files) or files created before
# calling this function will be installed.
function(configureModule MODULE_NAME)
    string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)

    # Configure compilation
    addModuleIncludeDirectories(${MODULE_NAME})

    # Configure headers
    file(GLOB CONFIGURABLE_HEADER_FILES "*.hpp.in" "*.tpp.in")
    foreach(CONFIGURABLE_HEADER_FILE IN LISTS CONFIGURABLE_HEADER_FILES)
        get_filename_component(CONFIGURED_HEADER_FILE ${CONFIGURABLE_HEADER_FILE} NAME)
        string(REGEX REPLACE "\\.[^.]*$" "" CONFIGURED_HEADER_FILE ${CONFIGURED_HEADER_FILE})
        set(CONFIGURED_HEADER_FILE "${PROJECT_BINARY_DIR}/src/${MODULE_NAME}/${CONFIGURED_HEADER_FILE}")
        configure_file("${CONFIGURABLE_HEADER_FILE}" "${CONFIGURED_HEADER_FILE}")
    endforeach()

    # Configure targets
    file(GLOB SOURCE_FILES "*.cpp")
    set(${UPPER_MODULE_NAME}_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of ${MODULE_NAME} module" FORCE)
    unset(SOURCE_FILES)

    if (NOT "${${UPPER_MODULE_NAME}_SOURCES}" STREQUAL "")
        add_library(${UPPER_MODULE_NAME}_TARGET_OBJECT OBJECT ${${UPPER_MODULE_NAME}_SOURCES})
        set_target_properties(${UPPER_MODULE_NAME}_TARGET_OBJECT PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})
        set_lto_property(${UPPER_MODULE_NAME}_TARGET_OBJECT)

        target_sources(${BITPIT_LIBRARY} PRIVATE $<TARGET_OBJECTS:${UPPER_MODULE_NAME}_TARGET_OBJECT>)
    endif ()

    # Configure installation
    file(GLOB HEADER_FILES "*.hpp" "*.tpp")
    file(GLOB CONFIGURED_HEADER_FILES "${PROJECT_BINARY_DIR}/src/${MODULE_NAME}/*.hpp" "${PROJECT_BINARY_DIR}/src/${MODULE_NAME}/*.tpp")
    set(${UPPER_MODULE_NAME}_HEADERS "${HEADER_FILES}" "${CONFIGURED_HEADER_FILES}" CACHE INTERNAL "Headers of ${MODULE_NAME} module" FORCE)
    unset(HEADER_FILES)
    unset(CONFIGURED_HEADER_FILES)

    if (NOT "${${UPPER_MODULE_NAME}_HEADERS}" STREQUAL "")
        install(FILES ${${UPPER_MODULE_NAME}_HEADERS} DESTINATION include/${PROJECT_NAME})
    endif ()

    # Configure source formatting
    if (CLANG_FORMAT_EXECUTABLE)
        set(MODULE_CLANG_FORMAT_TARGET "clang-format-${MODULE_NAME}")
        ClangFormatAddTarget("${MODULE_CLANG_FORMAT_TARGET}")
        ClangFormatTargetSources("${MODULE_CLANG_FORMAT_TARGET}" "${${UPPER_MODULE_NAME}_SOURCES}")
        ClangFormatTargetSources("${MODULE_CLANG_FORMAT_TARGET}" "${${UPPER_MODULE_NAME}_HEADERS}")
        add_dependencies(${CLANG_FORMAT_TARGET} "${MODULE_CLANG_FORMAT_TARGET}")
    endif()

endfunction()

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

# Modules
foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
    isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
    if (MODULE_ENABLED)
        add_subdirectory(${MODULE_NAME})
    endif ()
endforeach()

#------------------------------------------------------------------------------------#
# Targets
#------------------------------------------------------------------------------------#
add_custom_target(library DEPENDS ${BASE_LIBRARY})
add_custom_target(clean-library COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

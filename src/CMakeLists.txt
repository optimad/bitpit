# ---------------------------------------------------------------------------
#
# bitpit
#
# Copyright (C) 2015-2021 OPTIMAD engineering Srl
#
# -------------------------------------------------------------------------
# License
# This file is part of bitpit.
#
# bitpit is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License v3 (LGPL)
# as published by the Free Software Foundation.
#
# bitpit is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with bitpit. If not, see <http://www.gnu.org/licenses/>.
#
# ----------------------------------------------------------------------------

# ------------------------------------------------------------------------------------#
# Functions
# ------------------------------------------------------------------------------------#

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

    # Configure headers
    file(GLOB CONFIGURABLE_HEADER_FILES "*.hpp.in" "*.tpp.in")

    foreach(CONFIGURABLE_HEADER_FILE IN LISTS CONFIGURABLE_HEADER_FILES)
        get_filename_component(CONFIGURED_HEADER_FILE ${CONFIGURABLE_HEADER_FILE} NAME)
        string(REGEX REPLACE "\\.[^.]*$" "" CONFIGURED_HEADER_FILE ${CONFIGURED_HEADER_FILE})
        set(CONFIGURED_HEADER_FILE "${PROJECT_BINARY_DIR}/src/${MODULE_NAME}/${CONFIGURED_HEADER_FILE}")
        configure_file("${CONFIGURABLE_HEADER_FILE}" "${CONFIGURED_HEADER_FILE}")
    endforeach()

    # Get sources
    file(GLOB SOURCE_FILES "*.cpp")
    set(${UPPER_MODULE_NAME}_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of ${MODULE_NAME} module" FORCE)
    unset(SOURCE_FILES)

    # Get headers
    file(GLOB HEADER_FILES "*.hpp" "*.tpp")
    file(GLOB CONFIGURED_HEADER_FILES "${PROJECT_BINARY_DIR}/src/${MODULE_NAME}/*.hpp" "${PROJECT_BINARY_DIR}/src/${MODULE_NAME}/*.tpp")
    set(${UPPER_MODULE_NAME}_HEADERS "${HEADER_FILES}" "${CONFIGURED_HEADER_FILES}" CACHE INTERNAL "Headers of ${MODULE_NAME} module" FORCE)
    unset(HEADER_FILES)
    unset(CONFIGURED_HEADER_FILES)

    # Configure target
    #
    # For header-only no target is created, the headers are simply added to the list of bitpit
    # public headers.
    if(NOT "${${UPPER_MODULE_NAME}_SOURCES}" STREQUAL "")
        add_library(${UPPER_MODULE_NAME}_TARGET_OBJECT OBJECT)
        target_sources(${UPPER_MODULE_NAME}_TARGET_OBJECT PRIVATE ${${UPPER_MODULE_NAME}_SOURCES})

        target_compile_features(${UPPER_MODULE_NAME}_TARGET_OBJECT PUBLIC cxx_std_11)
        set_target_properties(${UPPER_MODULE_NAME}_TARGET_OBJECT PROPERTIES CXX_STANDARD 11)
        set_target_properties(${UPPER_MODULE_NAME}_TARGET_OBJECT PROPERTIES CXX_STANDARD_REQUIRED ON)
        set_target_properties(${UPPER_MODULE_NAME}_TARGET_OBJECT PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})
        set_lto_property(${UPPER_MODULE_NAME}_TARGET_OBJECT)
        addModuleIncludeDirectories(${UPPER_MODULE_NAME}_TARGET_OBJECT ${MODULE_NAME})

        get_target_property(BITPIT_COMPILE_DEFINITIONS ${BITPIT_LIBRARY} COMPILE_DEFINITIONS)
        target_compile_definitions(${UPPER_MODULE_NAME}_TARGET_OBJECT PUBLIC ${BITPIT_COMPILE_DEFINITIONS})

        target_sources(${BITPIT_LIBRARY} PRIVATE $<TARGET_OBJECTS:${UPPER_MODULE_NAME}_TARGET_OBJECT>)
    endif()

    if(NOT "${${UPPER_MODULE_NAME}_HEADERS}" STREQUAL "")
        if(NOT "${${UPPER_MODULE_NAME}_SOURCES}" STREQUAL "")
            target_sources(${UPPER_MODULE_NAME}_TARGET_OBJECT PRIVATE ${${UPPER_MODULE_NAME}_HEADERS})
            set_target_properties(${UPPER_MODULE_NAME}_TARGET_OBJECT PROPERTIES PUBLIC_HEADER "${${UPPER_MODULE_NAME}_HEADERS}")
        endif()

        get_target_property(BITPIT_PUBLIC_HEADERS ${BITPIT_LIBRARY} PUBLIC_HEADER)

        if(NOT "${BITPIT_PUBLIC_HEADERS}" STREQUAL "" AND NOT "${BITPIT_PUBLIC_HEADERS}" STREQUAL "BITPIT_PUBLIC_HEADERS-NOTFOUND")
            set(BITPIT_PUBLIC_HEADERS "${BITPIT_PUBLIC_HEADERS};${${UPPER_MODULE_NAME}_HEADERS}")
        else()
            set(BITPIT_PUBLIC_HEADERS "${${UPPER_MODULE_NAME}_HEADERS}")
        endif()

        set_target_properties(${BITPIT_LIBRARY} PROPERTIES PUBLIC_HEADER "${BITPIT_PUBLIC_HEADERS}")
    endif()

    # Configure source formatting
    if(CLANG_FORMAT_EXECUTABLE)
        set(MODULE_CLANG_FORMAT_TARGET "clang-format-${MODULE_NAME}")
        ClangFormatAddTarget("${MODULE_CLANG_FORMAT_TARGET}")
        ClangFormatTargetSources("${MODULE_CLANG_FORMAT_TARGET}" "${${UPPER_MODULE_NAME}_SOURCES}")
        ClangFormatTargetSources("${MODULE_CLANG_FORMAT_TARGET}" "${${UPPER_MODULE_NAME}_HEADERS}")
        add_dependencies(${CLANG_FORMAT_TARGET} "${MODULE_CLANG_FORMAT_TARGET}")
    endif()
endfunction()

# ------------------------------------------------------------------------------------#
# Subdirectories
# ------------------------------------------------------------------------------------#

# Modules
foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
    isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)

    if(MODULE_ENABLED)
        add_subdirectory(${MODULE_NAME})
    endif()
endforeach()

# ------------------------------------------------------------------------------------#
# Targets
# ------------------------------------------------------------------------------------#
add_custom_target(library DEPENDS ${BASE_LIBRARY})
add_custom_target(clean-library COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

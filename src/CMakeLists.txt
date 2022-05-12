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
# CMake settings
#------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.10)
initialize_lto_property()

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
        set(${UPPER_MODULE_NAME}_TARGET_OBJECT "${UPPER_MODULE_NAME}_TARGET_OBJECT")
        add_library(${${UPPER_MODULE_NAME}_TARGET_OBJECT} OBJECT ${${UPPER_MODULE_NAME}_SOURCES})

        initialize_lto_property()
        set_lto_property(${${UPPER_MODULE_NAME}_TARGET_OBJECT})
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
    set(MODULE_CLANG_FORMAT_TARGET "clang-format-${MODULE_NAME}")
    ClangFormatAddTarget("${MODULE_CLANG_FORMAT_TARGET}")
    ClangFormatTargetSources("${MODULE_CLANG_FORMAT_TARGET}" "${${UPPER_MODULE_NAME}_SOURCES}")
    ClangFormatTargetSources("${MODULE_CLANG_FORMAT_TARGET}" "${${UPPER_MODULE_NAME}_HEADERS}")
    add_dependencies(${CLANG_FORMAT_TARGET} "${MODULE_CLANG_FORMAT_TARGET}")

endfunction()

#------------------------------------------------------------------------------------#
# Library properties
#------------------------------------------------------------------------------------#

# Library type
if (BITPIT_BUILD_SHARED_LIBS)
    set(LIBRARY_TYPE "SHARED")
    set(PIC_FLAG TRUE)
else()
    set(LIBRARY_TYPE "STATIC")
    set(PIC_FLAG FALSE)
endif()

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

if (DEFINED MODULE_SOURCES)
    unset(${MODULE_SOURCES})
    unset(${MODULE_DEFINITIONS})
endif ()

# Modules
foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
    isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
    if (MODULE_ENABLED)
        string(TOUPPER ${MODULE_NAME} UPPERCASE_MODULE_NAME)
        add_subdirectory(${MODULE_NAME})
        if (TARGET ${UPPERCASE_MODULE_NAME}_TARGET_OBJECT)
            set(MODULE_SOURCES ${MODULE_SOURCES} $<TARGET_OBJECTS:${UPPERCASE_MODULE_NAME}_TARGET_OBJECT>)
            set_target_properties(${UPPERCASE_MODULE_NAME}_TARGET_OBJECT PROPERTIES POSITION_INDEPENDENT_CODE ${PIC_FLAG})
        endif ()
        if (DEFINED ${UPPERCASE_MODULE_NAME}_DEFINITIONS)
            set(MODULE_DEFINITIONS ${MODULE_DEFINITIONS} ${${UPPERCASE_MODULE_NAME}_DEFINITIONS})
        endif ()
    endif ()
endforeach()

#------------------------------------------------------------------------------------#
# Library
#------------------------------------------------------------------------------------#

# NOTE: to get lib files, enable both BITPIT_BUILD_SHARED_LIBS and CMAKE_GNUtoMS flags.

set(BITPIT_LIBRARY ${PROJECT_NAME} CACHE INTERNAL "bitpit Base Library" FORCE)

add_library(${BITPIT_LIBRARY} ${LIBRARY_TYPE} ${MODULE_SOURCES})

set_lto_property(${BITPIT_LIBRARY})

if (BITPIT_ENABLE_MPI)
    set_target_properties(${BITPIT_LIBRARY} PROPERTIES DEBUG_POSTFIX "_MPI_D")
    set_target_properties(${BITPIT_LIBRARY} PROPERTIES RELWITHDEBINFO_POSTFIX "_MPI_D")
    set_target_properties(${BITPIT_LIBRARY} PROPERTIES RELEASE_POSTFIX "_MPI")
else()
    set_target_properties(${BITPIT_LIBRARY} PROPERTIES DEBUG_POSTFIX "_D")
    set_target_properties(${BITPIT_LIBRARY} PROPERTIES RELWITHDEBINFO_POSTFIX "_D")
endif()

set_target_properties(${BITPIT_LIBRARY} PROPERTIES POSITION_INDEPENDENT_CODE ${PIC_FLAG})
set_target_properties(${BITPIT_LIBRARY} PROPERTIES VERSION "${BITPIT_MAJOR_VERSION}.${BITPIT_MINOR_VERSION}.${BITPIT_PATCH_VERSION}"
                                                 SOVERSION "${BITPIT_MAJOR_VERSION}")

# Windows shared libraries have to be linked against their dependencies.
if (MINGW AND BITPIT_BUILD_SHARED_LIBS)
    target_link_libraries(${BITPIT_LIBRARY} ${BITPIT_EXTERNAL_LIBRARIES})
endif()

# Install target
#
# On Windows, shared libraries have to be installed in the binary directory.
if (WIN32)
    install(TARGETS ${BITPIT_LIBRARY} ARCHIVE DESTINATION lib LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
else()
    install(TARGETS ${BITPIT_LIBRARY} DESTINATION lib)
endif()

# Definitions
set(BITPIT_DEFINITIONS ${MODULE_DEFINITIONS} CACHE INTERNAL "Pre-processor definitions needed by the Base Library" FORCE)

# Targets
add_custom_target(library DEPENDS ${BASE_LIBRARY})
add_custom_target(clean-library COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

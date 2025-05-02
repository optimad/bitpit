#---------------------------------------------------------------------------
#
#  bitpit
#
#  Copyright (C) 2015-2025 OPTIMAD engineering Srl
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

# Collects Windows runtime DLLs for the specified target
#
# \param TARGET_NAME is the name of the target for which DLLs will be collected
# \param RUNTIME_DLLS_EXCLUDE_PATTERNS is a list of patterns to exclude from the DLL collection
# \param RUNTIME_DLLS_INSTALL_DIR is the directory where the runtime DLLs will be installed,
# the directory can be an absolute path or a path relative to the CMake install prefix
# \param COLLECTION_TIME specifies when the command should be executed: "POST_BUILD" or "INSTALL"
function(windows_runtime_collect_dlls TARGET_NAME RUNTIME_DLLS_EXCLUDE_PATTERNS RUNTIME_DLLS_INSTALL_DIR COLLECTION_TIME)
    if(NOT IS_ABSOLUTE "${RUNTIME_DLLS_INSTALL_DIR}")
        get_filename_component(ABSOLUTE_RUNTIME_DLLS_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${RUNTIME_DLLS_INSTALL_DIR}" ABSOLUTE)
    else()
        set(ABSOLUTE_RUNTIME_DLLS_INSTALL_DIR "${RUNTIME_DLLS_INSTALL_DIR}")
    endif()

    SET(COLLECT_RUNTIME_DLLS_COMMAND
        ${CMAKE_COMMAND} -E echo "Collecting runtime DLLs for target ${TARGET_NAME}" &&
        ${CMAKE_COMMAND} -E echo "Runtime DLLs will be copied to ${ABSOLUTE_RUNTIME_DLLS_INSTALL_DIR}" &&
        ${CMAKE_COMMAND}
        -D filepath_list="$<JOIN:$<TARGET_RUNTIME_DLLS:${TARGET_NAME}>,$<SEMICOLON>>"
        -D destination_directory="${ABSOLUTE_RUNTIME_DLLS_INSTALL_DIR}"
        -D exclude_filename_regex_list="${RUNTIME_DLLS_EXCLUDE_PATTERNS}"
        -P ${PROJECT_SOURCE_DIR}/cmake/CopyFiles.cmake
    )

    if(COLLECTION_TIME STREQUAL "POST_BUILD")
        add_custom_command(TARGET ${TARGET_NAME} POST_BUILD COMMAND ${COLLECT_RUNTIME_DLLS_COMMAND})
    elseif(COLLECTION_TIME STREQUAL "INSTALL")
        set(COLLECT_RUNTIME_DLLS_TARGET ${TARGET_NAME}_windows_runtime_dlls_collection)
        add_custom_target(${COLLECT_RUNTIME_DLLS_TARGET} COMMAND ${COLLECT_RUNTIME_DLLS_COMMAND})
        install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} --build . --config $<CONFIG> --target ${COLLECT_RUNTIME_DLLS_TARGET})")
    else()
        message(FATAL_ERROR "Invalid COLLECTION_TIME specified. Use 'POST_BUILD' or 'INSTALL'.")
    endif()
endfunction()

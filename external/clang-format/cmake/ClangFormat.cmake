#---------------------------------------------------------------------------
#
#  bitpit
#
#  Copyright (C) 2015-2022 OPTIMAD engineering Srl
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

#[[

This file allows to format source files using clang-format.

This module provides the following function:

::

  ClangFormatAddTarget(<target>)

This function will create an empty clang-format target. Once the target has been
created, it is possible to add sources files to be formatted by the target using
the function ClangFormatTargetSources. Files can be formatted invoke the target
with ``make <target>``.

  ClangFormatAddTarget(<target> <files>)

This function will add the specified source files to clanf-format target
previously created.

Example usage:

.. code-block:: cmake

  include(ClangFormat)
  file(GLOB_RECURSE SOURCE_FILES *.cpp *.h *.hpp *.c)
  ClangFormatAddTarget("clang-format")
  kde_clang_format("clang-format" ${SOURCE_FILES})

The ``.clang-format`` file should be provided by the repository.

To exclude directories from the formatting add a ``.clang-format``
file in the directory with the following contents:

.. code-block:: yaml

   DisableFormat: true
   SortIncludes: false

]]

# Find clang-format executabel
find_program(CLANG_FORMAT_EXECUTABLE clang-format)

# Add clang-format target
function(ClangFormatAddTarget TARGET)
    # Early return if clang-formt is not available
    if (NOT CLANG_FORMAT_EXECUTABLE)
        return()
    endif()

    # Create format target
    add_custom_target(${TARGET} COMMENT "Formatting sources using ${CLANG_FORMAT_EXECUTABLE}...")
endfunction()

# Specifies sources to be formatted by the target
function(ClangFormatTargetSources TARGET FILE_PATHS)
    # Early return if clang-formt is not available
    if (NOT TARGET ${TARGET})
        message(FATAL_ERROR "Target ${TARGET} is not valid.")
    endif()

    # Early return if clang-formt is not available
    if (NOT CLANG_FORMAT_EXECUTABLE)
        return()
    endif()

    # Add files to the target
    get_filename_component(FULL_BINARY_DIR ${CMAKE_BINARY_DIR} REALPATH)
    foreach(FILE_PATH ${FILE_PATHS})
        # Skip files inside the buil directory
        get_filename_component(FULL_FILE_PATH ${FILE_PATH} REALPATH)
        string(FIND ${FULL_FILE_PATH} ${FULL_BINARY_DIR} BINARY_DIR_INDEX)
        if(BINARY_DIR_INDEX EQUAL 0)
            continue()
        endif()

        # Add file
        add_custom_command(
            TARGET ${TARGET}
            COMMAND ${CLANG_FORMAT_EXECUTABLE} -style=file -i ${FULL_FILE_PATH}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Formatting ${FULL_FILE_PATH}..."
        )
    endforeach()
endfunction()

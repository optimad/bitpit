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

#---------------------------------------------------------------------------
#       CMake script used to copy dependencies into a designated folder
#---------------------------------------------------------------------------
#       What does this script do?
#
#           It copies the dependencies of the indicated executables and libraries into the
#           designated folder.
#
#           Should one specify some files to be excluded, the dependencies with the same name
#           as one of the specified files will be ignored and not copied (e.g. if, as input,
#           this module receives a path to boost_unit_test_framework-vc142-mt-gd-x64-1_77.dll,
#           then it will ignore any other versions of that same library it finds during the
#           search).
#
#           In addition to searching the usual OS-specific directories for libraries, it also
#           searches the directories specified in the additional_search_directories variable.
#
#           This module cycles through all of the excluded file paths, and appends the folder
#           that contains each one to the additional_search_directories.
#
#           It automatically excludes all libraries. You must explicitly include libraries you
#           wish to copy. You must also include libraries for which you wish to find transitive
#           dependencies otherwise they will not be considered. When the search encounters a
#           library/exe not in the include list, it ignores that file and continues the search.
#           Therefore all of its dependencies are not included.
#
#       When to use it?
#
#           During install, as a custom command, or as a custom target. Because it makes use of
#           GET_RUNTIME_DEPENDENCIES, it cannot be used during cmake configuration.
#           Example:
#               add_custom_target(copy_dependencies COMMAND
#                       ${CMAKE_COMMAND}
#                       -D executable_list="$<TARGET_FILE:Executable>"
#                       -D library_list="$<TARGET_FILE:Library>"
#                       -D copy_to_directory="${CMAKE_BINARY_DIR}/$<CONFIG>"
#                       -D additional_search_directories="${OMP_DLL_DIR}"
#                       -D exclude_filepath_list="$<TARGET_RUNTIME_DLLS:RomboxTests>"
#                       -D include_filename_patterns=".*rombox.*$<SEMICOLON>.*romby.*$<SEMICOLON>
#                                                   .*boost.*$<SEMICOLON>.*vtk.*"
#                       -P ${PROJECT_SOURCE_DIR}/cmake/copy_dependencies.cmake
#                )
#
#       What are the inputs?
#
#           NOTE: Almost all these inputs are lists. Be sure to use $<SEMICOLON> instead of ";"
#                 so that the lists get resolved properly during the build phase. You can always
#                 perform a replace with string (REPLACE ";" "$<SEMICOLON>" my_list "${my_list}")
#
#           NOTE: Include in the include_filename_patterns your libraries that you need to include in
#                 the search. e.g. Calling this script on a unit tests exe doesn't really need the
#                 library since it's being built as well, but if you don't include the library,
#                 then that library's dependencies will not be included.
#
#           executable_list -               A list of executable files to find dependencies for
#           library_list -                  A list of libraries to find dependencies for
#           include_filename_patterns -     A list of regex patterns that identify dependency names
#                                           that should be included in the search. At least one
#                                           pattern should be specified otherwise an error is
#                                           raised
#           exclude_filepath_list -         A list of absolute paths of dependencies that are
#                                           already resolved and should not be copied if their
#                                           filename matches one of the requested search patterns.
#                                           Any other dependency with the same name found during
#                                           the search will be ignored.
#           copy_to_directory -             The directory where the dependencies should be copied
#                                           If not specified, position of first target.
#           additional_search_directories - Extra directories where cmake should search for
#                                           dependencies
#
#---------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.16)

if (NOT include_filename_patterns)
    message(FATAL_ERROR "You have specified no dependencies to include. That means that no dependencies will be copied. Calling this script is pointless")
endif()
if (NOT additional_search_directories)
    set(additional_search_directories " ")
endif()
if((NOT executable_list) AND (NOT library_list))
    message(FATAL_ERROR "You have specified neither a executable nor a library list. This script is therefore pointless.")
endif()

if(NOT copy_to_directory)
    if (executable_list)
        foreach(executable ${executable_list})
            get_filename_component(folder ${executable} DIRECTORY)
            continue()
        endforeach()
    else()
        foreach(executable ${library_list})
            get_filename_component(folder ${executable} DIRECTORY)
            continue()
        endforeach()
    endif()
    set(copy_to_directory "${folder}/")
    message(STATUS "Because you didn't specify a directory, the dependencies will be copied to folder ${copy_to_directory}")
endif()

set(exclude_filename_list "")
foreach (filepath ${exclude_filepath_list})
    get_filename_component(directory ${filepath} DIRECTORY)
    list(APPEND additional_search_directories ${directory})

    get_filename_component(filename ${filepath} NAME)
    string(TOLOWER "${filename}" filename)
    list(APPEND exclude_filename_list "${filename}")
endforeach()

if (executable_list AND library_list)
    file(GET_RUNTIME_DEPENDENCIES
            EXECUTABLES ${executable_list}
            LIBRARIES ${library_list}
            RESOLVED_DEPENDENCIES_VAR resolved_dependencies
            UNRESOLVED_DEPENDENCIES_VAR unresolved_dependencies
            PRE_EXCLUDE_REGEXES ".*"
            PRE_INCLUDE_REGEXES ${include_filename_patterns}
            POST_EXCLUDE_REGEXES ".*"
            POST_INCLUDE_REGEXES ${include_filename_patterns}
            DIRECTORIES ${additional_search_directories}
            )
elseif(executable_list)
    file(GET_RUNTIME_DEPENDENCIES
            EXECUTABLES ${executable_list}
            RESOLVED_DEPENDENCIES_VAR resolved_dependencies
            UNRESOLVED_DEPENDENCIES_VAR unresolved_dependencies
            PRE_EXCLUDE_REGEXES ".*"
            PRE_INCLUDE_REGEXES ${include_filename_patterns}
            POST_EXCLUDE_REGEXES ".*"
            POST_INCLUDE_REGEXES ${include_filename_patterns}
            DIRECTORIES ${additional_search_directories}
            )
else()
    file(GET_RUNTIME_DEPENDENCIES
            LIBRARIES ${library_list}
            RESOLVED_DEPENDENCIES_VAR resolved_dependencies
            UNRESOLVED_DEPENDENCIES_VAR unresolved_dependencies
            PRE_EXCLUDE_REGEXES ".*"
            PRE_INCLUDE_REGEXES ${include_filename_patterns}
            POST_EXCLUDE_REGEXES ".*"
            POST_INCLUDE_REGEXES ${include_filename_patterns}
            DIRECTORIES ${additional_search_directories}
            )
endif()

foreach (filepath ${resolved_dependencies})
    get_filename_component(filename ${filepath} NAME)
    string(TOLOWER "${filename}" filename)
    if (NOT filename IN_LIST exclude_filename_list)
        message("Copying resolved dependency ${filepath}")
        file(COPY "${filepath}" DESTINATION "${copy_to_directory}" FOLLOW_SYMLINK_CHAIN)
    endif()
endforeach ()

list(LENGTH unresolved_dependencies unresolved_dependency_count)
if("${unresolved_dependency_count}" GREATER 0)
    message(WARNING "Unresolved dependencies detected: ${unresolved_dependencies}")
endif()

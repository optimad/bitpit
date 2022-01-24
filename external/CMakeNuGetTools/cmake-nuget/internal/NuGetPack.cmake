## Include implementation
include("${CMAKE_CURRENT_LIST_DIR}/NuGetPack.nuspec.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/NuGetPack.nupkg.cmake")

## Internal cache variables
set(NUGET_NUSPEC_INDENT_SIZE "  " CACHE INTERNAL
    "Specifies the size of a single indentation level for generated .nuspec files."
)

## Public interface.
function(nuget_generate_nuspec_files)
    # Sanity checks
    if("${NUGET_COMMAND}" STREQUAL "")
        message(STATUS "NUGET_COMMAND is empty: CMakeNuGetTools is disabled, no .nuspec files are written.")
        return()
    endif()
    if("${ARGV}" STREQUAL "")
        message(FATAL_ERROR "No arguments provided.")
        return()
    endif()
    message(STATUS "Writing .nuspec file(s)...")
    # Begin .nuspec XML file
    set(NUSPEC_CONTENT "<?xml version=\"1.0\" encoding=\"utf-8\"?>")
    string(APPEND NUSPEC_CONTENT "\n<package xmlns=\"http://schemas.microsoft.com/packaging/2011/10/nuspec.xsd\">")
    # Process arguments
    set(options "")
    set(oneValueArgs CMAKE_OUTPUT_DIR)
    set(multiValueArgs CMAKE_CONFIGURATIONS METADATA FILES)
    cmake_parse_arguments(NUARG
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGV}
    )
    nuget_internal_helper_error_if_unparsed_args(
        "${NUARG_UNPARSED_ARGUMENTS}"
        "${NUARG_KEYWORDS_MISSING_VALUES}"
    )
    # METADATA arg: /package/metadata in .nuspec XML file (METADATA as section identifier CMake argument)
    nuget_internal_helper_error_if_empty("${NUARG_METADATA}" "METADATA identifier is not found: it is a required element (/package/metadata) of a .nuspec XML file.")
    nuget_internal_nuspec_process_metadata_args("${NUGET_NUSPEC_INDENT_SIZE}" "${NUSPEC_CONTENT}" NUSPEC_CONTENT PACKAGE_ID ${NUARG_METADATA})
    # FILES arg: /package/files in .nuspec XML file (FILES as section identifier CMake argument)
    nuget_internal_helper_error_if_empty("${NUARG_FILES}"
        "FILES must not be empty: although the files node is not a required element (/package/files) of a .nuspec XML file, "
        "the implementation of the nuget_generate_nuspec_files() CMake command requires you to generate a non-empty files node."
    )
    nuget_internal_nuspec_process_files_args("${NUGET_NUSPEC_INDENT_SIZE}" "${NUSPEC_CONTENT}" NUSPEC_CONTENT ${NUARG_FILES})
    # End .nuspec XML file
    string(APPEND NUSPEC_CONTENT "\n</package>")
    # Write output file(s)
    # Top-level CMAKE_* args: CMake-specific (without special section identifier CMake argument)
    nuget_internal_nuspec_generate_output("${NUSPEC_CONTENT}" "${PACKAGE_ID}" "${NUARG_CMAKE_OUTPUT_DIR}" "${NUARG_CMAKE_CONFIGURATIONS}")
endfunction()

## Public interface.
function(nuget_merge_nuspec_files)
    # Inputs
    set(options "")
    set(oneValueArgs OUTPUT)
    set(multiValueArgs INPUTS)
    cmake_parse_arguments(NUARG
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGV}
    )
    nuget_internal_helper_error_if_unparsed_args(
        "${NUARG_UNPARSED_ARGUMENTS}"
        "${NUARG_KEYWORDS_MISSING_VALUES}"
    )
    nuget_internal_helper_error_if_empty("${NUARG_OUTPUT}" "You must provide a filepath as an OUTPUT argument.")
    nuget_internal_helper_error_if_empty("${NUARG_INPUTS}" "You must provide at least one filepath to a .nuspec file as an INPUTS argument.")
    # Actual functionality
    nuget_internal_merge_n_nuspec_files("${NUARG_OUTPUT}" ${NUARG_INPUTS})
endfunction()

## Public interface.
function(nuget_pack)
    # Sanity checks
    if("${NUGET_COMMAND}" STREQUAL "")
        message(STATUS "NUGET_COMMAND is empty: CMakeNuGetTools is disabled, no packages are packed.")
        return()
    endif()
    # Inputs
    set(options "")
    set(oneValueArgs NUSPEC_FILEPATH OUTPUT_DIRECTORY VERSION_OVERRIDE)
    set(multiValueArgs "")
    cmake_parse_arguments(NUARG
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGV}
    )
    nuget_internal_helper_error_if_unparsed_args(
        "${NUARG_UNPARSED_ARGUMENTS}"
        "${NUARG_KEYWORDS_MISSING_VALUES}"
    )
    nuget_internal_helper_error_if_empty("${NUARG_NUSPEC_FILEPATH}" "You must provide a filepath to a .nuspec file as a NUSPEC_FILEPATH argument.")
    nuget_internal_helper_error_if_empty("${NUARG_OUTPUT_DIRECTORY}" "You must provide an output directory where the .nupkg is created as an OUTPUT_DIRECTORY argument.")
    # Actual functionality
    nuget_internal_pack("${NUARG_NUSPEC_FILEPATH}" "${NUARG_OUTPUT_DIRECTORY}" "${NUARG_VERSION_OVERRIDE}")
endfunction()

## Public interface.
function(nuget_pack_install)
    # Sanity checks
    if("${NUGET_COMMAND}" STREQUAL "")
        message(STATUS "NUGET_COMMAND is empty: CMakeNuGetTools is disabled, no packages are packed.")
        return()
    endif()
    # Inputs
    set(options "")
    set(oneValueArgs PACKAGE VERSION OUTPUT_DIRECTORY SOURCE)
    set(multiValueArgs "")
    cmake_parse_arguments(NUARG
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGV}
    )
    nuget_internal_helper_error_if_unparsed_args(
        "${NUARG_UNPARSED_ARGUMENTS}"
        "${NUARG_KEYWORDS_MISSING_VALUES}"
    )
    nuget_internal_helper_error_if_empty("${NUARG_PACKAGE}" "PACKAGE_ID must be non-empty.")
    nuget_internal_helper_error_if_empty("${NUARG_VERSION}" "PACKAGE_VERSION must be non-empty.")
    nuget_internal_helper_error_if_empty("${NUARG_OUTPUT_DIRECTORY}" "OUTPUT_DIRECTORY must be non-empty.")
    nuget_internal_helper_error_if_empty("${NUARG_SOURCE}" "SOURCE must be non-empty.")
    # Actual functionality
    nuget_internal_pack_install("${NUARG_PACKAGE}" "${NUARG_VERSION}" "${NUARG_OUTPUT_DIRECTORY}" "${NUARG_SOURCE}")
endfunction()

## Public interface.
function(nuget_pack_push)
    # Sanity checks
    if("${NUGET_COMMAND}" STREQUAL "")
        message(STATUS "NUGET_COMMAND is empty: CMakeNuGetTools is disabled, no packages are packed.")
        return()
    endif()
    # Inputs
    set(options "")
    set(oneValueArgs PACKAGE_PATH API_KEY SOURCE)
    set(multiValueArgs "")
    cmake_parse_arguments(NUARG
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGV}
    )
    nuget_internal_helper_error_if_unparsed_args(
        "${NUARG_UNPARSED_ARGUMENTS}"
        "${NUARG_KEYWORDS_MISSING_VALUES}"
    )
    nuget_internal_helper_error_if_empty("${NUARG_PACKAGE_PATH}" "PACKAGE_PATH must be non-empty.")
    nuget_internal_helper_error_if_empty("${NUARG_SOURCE}" "SOURCE must be non-empty.")
    # Actual functionality
    nuget_internal_pack_push("${NUARG_PACKAGE_PATH}" "${NUARG_API_KEY}" "${NUARG_SOURCE}")
endfunction()

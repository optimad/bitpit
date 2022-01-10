# Internal. See arguments below.
function(nuget_internal_single_register_as_interface)
    # Inputs
    set(options INTERFACE)
    set(oneValueArgs PACKAGE VERSION)
    set(multiValueArgs "")
    cmake_parse_arguments(NUARG
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGV}
    )
    nuget_internal_helper_error_if_unparsed_args(
        "${NUARG_UNPARSED_ARGUMENTS}"
        "${NUARG_KEYWORDS_MISSING_VALUES}"
    )
    # Actual functionality
    nuget_internal_core_register("${NUARG_PACKAGE}" "${NUARG_VERSION}" INTERFACE)
endfunction()

# Internal. See arguments below. This should be called if IMPORT_DOT_TARGETS or IMPORT_DOT_TARGETS_AS
# is explicitly specified by the user. Default usage requirement is PRIVATE.
function(nuget_internal_single_import_dot_targets)
    # Inputs
    set(options PUBLIC PRIVATE IMPORT_DOT_TARGETS)
    set(oneValueArgs PACKAGE VERSION IMPORT_DOT_TARGETS_AS)
    set(multiValueArgs INCLUDE_DIRS)
    cmake_parse_arguments(NUARG
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGV}
    )
    nuget_internal_helper_error_if_unparsed_args(
        "${NUARG_UNPARSED_ARGUMENTS}"
        "${NUARG_KEYWORDS_MISSING_VALUES}"
    )
    # Stricter validation
    if(NUARG_PUBLIC)
        if(NUARG_PRIVATE)
            message(FATAL_ERROR "PUBLIC and PRIVATE usage requirement options are mutually exclusive.")
        endif()
        set(USAGE_REQUIREMENT PUBLIC)
    else()
        set(USAGE_REQUIREMENT PRIVATE) # Default
    endif()
    # Actual functionality
    nuget_internal_core_register("${NUARG_PACKAGE}" "${NUARG_VERSION}" "${USAGE_REQUIREMENT}")
    nuget_internal_core_install("${NUARG_PACKAGE}" "${NUARG_VERSION}")
    nuget_internal_core_import_dot_targets(
        "${NUARG_PACKAGE}"
        "${NUARG_VERSION}"
        "${NUARG_IMPORT_DOT_TARGETS_AS}"
        "${NUARG_INCLUDE_DIRS}"
    )
endfunction()

# Internal. See arguments below. This should be called by default if no import method is specified explicitly.
# NOTE: the IMPORT_CMAKE_EXPORTS option is only cosmetics for the user. The presence of CMAKE_PREFIX_PATHS is not treated
# as a differentiator indicating an "IMPORT_CMAKE_EXPORTS" import method either: if the user explicitly provided a
# different import method, the dispatcher logic dispatches the call accordingly and CMAKE_PREFIX_PATHS is not taken into
# account. E.g. providing IMPORT_DOT_TARGETS and CMAKE_PREFIX_PATHS in nuget_internal_single_dependencies() would dispatch the call
# to nuget_internal_single_import_dot_targets() but CMAKE_PREFIX_PATHS does not have any meaning there, so that function should
# raise a CMake Error. Default usage requirement is PRIVATE.
function(nuget_internal_single_import_cmake_exports)
    # Inputs
    set(options PUBLIC PRIVATE IMPORT_CMAKE_EXPORTS CMAKE_APPEND_PATHS)
    set(oneValueArgs PACKAGE VERSION CMAKE_TOOLCHAIN_FILE)
    set(multiValueArgs CMAKE_PREFIX_PATHS CMAKE_MODULE_PATHS)
    cmake_parse_arguments(NUARG
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGV}
    )
    nuget_internal_helper_error_if_unparsed_args(
        "${NUARG_UNPARSED_ARGUMENTS}"
        "${NUARG_KEYWORDS_MISSING_VALUES}"
    )
    # Stricter validation
    if(NUARG_PUBLIC)
        if(NUARG_PRIVATE)
            message(FATAL_ERROR "PUBLIC and PRIVATE usage requirement options are mutually exclusive.")
        endif()
        set(USAGE_REQUIREMENT PUBLIC)
    else()
        set(USAGE_REQUIREMENT PRIVATE) # Default
    endif()
    # Actual functionality
    nuget_internal_core_register("${NUARG_PACKAGE}" "${NUARG_VERSION}" "${USAGE_REQUIREMENT}")
    nuget_internal_core_install("${NUARG_PACKAGE}" "${NUARG_VERSION}")
    nuget_internal_core_import_cmake_exports(
        "${NUARG_PACKAGE}"
        "${NUARG_VERSION}"
        "${NUARG_CMAKE_PREFIX_PATHS}"
        "${NUARG_CMAKE_MODULE_PATHS}"
        "${NUARG_CMAKE_APPEND_PATHS}"
        "${NUARG_CMAKE_TOOLCHAIN_FILE}"
    )
endfunction()

# Internal. Dispatcher to above functions.
function(nuget_internal_single_dependencies)
    # Case: INTERFACE
    list(FIND ARGV INTERFACE INTERFACE_IDX)
    if(NOT ${INTERFACE_IDX} EQUAL -1)
        nuget_internal_single_register_as_interface(${ARGV})
        return()
    endif()
    # Case: IMPORT_DOT_TARGETS or IMPORT_DOT_TARGETS_AS
    list(FIND ARGV IMPORT_DOT_TARGETS IMPORT_DOT_TARGETS_IDX)
    list(FIND ARGV IMPORT_DOT_TARGETS_AS IMPORT_DOT_TARGETS_AS_IDX)
    if(NOT ${IMPORT_DOT_TARGETS_IDX} EQUAL -1 OR NOT ${IMPORT_DOT_TARGETS_AS_IDX} EQUAL -1)
        nuget_internal_single_import_dot_targets(${ARGV})
        return()
    endif()
    # Default: no explicit import method (other than IMPORT_CMAKE_EXPORTS) provided
    nuget_internal_single_import_cmake_exports(${ARGV}) # Default
endfunction()

# Internal. Process each PACKAGE argument pack one-by-one. Deliberately a *function*.
function(nuget_internal_foreach_dependencies)
    set(ARGS_HEAD "")
    set(ARGS_TAIL ${ARGV})
    while(NOT "${ARGS_TAIL}" STREQUAL "")
        nuget_internal_helper_cut_arg_list(PACKAGE "${ARGS_TAIL}" ARGS_HEAD ARGS_TAIL)
        nuget_internal_single_dependencies(${ARGS_HEAD})
    endwhile()
endfunction()

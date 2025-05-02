
# Initialize dependency resolution
#
# This function sets up the collection of third-party libraries and Windows runtime DLLs.
#
# \param PREFIX is the prefix used to create unique variable names for dependency resolution
# \param DEFAULT_THIRD_PARTY_PATTERNS is a semicolon separated list that is used to initialize
# the patterns to identify the third-party libraries
# \param DEFAULT_THIRD_PARTY_DIR is the value that will be used to initialize the path of the
# directory where the third-party libraries will be installed. This can be an absolute path or
# a path relative to the CMake install prefix
# \param DEFAULT_RUNTIME_DLL_EXCLUDE_PATTERNS is a semicolon separated list that is used to
# initialize the patterns to exclude certain DLLs from being collected.
# \param DEFAULT_RUNTIME_DLL_DIR is the value that will be used to initialize the path of the
# directory where the runtime DLLs will be installed. This can be an absolute path or a path
# relative to the CMake install prefix
function(initialize_dependency_resolution PREFIX DEFAULT_THIRD_PARTY_PATTERNS DEFAULT_THIRD_PARTY_DIR DEFAULT_RUNTIME_DLL_EXCLUDE_PATTERNS DEFAULT_RUNTIME_DLL_DIR)

    string(TOUPPER ${PREFIX} UPPER_PREFIX)
    set(DEPENDENCY_RESOLUTION_PREFIX ${UPPER_PREFIX} CACHE INTERNAL "Prefix prepended to the names of the variables related to dependency resolution")

    option(${DEPENDENCY_RESOLUTION_PREFIX}_ENABLE_SELF_CONTAINMENT "If set, a curated list of third-party libraries will be copied in the installation \
    directory and the library will use their symbols in preference to global symbols with the same name." OFF)
    set(${DEPENDENCY_RESOLUTION_PREFIX}_SELF_CONTAINMENT_THIRD_PARTY_PATTERNS "${DEFAULT_THIRD_PARTY_PATTERNS}" CACHE STRING "A semicolon-separated \
    list of patterns that will be used to identify the libraries needed for self-containment")
    mark_as_advanced(${DEPENDENCY_RESOLUTION_PREFIX}_SELF_CONTAINMENT_THIRD_PARTY_PATTERNS)
    set(${DEPENDENCY_RESOLUTION_PREFIX}_SELF_CONTAINMENT_THIRD_PARTY_DIR "${DEFAULT_THIRD_PARTY_DIR}" CACHE STRING "The path of the directory where \
    third-party libraries needed for self-containment will be installed, it can be an absolute path or a path relative to the CMake install prefix")
    mark_as_advanced(${DEPENDENCY_RESOLUTION_PREFIX}_SELF_CONTAINMENT_THIRD_PARTY_DIR)
    if (${DEPENDENCY_RESOLUTION_PREFIX}_ENABLE_SELF_CONTAINMENT)
        message(STATUS "Self containment is enabled")
        message(STATUS "Third-party libraries will be selected using the following patterns \"${${DEPENDENCY_RESOLUTION_PREFIX}_SELF_CONTAINMENT_THIRD_PARTY_PATTERNS}\"")
        message(STATUS "Third-party libraries needed for self containment will be installed into \"${${DEPENDENCY_RESOLUTION_PREFIX}_SELF_CONTAINMENT_THIRD_PARTY_DIR}\"")
    else()
        message(STATUS "Self containment is disabled")
    endif()

    if (${WIN32})
        set(${DEPENDENCY_RESOLUTION_PREFIX}_WINDOWS_RUNTIME_DLLS_EXCLUDE_PATTERNS "${DEFAULT_RUNTIME_DLL_EXCLUDE_PATTERNS}" CACHE STRING "A semicolon-separated \
        list of patterns used to exclude runtime DLLs from being copied")
        mark_as_advanced(${DEPENDENCY_RESOLUTION_PREFIX}_WINDOWS_RUNTIME_DLLS_EXCLUDE_PATTERNS)
        message(STATUS  "Windows runtime DLLs that match the following patterns \"${${DEPENDENCY_RESOLUTION_PREFIX}_WINDOWS_RUNTIME_DLLS_EXCLUDE_PATTERNS}\" \
        will not be copied")
        set(${DEPENDENCY_RESOLUTION_PREFIX}_WINDOWS_RUNTIME_DLLS_DIR "${DEFAULT_RUNTIME_DLL_DIR}" CACHE STRING "The path of the directory \
        where runtime DLLs will be installed, it can be an absolute path or a path relative to the CMake install prefix")
        mark_as_advanced(${DEPENDENCY_RESOLUTION_PREFIX}_WINDOWS_RUNTIME_DLLS_DIR)
        message(STATUS  "Windows runtime DLLs will be installed into \"${${DEPENDENCY_RESOLUTION_PREFIX}_WINDOWS_RUNTIME_DLLS_DIR}\"")
    endif()

endfunction()

# Configure dependency resolution for the specified target
#
# \param TARGET_NAME is the name of the target for which dependency resolution will be configured
# \param TARGET_INSTALL_DIR is the directory where the target will be installed, the directory
# can be an absolute path or a path relative to the CMake install prefix
function(configure_dependency_resolution TARGET_NAME TARGET_INSTALL_DIR)

    # Collect Windows runtime DLLs
    if (${WIN32})
        include(WindowsRuntime)
        windows_runtime_collect_dlls(${TARGET_NAME}
                                     "${${DEPENDENCY_RESOLUTION_PREFIX}_WINDOWS_RUNTIME_DLLS_EXCLUDE_PATTERNS}"
                                     "${${DEPENDENCY_RESOLUTION_PREFIX}_WINDOWS_RUNTIME_DLLS_DIR}"
                                     "INSTALL")
    endif()

    # Configure dependency search path
    if (${DEPENDENCY_RESOLUTION_PREFIX}_ENABLE_SELF_CONTAINMENT)
        include(SelfContainment)
        configure_self_containment(${TARGET_NAME} "${TARGET_INSTALL_DIR}"
                                   "${${DEPENDENCY_RESOLUTION_PREFIX}_SELF_CONTAINMENT_THIRD_PARTY_PATTERNS}"
                                   "${${DEPENDENCY_RESOLUTION_PREFIX}_SELF_CONTAINMENT_THIRD_PARTY_DIR}")
    else()
       # Append to the runtime search path (rpath) of installed binaries any directories outside the
       # project that are in the linker search path or contain linked library files.
       set_target_properties(${TARGET_NAME} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
    endif()

endfunction()

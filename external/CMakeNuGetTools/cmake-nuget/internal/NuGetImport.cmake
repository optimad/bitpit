## Include implementation
include("${CMAKE_CURRENT_LIST_DIR}/NuGetImport.core.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/NuGetImport.single.cmake")

## Public interface. Needs to be called once before every other nuget_*() command. Otherwise the
## result of those commands are all considered undefined, e.g. we might mistakenly detect the
## same nuget package registered twice but with different versions if you update the version of
## the given package in a nuget_add_dependencies() call.
function(nuget_initialize)
    unset("NUGET_GLOBAL_PACKAGES_DIR" CACHE)
    # Check whether the effective packages directory has been changed
    nuget_internal_helper_get_packages_dir(EFFECTIVE_PACKAGES_DIR)
    if(NUGET_PREV_EFFECTIVE_PACKAGES_DIR AND NOT "${NUGET_PREV_EFFECTIVE_PACKAGES_DIR}" STREQUAL "${EFFECTIVE_PACKAGES_DIR}")
        message(STATUS "Effective packages directory has been changed to \"${EFFECTIVE_PACKAGES_DIR}\" from \"${NUGET_PREV_EFFECTIVE_PACKAGES_DIR}\".")
        set(NUGET_EFFECTIVE_PACKAGES_DIR_CHANGED "TRUE" CACHE INTERNAL "")
    else()
        set(NUGET_EFFECTIVE_PACKAGES_DIR_CHANGED "FALSE" CACHE INTERNAL "")
    endif()
    # Unset cache variables we need to use as global variables within a single CMake run only
    nuget_internal_helper_get_internal_cache_variables_with_prefix(NUGET_DEPENDENCY_ NUGET_DEPENDENCY_VARIABLES)
    nuget_internal_helper_get_internal_cache_variables_with_prefix(NUGET_LAST_DEPENDENCY_ NUGET_LAST_DEPENDENCY_VARIABLES)
    foreach(DEPENDENCY_VAR IN LISTS NUGET_DEPENDENCY_VARIABLES NUGET_LAST_DEPENDENCY_VARIABLES)
        unset("${DEPENDENCY_VAR}" CACHE)
    endforeach()
    # Save current effective packages directory for recognizing whether it has been changed
    set(NUGET_PREV_EFFECTIVE_PACKAGES_DIR "${EFFECTIVE_PACKAGES_DIR}" CACHE INTERNAL "")
    # Helps in giving an error if you never-ever called nuget_initialize() before other public nuget_*() calls
    set(NUGET_IS_INITIALIZED_ONCE_CACHED TRUE CACHE INTERNAL "")
endfunction()

## Public interface. Needs to be macro for properly setting CMAKE_MODULE_PATH
## and CMAKE_PREFIX_PATH. It is assumed to be called from directory scope
## (or from another macro that is in dir. scope etc.).
macro(nuget_add_dependencies)
    # Sanity checks
    if("${NUGET_COMMAND}" STREQUAL "")
        message(STATUS "NUGET_COMMAND is empty: CMakeNuGetTools is disabled, no packages are restored.")
        return()
    endif()
    if(NOT NUGET_IS_INITIALIZED_ONCE_CACHED)
        message(FATAL_ERROR
            "NuGetTools for CMake has never been initialized before. "
            "Please call nuget_initialize() once before any other nuget_*() calls."
        )
    endif()
    if("${ARGV}" STREQUAL "")
        message(FATAL_ERROR "No arguments provided.")
        return()
    endif()
    message(STATUS "Importing NuGet package dependencies...")
    # Reset last registered packages list. This is about to be filled in with
    # packages registered via only this single nuget_add_dependencies() call.
    set(NUGET_LAST_DEPENDENCIES_REGISTERED "" CACHE INTERNAL "")
    # Process each PACKAGE argument pack one-by-one. This is a *function* call.
    nuget_internal_foreach_dependencies(${ARGV})
    # Foreach's loop_var should not introduce a new real variable: we are safe macro-wise.
    foreach(PACKAGE_ID IN LISTS NUGET_LAST_DEPENDENCIES_REGISTERED)
        # Set CMAKE_MODULE_PATH and CMAKE_PREFIX_PATH via a *macro* call. Since
        # nuget_add_dependencies() is a macro as well, no new scopes are introduced
        # between the call of nuget_add_dependencies() and setting those variables.
        # I.e. CMake's find_package() will respect those set variables within the
        # same scope (or below directory scopes for example).
        nuget_internal_core_import_cmake_exports_set_cmake_paths("${PACKAGE_ID}")
    endforeach()
    # NOTE: Make sure we did not introduce new normal variables here. Then we are safe macro-wise.
    # (NUGET_LAST_DEPENDENCIES_REGISTERED is an internal *cache* variable so that does not count.)
endmacro()

## Public interface. Returns the list of NuGet package IDs of registered dependencies.
function(nuget_get_dependencies OUT_DEPENDENCIES)
    nuget_internal_helper_get_internal_cache_variables_with_prefix(NUGET_DEPENDENCY_VERSION_ NUGET_PREFIXED_DEPENDENCIES)
    string(REGEX REPLACE "(^|;)NUGET_DEPENDENCY_VERSION_" "\\1" NUGET_DEPENDENCIES "${NUGET_PREFIXED_DEPENDENCIES}")
    set("${OUT_DEPENDENCIES}" "${NUGET_DEPENDENCIES}" PARENT_SCOPE)
endfunction()

## Public interface. Returns the list of NuGet package IDs of installed dependencies. Please
## note that packages with INTERFACE usage requirement are never installed. Please also note
## that a package is known to be installed if nuget_is_dependency_installed() returns TRUE
## for the given package.
function(nuget_get_installed_dependencies OUT_PACKAGE_IDS)
    set(PACKAGE_IDS "")
    nuget_get_dependencies(DEPENDENCIES)
    foreach(DEPENDENCY IN LISTS DEPENDENCIES)
        nuget_is_dependency_installed("${DEPENDENCY}" IS_INSTALLED)
        if(NOT IS_INSTALLED)
            continue()
        endif()
        list(APPEND PACKAGE_IDS "${DEPENDENCY}")
    endforeach()
    set("${OUT_PACKAGE_IDS}" "${PACKAGE_IDS}" PARENT_SCOPE)
endfunction()

## Public interface. Returns the version of the given registered NuGet package ID.
function(nuget_get_dependency_version PACKAGE_ID OUT_VERSION)
    if(NOT DEFINED "NUGET_DEPENDENCY_VERSION_${PACKAGE_ID}")
        message(FATAL_ERROR "\"${PACKAGE_ID}\" is not registered.")
    endif()
    set("${OUT_VERSION}" "${NUGET_DEPENDENCY_VERSION_${PACKAGE_ID}}" PARENT_SCOPE)
endfunction()

## Public interface. Returns the usage requirement of the given registered NuGet package ID.
function(nuget_get_dependency_usage PACKAGE_ID OUT_USAGE)
    if(NOT DEFINED "NUGET_DEPENDENCY_VERSION_${PACKAGE_ID}")
        message(FATAL_ERROR "\"${PACKAGE_ID}\" is not registered.")
    endif()
    set("${OUT_USAGE}" "${NUGET_DEPENDENCY_USAGE_${PACKAGE_ID}}" PARENT_SCOPE)
endfunction()

## Public interface. Returns the absolute root directory path to the given installed NuGet package ID.
## Please note that packages with INTERFACE usage requirement are never installed. Only returns a non-
## empty path if nuget_is_dependency_installed() returns TRUE.
function(nuget_get_dependency_dir PACKAGE_ID OUT_DIR)
    if(NOT DEFINED "NUGET_DEPENDENCY_VERSION_${PACKAGE_ID}")
        message(FATAL_ERROR "\"${PACKAGE_ID}\" is not registered.")
    endif()
    set("${OUT_DIR}" "${NUGET_DEPENDENCY_DIR_${PACKAGE_ID}}" PARENT_SCOPE)
endfunction()

## Public interface. Returns the absolute root directory paths to all installed NuGet packages.
## Please note that packages with INTERFACE usage requirement are never installed. Please also
## note that a package is known to be installed if nuget_is_dependency_installed() returns TRUE
## for the given package.
function(nuget_get_installed_dependencies_dirs OUT_DIRS)
    set(DIRS "")
    nuget_get_dependencies(DEPENDENCIES)
    foreach(DEPENDENCY IN LISTS DEPENDENCIES)
        nuget_get_dependency_dir("${DEPENDENCY}" PACKAGE_DIR)
        nuget_is_dependency_installed("${DEPENDENCY}" IS_INSTALLED)
        if(NOT IS_INSTALLED)
            continue()
        endif()
        list(APPEND DIRS "${PACKAGE_DIR}")
    endforeach()
    set("${OUT_DIRS}" "${DIRS}" PARENT_SCOPE)
endfunction()

## Public interface. Returns whether the given registered NuGet package ID is or is not installed.
## Please note that packages with INTERFACE usage requirement are never installed.
function(nuget_is_dependency_installed PACKAGE_ID OUT_IS_INSTALLED)
    if(NOT DEFINED "NUGET_DEPENDENCY_VERSION_${PACKAGE_ID}")
        message(FATAL_ERROR "\"${PACKAGE_ID}\" is not registered.")
    endif()
    set("${OUT_IS_INSTALLED}" "${NUGET_DEPENDENCY_INSTALLED_${PACKAGE_ID}}" PARENT_SCOPE)
endfunction()

# Internal. Registers the given package, nothing more. The nuget_internal_core_install()
# function should be called for actually installing the package. USAGE_REQUIREMENT
# is expected to be either PUBLIC, INTERFACE, or PRIVATE.
function(nuget_internal_core_register
    PACKAGE_ID
    PACKAGE_VERSION
    USAGE_REQUIREMENT
)
    # Inputs
    nuget_internal_helper_error_if_empty("${PACKAGE_ID}" "Package ID must be provided.")
    nuget_internal_helper_error_if_empty("${PACKAGE_VERSION}" "Package version must be provided.")
    nuget_internal_helper_error_if_empty("${USAGE_REQUIREMENT}" "Usage requirement for package must be provided.")
    if(DEFINED "NUGET_DEPENDENCY_VERSION_${PACKAGE_ID}")
        # Safety check. Do not allow more than one package to be used within the same CMake build
        # system with the same package ID but with different version numbers. The same version
        # can be installed as many times as you want -- only the first successful NuGet install
        # invocation has the effect of fetching the package to the OutputDirectory. Additional
        # calls to NuGet install are allowed but they would not do anything as the package is
        # already installed.
        #
        # NOTE. This check fails to catch additional dependencies installed by NuGet in
        # execute_process() in nuget_internal_core_install(). This problem can be avoided if
        # nuget_internal_core_install() is invoked for all the packages that are your transitive
        # dependencies.
        if(NOT "${NUGET_DEPENDENCY_VERSION_${PACKAGE_ID}}" STREQUAL "${PACKAGE_VERSION}")
            message(FATAL_ERROR
                "NuGet package \"${PACKAGE_ID}\" is already registered "
                "with version ${NUGET_DEPENDENCY_VERSION_${PACKAGE_ID}}. "
                "You are trying to register version ${PACKAGE_VERSION}."
            )
        endif()
        # Registering the same package multiple times with different usage requirements
        # is not allowed.
        if(NOT "${NUGET_DEPENDENCY_USAGE_${PACKAGE_ID}}" STREQUAL "${USAGE_REQUIREMENT}")
            message(FATAL_ERROR
                "NuGet package \"${PACKAGE_ID}\" is already registered "
                "with usage requirement ${NUGET_DEPENDENCY_USAGE_${PACKAGE_ID}}. "
                "You are trying to register with usage ${USAGE_REQUIREMENT}."
            )
        endif()
    endif()
    # Set internal cache variables
    set("NUGET_DEPENDENCY_VERSION_${PACKAGE_ID}" ${PACKAGE_VERSION} CACHE INTERNAL
        "The version of the registered package \"${PACKAGE_ID}\"."
    )
    set("NUGET_DEPENDENCY_USAGE_${PACKAGE_ID}" ${USAGE_REQUIREMENT} CACHE INTERNAL
        "The usage requirement of the registered package \"${PACKAGE_ID}\"."
    )
    # The same package is not allowed to be registered within the same nuget_add_dependencies() call
    list(FIND NUGET_LAST_DEPENDENCIES_REGISTERED "${PACKAGE_ID}" LAST_DEPENDENCY_IDX)
    if(NOT ${LAST_DEPENDENCY_IDX} EQUAL -1)
        message(FATAL_ERROR
            "The \"${PACKAGE_ID}\" package is already registered by the same nuget_add_dependencies() "
            "call. A PACKAGE argument pack with the same package ID can only occur once in the "
            "argument list of the same nuget_add_dependencies() call."
        )
    endif()
    # Spec. cache var. for listing packages registered by a single nuget_add_dependencies() call
    set(LAST_DEPENDENCIES_REGISTERED ${NUGET_LAST_DEPENDENCIES_REGISTERED})
    list(APPEND LAST_DEPENDENCIES_REGISTERED "${PACKAGE_ID}")
    set(NUGET_LAST_DEPENDENCIES_REGISTERED "${LAST_DEPENDENCIES_REGISTERED}" CACHE INTERNAL
        "List of packages registered by the last nuget_add_dependencies() call."
    )
endfunction()

## Internal. Runs NuGet install with PACKAGE_ID and PACKAGE_VERSION.
## The OutputDirectory is set by the CACHE variable NUGET_PACKAGES_DIR by default.
function(nuget_internal_core_install
    PACKAGE_ID
    PACKAGE_VERSION
)
    # Inputs
    nuget_internal_helper_error_if_empty("${NUGET_COMMAND}"
        "No NuGet executable was provided; this means NuGetTools should have been disabled, and "
        "we should not ever reach a call to nuget_internal_core_install()."
    )
    if(NUGET_NO_CACHE)
        set(NUGET_INSTALL_NO_CACHE_OPTION "-NoCache")
    endif()
    if(NUGET_DIRECT_DOWNLOAD)
        set(NUGET_INSTALL_DIRECT_DOWNLOAD_OPTION "-DirectDownload")
    endif()
    if(NOT "${NUGET_PACKAGE_SAVE_MODE}" STREQUAL "")
        set(NUGET_PACKAGE_SAVE_MODE_OPTION -PackageSaveMode "${NUGET_PACKAGE_SAVE_MODE}")
    endif()
    if(NUGET_EXCLUDE_VERSION)
        set(NUGET_EXCLUDE_VERSION_OPTION "-ExcludeVersion")
    endif()
    # Execute install
    #
    # NOTE. Output is not parsed for additionally installed dependencies. It does not worth
    # the coding effort. Just make sure nuget_internal_core_install() is explicitly called for all
    # your PUBLIC and PRIVATE transitive dependencies. I.e. the user should explicitly list
    # all dependencies that should be installed.
    nuget_internal_helper_get_packages_dir(PACKAGES_DIR)
    execute_process(
        COMMAND "${NUGET_COMMAND}" install ${PACKAGE_ID}
            -Version ${PACKAGE_VERSION}
            -OutputDirectory "${PACKAGES_DIR}"
            -NonInteractive
            ${NUGET_INSTALL_NO_CACHE_OPTION}
            ${NUGET_INSTALL_DIRECT_DOWNLOAD_OPTION}
            ${NUGET_PACKAGE_SAVE_MODE_OPTION}
            ${NUGET_EXCLUDE_VERSION_OPTION}
        ERROR_VARIABLE
            NUGET_INSTALL_ERROR_VAR
        RESULT_VARIABLE
            NUGET_INSTALL_RESULT_VAR
    )
    nuget_internal_helper_error_if_not_empty(
        "${NUGET_INSTALL_ERROR_VAR}"
        "Running NuGet package install encountered some errors: "
    )
    if(NOT ${NUGET_INSTALL_RESULT_VAR} EQUAL 0)
        message(FATAL_ERROR "NuGet package install returned with: \"${NUGET_INSTALL_RESULT_VAR}\"")
    endif()
    # Mark package as succesfully installed
    set("NUGET_DEPENDENCY_INSTALLED_${PACKAGE_ID}" TRUE CACHE INTERNAL
        "True if the package \"${PACKAGE_ID}\" is successfully installed."
    )
    nuget_internal_helper_get_package_dir(${PACKAGE_ID} ${PACKAGE_VERSION} PACKAGE_DIR)
    set("NUGET_DEPENDENCY_DIR_${PACKAGE_ID}" "${PACKAGE_DIR}" CACHE INTERNAL
        "Absolute path to the directory of the installed package \"${PACKAGE_ID}\"."
    )
    # Check whether the effective package directory for this particular package has been changed
    if(NUGET_PREV_DEPENDENCY_DIR_${PACKAGE_ID} AND NOT "${NUGET_PREV_DEPENDENCY_DIR_${PACKAGE_ID}}" STREQUAL "${PACKAGE_DIR}")
        message(STATUS "Effective package directory of package \"${PACKAGE_ID}.${PACKAGE_VERSION}\" has been changed to \"${PACKAGE_DIR}\" from \"${NUGET_PREV_DEPENDENCY_DIR_${PACKAGE_ID}}\".")
        if(NUGET_AUTO_UNSET_CACHE_VARS_WHEN_PACKAGE_DIR_CHANGED)
            nuget_internal_helper_unset_cache_vars_containing("${NUGET_PREV_DEPENDENCY_DIR_${PACKAGE_ID}}" "NUGET_")
        endif()
    endif()
    set("NUGET_PREV_DEPENDENCY_DIR_${PACKAGE_ID}" "${PACKAGE_DIR}" CACHE INTERNAL "")
endfunction()

## Internal. Only Visual Studio generators are compatible. Creates a CMake build target
## named IMPORT_AS (if non-empty, otherwise: PACKAGE_ID) from the PACKAGE_ID package with
## PACKAGE_VERSION version using the "build/native/${PACKAGE_ID}.targets" within the package.
## INCLUDE_DIRS is required for better Visual Studio editing experience.
function(nuget_internal_core_import_dot_targets
    PACKAGE_ID
    PACKAGE_VERSION
    IMPORT_AS
    INCLUDE_DIRS
)
    # Compatibility check
    if(NOT CMAKE_GENERATOR MATCHES "^Visual Studio")
        message(FATAL_ERROR
            "You are trying to import the \"${PACKAGE_ID}\" NuGet package via a "
            ".targets file whereas your CMake generator is \"${CMAKE_GENERATOR}\". "
            "Only Visual Studio generators are compatible with this import method."
        )
    endif()
    # Inputs
    if("${IMPORT_AS}" STREQUAL "")
        set(IMPORT_AS "${PACKAGE_ID}")
    endif()
    nuget_internal_helper_list_transform_prepend(
        "${INCLUDE_DIRS}" "${NUGET_DEPENDENCY_DIR_${PACKAGE_ID}}/" INCLUDE_DIRS
    )
    if(TARGET ${IMPORT_AS})
        message(FATAL_ERROR
            "You are trying to import the \"${PACKAGE_ID}\" NuGet package "
            "with an already existing target name \"${IMPORT_AS}\"."
        )
    endif()
    # NOTE. The fallback "build/${PACKAGE_ID}.targets" is deliberately not tried.
    # One should always place the native .target file under "build/native/".
    set(DOT_TARGETS_FILE
        "${NUGET_DEPENDENCY_DIR_${PACKAGE_ID}}/build/native/${PACKAGE_ID}.targets" # Default (and non-settable)
    )
    if(NOT EXISTS "${DOT_TARGETS_FILE}")
        message(FATAL_ERROR "The file \"${DOT_TARGETS_FILE}\" does not exist.")
    endif()

    # Create build target
    add_library(${IMPORT_AS} INTERFACE IMPORTED GLOBAL)
    # See https://gitlab.kitware.com/cmake/cmake/issues/16340 -- since CMake version 2.8.12 (?)
    set_property(TARGET ${IMPORT_AS} PROPERTY INTERFACE_LINK_LIBRARIES "${DOT_TARGETS_FILE}")
    if(NOT "${INCLUDE_DIRS}" STREQUAL "")
        # Experience shows that the Visual Studio editor does not recognize anything included
        # unless you set this property. Building your target should work regardless of setting
        # this property if the imported .targets file is working properly.
        set_property(TARGET ${IMPORT_AS} PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${INCLUDE_DIRS})
    endif()
endfunction()

## Internal. Preparations for prepending the CMAKE_PREFIX_PATHS and CMAKE_MODULE_PATHS relative-to-package-directory
## paths to the CMAKE_PREFIX_PATH and CMAKE_MODULE_PATHS if CMAKE_APPEND_PATHS is FALSE. If CMAKE_APPEND_PATHS is TRUE,
## the operation becomes an append instead of a prepend later on.
function(nuget_internal_core_import_cmake_exports
    PACKAGE_ID
    PACKAGE_VERSION
    CMAKE_PREFIX_PATHS
    CMAKE_MODULE_PATHS
    CMAKE_APPEND_PATHS
    CMAKE_TOOLCHAIN_FILE
)
    # Inputs
    if("${CMAKE_PREFIX_PATHS}" STREQUAL "" AND "${CMAKE_MODULE_PATHS}" STREQUAL "" AND "${CMAKE_TOOLCHAIN_FILE}" STREQUAL "")
        message(FATAL_ERROR "At least one of CMAKE_PREFIX_PATHS, CMAKE_MODULE_PATHS, or CMAKE_TOOLCHAIN_FILE should be non-empty.")
    endif()
    nuget_internal_helper_list_transform_prepend(
        "${CMAKE_PREFIX_PATHS}" "${NUGET_DEPENDENCY_DIR_${PACKAGE_ID}}/" CMAKE_PREFIX_PATHS
    )
    nuget_internal_helper_list_transform_prepend(
        "${CMAKE_MODULE_PATHS}" "${NUGET_DEPENDENCY_DIR_${PACKAGE_ID}}/" CMAKE_MODULE_PATHS
    )
    nuget_internal_helper_list_transform_prepend(
        "${CMAKE_TOOLCHAIN_FILE}" "${NUGET_DEPENDENCY_DIR_${PACKAGE_ID}}/" CMAKE_TOOLCHAIN_FILE
    )

    # Save settings: we do not actually set CMAKE_PREFIX_PATH or CMAKE_MODULE_PATH here,
    # see the call point of nuget_internal_core_import_cmake_exports_set_cmake_paths() for that.
    # Since we are in a new (function) scope here, setting those variables here would not
    # have any effect.
    set("NUGET_LAST_DEPENDENCY_CMAKE_PREFIX_PATHS_${PACKAGE_ID}" "${CMAKE_PREFIX_PATHS}" CACHE INTERNAL "")
    set("NUGET_LAST_DEPENDENCY_CMAKE_MODULE_PATHS_${PACKAGE_ID}" "${CMAKE_MODULE_PATHS}" CACHE INTERNAL "")
    set("NUGET_LAST_DEPENDENCY_CMAKE_APPEND_PATHS_${PACKAGE_ID}" "${CMAKE_APPEND_PATHS}" CACHE INTERNAL "")
    set("NUGET_LAST_DEPENDENCY_CMAKE_TOOLCHAIN_FILE_${PACKAGE_ID}" "${CMAKE_TOOLCHAIN_FILE}" CACHE INTERNAL "")
endfunction()

## Internal. Needs to be macro for properly setting CMAKE_MODULE_PATH or CMAKE_PREFIX_PATH.
## ASSUME: called from directory scope (or from another macro that is in dir. scope etc.).
## E.g. we want to achieve that CMake's find_package() respects our CMAKE_PREFIX_PATH
## modification here.
macro(nuget_internal_core_import_cmake_exports_set_cmake_paths PACKAGE_ID)
    # Modify prefix or module path
    # See https://cmake.org/cmake/help/latest/command/find_package.html#search-procedure
    if(NOT "${NUGET_LAST_DEPENDENCY_CMAKE_PREFIX_PATHS_${PACKAGE_ID}}" STREQUAL "")
        if("${NUGET_LAST_DEPENDENCY_CMAKE_APPEND_PATHS_${PACKAGE_ID}}")
            list(APPEND CMAKE_PREFIX_PATH "${NUGET_LAST_DEPENDENCY_CMAKE_PREFIX_PATHS_${PACKAGE_ID}}")
        else()
            list(INSERT CMAKE_PREFIX_PATH 0 "${NUGET_LAST_DEPENDENCY_CMAKE_PREFIX_PATHS_${PACKAGE_ID}}")
        endif()
    endif()
    if(NOT "${NUGET_LAST_DEPENDENCY_CMAKE_MODULE_PATHS_${PACKAGE_ID}}" STREQUAL "")
        if("${NUGET_LAST_DEPENDENCY_CMAKE_APPEND_PATHS_${PACKAGE_ID}}")
            list(APPEND CMAKE_MODULE_PATH "${NUGET_LAST_DEPENDENCY_CMAKE_MODULE_PATHS_${PACKAGE_ID}}")
        else()
            list(INSERT CMAKE_MODULE_PATH 0 "${NUGET_LAST_DEPENDENCY_CMAKE_MODULE_PATHS_${PACKAGE_ID}}")
        endif()
    endif()
    if(NOT "${NUGET_LAST_DEPENDENCY_CMAKE_TOOLCHAIN_FILE_${PACKAGE_ID}}" STREQUAL "")
        if(NOT "${NUGET_LAST_DEPENDENCY_CMAKE_TOOLCHAIN_FILE}" STREQUAL "")
            message(FATAL_ERROR "CMAKE_TOOLCHAIN_FILE is already set by another package.")
        endif()
        if(NOT "${CMAKE_TOOLCHAIN_FILE}" STREQUAL "" AND
           NOT "${CMAKE_TOOLCHAIN_FILE}" STREQUAL "${NUGET_LAST_DEPENDENCY_CMAKE_TOOLCHAIN_FILE_${PACKAGE_ID}}"
        )
            message(FATAL_ERROR "CMAKE_TOOLCHAIN_FILE is already set by the user.")
        endif()
        if(DEFINED PROJECT_NAME)
            message(FATAL_ERROR "Value of CMAKE_TOOLCHAIN_FILE is ignored after the very first project() command call.")
        endif()
        set(NUGET_LAST_DEPENDENCY_CMAKE_TOOLCHAIN_FILE "${NUGET_LAST_DEPENDENCY_CMAKE_TOOLCHAIN_FILE_${PACKAGE_ID}}" CACHE INTERNAL "")
        set(CMAKE_TOOLCHAIN_FILE "${NUGET_LAST_DEPENDENCY_CMAKE_TOOLCHAIN_FILE_${PACKAGE_ID}}")
    endif()
    # NOTE: Make sure we did not introduce new normal variables here. Then we are safe macro-wise.
endmacro()

cmake_minimum_required(VERSION 3.8.2 FATAL_ERROR)

## User-settable cache variables
set(NUGET_COMMAND "" CACHE STRING "NuGet executable used for package installs. Empty means NuGetTools is disabled. Deliberately not a FILEPATH cache variable: you can set it simply to \"nuget.exe\" if the executable is within your PATH environment.")
set(NUGET_PACKAGES_DIR "${CMAKE_BINARY_DIR}/nuget-packages" CACHE PATH "Path to the directory used by NuGet to store installed packages.")

## Advanced cache variables
option(NUGET_NO_CACHE "Advanced. Add -NoCache option to NuGet install commands: \"Disable using the machine cache as the first package source.\"" FALSE)
option(NUGET_DIRECT_DOWNLOAD "Advanced. Add -DirectDownload option to NuGet install commands: \"Download directly without populating any caches with metadata or binaries.\"" FALSE)
set(NUGET_PACKAGE_SAVE_MODE "" CACHE STRING "Advanced. Add -PackageSaveMode option with passed value if non-empty to NuGet install commands: \"Specifies types of files to save after package installation: nuspec, nupkg, nuspec;nupkg.\"")
option(NUGET_EXCLUDE_VERSION "Advanced. Add -ExcludeVersion option to NuGet install commands: \"If set, the destination folder will contain only the package name, not the version number.\" This is required if you want to avoid very inconvenient error messages like cannot include CMAKE_TOOLCHAIN_FILE or having invalid find_package() cache entries after changing the version of a package in a nuget_add_dependencies() call. Also see https://github.com/katusk/CMakeNuGetTools/issues/2" TRUE)
option(NUGET_GLOBAL_PACKAGES_AS_INSTALL_OUTPUT_DIR "Advanced. Use machine-wide global packages directory as output directory for NuGet install commands instead of NUGET_PACKAGES_DIR." FALSE)
set(NUGET_GLOBAL_PACKAGES_DIR_SUFFIX ".installed" CACHE PATH "Advanced. Used if NUGET_GLOBAL_PACKAGES_AS_INSTALL_OUTPUT_DIR is set to TRUE: given relative path is appended to the path of the global packages directory for package installs.")
option(NUGET_AUTO_UNSET_CACHE_VARS_WHEN_PACKAGE_DIR_CHANGED "Advanced. Automatically unset CMake cache variables that contain the previous packages directory or a package version has been changed and we use package directories with the version number included." FALSE)
mark_as_advanced(FORCE NUGET_NO_CACHE)
mark_as_advanced(FORCE NUGET_DIRECT_DOWNLOAD)
mark_as_advanced(FORCE NUGET_PACKAGE_SAVE_MODE)
mark_as_advanced(FORCE NUGET_EXCLUDE_VERSION)
mark_as_advanced(FORCE NUGET_GLOBAL_PACKAGES_AS_INSTALL_OUTPUT_DIR)
mark_as_advanced(FORCE NUGET_GLOBAL_PACKAGES_DIR_SUFFIX)
mark_as_advanced(FORCE NUGET_AUTO_UNSET_CACHE_VARS_WHEN_PACKAGE_DIR_CHANGED)
set(NUGET_IGNORE_ME ${NUGET_AUTO_UNSET_CACHE_VARS_WHEN_PACKAGE_DIR_CHANGED})

## Includes
include("${CMAKE_CURRENT_LIST_DIR}/internal/NuGetTools.helper.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/internal/NuGetImport.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/internal/NuGetPack.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/internal/NuGetGit.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/internal/NuGetSemVer.cmake")

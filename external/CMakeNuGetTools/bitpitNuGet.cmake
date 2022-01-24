# Specify the version being used as well as the language
cmake_minimum_required(VERSION 3.10)

#store the current directory where the bitpitNuGet.cmake is 
set(BITPIT_NUGET_CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_DIR}) 

# macro to retrieve some bitpit packages directly from nuget (VisualStudio MSVC context)
macro(bitpitNuGet)
    # download directly libxml2, rapidjson, mkl blas, mkl lapack and mkl lapacke 
    # through nuget system and attach it to current VS project
    include("${BITPIT_NUGET_CURRENT_LIST_DIR}/cmake-nuget/NuGetTools.cmake")

    #find the nuget.exe executable 
    find_program(NUGET_COMMAND nuget)
    if(NOT NUGET_COMMAND)
        message(FATAL "Cannot find nuget command line tool that should come with VS!")
    endif()

    # Call this once before any other nuget_* calls.
    nuget_initialize()
    # NuGet will install .
    nuget_add_dependencies(
        PACKAGE libxml2 VERSION 2.7.8.7 CMAKE_PREFIX_PATHS build/native;build/native/lib/v110/x64/Release/dynamic/cdecl
        PACKAGE libiconv VERSION 1.14.0.11 CMAKE_PREFIX_PATHS build/native
        PACKAGE zlib_native.redist VERSION 1.2.11 CMAKE_PREFIX_PATHS build/native
        PACKAGE tencent.rapidjson VERSION 1.1.1 CMAKE_PREFIX_PATHS lib/native
        PACKAGE intelmkl.devel.win-x64 VERSION 2022.0.0.115 CMAKE_PREFIX_PATHS lib/native;lib/native/win-x64
    )

    ## add C library missing GNU libiconv header
    include_directories("${NUGET_PACKAGES_DIR}/libiconv/build/native/include")

    ## set temporary exception to include intelmkl include dirs 
    set(INTELMKL_INCLUDE_DIRS "${NUGET_PACKAGES_DIR}/intelmkl.devel.win-x64/lib/native/include")
 
    #create ad-hoc CMakeConfig for tencent.rapidjson (they do not provide it in nuget package)
    set(CONFIGCMAKEFILENAME "${NUGET_PACKAGES_DIR}/tencent.rapidjson/lib/native/RapidJSONConfig.cmake")
    file(WRITE ${CONFIGCMAKEFILENAME} "get_filename_component(RAPIDJSON_CMAKE_DIR \"\${CMAKE_CURRENT_LIST_FILE}\" PATH)\n")
    file(APPEND ${CONFIGCMAKEFILENAME} "set(RAPIDJSON_INCLUDE_DIRS \"\${RAPIDJSON_CMAKE_DIR}/include\")\n")
    file(APPEND ${CONFIGCMAKEFILENAME} "message(STATUS \"RapidJSON found. Headers: \${RAPIDJSON_INCLUDE_DIRS}\")")
    
    set(CONFIGCMAKEFILEVERSIONNAME "${NUGET_PACKAGES_DIR}/tencent.rapidjson/lib/native/RapidJSONConfigVersion.cmake")
    file(WRITE ${CONFIGCMAKEFILEVERSIONNAME} "SET(PACKAGE_VERSION \"1.1.0\")\n\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "IF (PACKAGE_FIND_VERSION VERSION_EQUAL PACKAGE_VERSION)\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "    SET(PACKAGE_VERSION_EXACT \"true\")\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "ENDIF (PACKAGE_FIND_VERSION VERSION_EQUAL PACKAGE_VERSION)\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "IF (NOT PACKAGE_FIND_VERSION VERSION_GREATER PACKAGE_VERSION)\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "    SET(PACKAGE_VERSION_COMPATIBLE \"true\")\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "ELSE (NOT PACKAGE_FIND_VERSION VERSION_GREATER PACKAGE_VERSION)\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "    SET(PACKAGE_VERSION_UNSUITABLE \"true\")\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "ENDIF (NOT PACKAGE_FIND_VERSION VERSION_GREATER PACKAGE_VERSION)")

    list(APPEND CMAKE_MODULE_PATH "${NUGET_PACKAGES_DIR}/tencent.rapidjson/lib/native/")
endmacro()

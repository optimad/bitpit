## Internal.
function(nuget_internal_pack NUSPEC_FILEPATH OUTPUT_DIRECTORY VERSION_OVERRIDE)
    # Inputs
    nuget_internal_helper_error_if_empty("${NUGET_COMMAND}"
        "No NuGet executable was provided; this means NuGetTools should have been disabled, and "
        "we should not ever reach a call to nuget_internal_pack()."
    )
    if(NOT "${VERSION_OVERRIDE}" STREQUAL "")
        set(VERSION_OVERRIDE "-Version ${VERSION_OVERRIDE}")
    endif()
    # Execute pack
    execute_process(
        COMMAND "${NUGET_COMMAND}" pack "${NUSPEC_FILEPATH}"
            -OutputDirectory "${OUTPUT_DIRECTORY}"
            ${VERSION_OVERRIDE}
            -NonInteractive
            -NoDefaultExcludes
            # Consider -NoPackageAnalysis as option (default FALSE)
        ERROR_VARIABLE
            NUGET_PACK_ERROR_VAR
        RESULT_VARIABLE
            NUGET_PACK_RESULT_VAR
    )
    nuget_internal_helper_error_if_not_empty(
        "${NUGET_PACK_ERROR_VAR}"
        "Running NuGet pack based on \"${NUSPEC_FILEPATH}\" encountered some errors: "
    )
    if(NOT ${NUGET_PACK_RESULT_VAR} EQUAL 0)
        message(FATAL_ERROR "NuGet pack returned with: \"${NUGET_PACK_RESULT_VAR}\"")
    endif()
endfunction()

## Internal.
function(nuget_internal_pack_install
    PACKAGE_ID
    PACKAGE_VERSION
    OUTPUT_DIRECTORY
    SOURCE
)
    # Inputs
    nuget_internal_helper_error_if_empty("${NUGET_COMMAND}"
        "No NuGet executable was provided; this means NuGetTools should have been disabled, and "
        "we should not ever reach a call to nuget_internal_pack_install()."
    )
    # Execute
    execute_process(
        COMMAND "${NUGET_COMMAND}" install ${PACKAGE_ID}
            -Version ${PACKAGE_VERSION}
            -OutputDirectory "${OUTPUT_DIRECTORY}"
            -Source "${SOURCE}"
            -NonInteractive
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
endfunction()

## Internal.
function(nuget_internal_pack_push
    PACKAGE_PATH
    API_KEY
    SOURCE
)
    # Inputs
    nuget_internal_helper_error_if_empty("${NUGET_COMMAND}"
        "No NuGet executable was provided; this means NuGetTools should have been disabled, and "
        "we should not ever reach a call to nuget_internal_pack_push()."
    )
    if(NOT "${API_KEY}" STREQUAL "")
        set(API_KEY -ApiKey "${API_KEY}")
    endif()
    # Execute
    execute_process(
        COMMAND "${NUGET_COMMAND}" push "${PACKAGE_PATH}"
            ${API_KEY}
            -Source "${SOURCE}"
            -NonInteractive
        ERROR_VARIABLE
            NUGET_INSTALL_ERROR_VAR
        RESULT_VARIABLE
            NUGET_INSTALL_RESULT_VAR
    )
    nuget_internal_helper_error_if_not_empty(
        "${NUGET_INSTALL_ERROR_VAR}"
        "Running NuGet package push encountered some errors: "
    )
    if(NOT ${NUGET_INSTALL_RESULT_VAR} EQUAL 0)
        message(FATAL_ERROR "NuGet package push returned with: \"${NUGET_INSTALL_RESULT_VAR}\"")
    endif()
endfunction()

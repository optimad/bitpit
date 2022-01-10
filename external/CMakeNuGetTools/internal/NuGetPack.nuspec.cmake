## Internal. Section: /package/metadata in .nuspec XML file (METADATA as section identifier CMake argument).
function(nuget_internal_nuspec_process_metadata_args NUSPEC_INDENT_SIZE NUSPEC_CONTENT OUT_NUSPEC_CONTENT OUT_PACKAGE_ID)
    # Inputs
    set(options "")
    set(oneValueArgs PACKAGE VERSION DESCRIPTION PROJECT_URL ICON COPYRIGHT
        REPOSITORY_TYPE REPOSITORY_URL REPOSITORY_BRANCH REPOSITORY_COMMIT
    )
    set(multiValueArgs AUTHORS)
    cmake_parse_arguments(NUARG
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN}
    )
    nuget_internal_helper_error_if_unparsed_args(
        "${NUARG_UNPARSED_ARGUMENTS}"
        "${NUARG_KEYWORDS_MISSING_VALUES}"
    )
    # See https://docs.microsoft.com/en-us/nuget/reference/nuspec#general-form-and-schema for below requirements
    nuget_internal_helper_error_if_empty("${NUARG_PACKAGE}" "PACKAGE must not be empty: it is a required element (/package/metadata/id) of a .nuspec XML file.")
    nuget_internal_helper_error_if_empty("${NUARG_VERSION}" "VERSION must not be empty: it is a required element (/package/metadata/version) of a .nuspec XML file.")
    nuget_internal_helper_error_if_empty("${NUARG_DESCRIPTION}" "DESCRIPTION must not be empty: it is a required element (/package/metadata/description) of a .nuspec XML file.")
    nuget_internal_helper_error_if_empty("${NUARG_AUTHORS}" "AUTHORS must not be empty: it is a required element (/package/metadata/authors) of a .nuspec XML file.")
    # Actual functionality
    # Begin /package/metadata
    string(APPEND NUSPEC_CONTENT "\n${NUSPEC_INDENT_SIZE}<metadata>")
    set(NUSPEC_SUBELEMENT_INDENT_SIZE "${NUSPEC_INDENT_SIZE}${NUGET_NUSPEC_INDENT_SIZE}")
    # Required metadata subelements
    string(APPEND NUSPEC_CONTENT "\n${NUSPEC_SUBELEMENT_INDENT_SIZE}<id>${NUARG_PACKAGE}</id>")
    string(APPEND NUSPEC_CONTENT "\n${NUSPEC_SUBELEMENT_INDENT_SIZE}<version>${NUARG_VERSION}</version>")
    string(APPEND NUSPEC_CONTENT "\n${NUSPEC_SUBELEMENT_INDENT_SIZE}<description>${NUARG_DESCRIPTION}</description>")
    string(REPLACE ";" "," AUTHORS "${NUARG_AUTHORS}")
    string(APPEND NUSPEC_CONTENT "\n${NUSPEC_SUBELEMENT_INDENT_SIZE}<authors>${AUTHORS}</authors>")
    # Optional simple metadata subelements
    if(NOT "${NUARG_PROJECT_URL}" STREQUAL "")
        string(APPEND NUSPEC_CONTENT "\n${NUSPEC_SUBELEMENT_INDENT_SIZE}<projectUrl>${NUARG_PROJECT_URL}</projectUrl>")
    endif()
    if(NOT "${NUARG_ICON}" STREQUAL "")
        string(APPEND NUSPEC_CONTENT "\n${NUSPEC_SUBELEMENT_INDENT_SIZE}<icon>${NUARG_ICON}</icon>")
    endif()
    if(NOT "${NUARG_COPYRIGHT}" STREQUAL "")
        string(APPEND NUSPEC_CONTENT "\n${NUSPEC_SUBELEMENT_INDENT_SIZE}<copyright>${NUARG_COPYRIGHT}</copyright>")
    endif()
    # Optional complex metadata subelements
    # Begin /package/metadata/repository
    set(NUSPEC_REPOSITORY_CONTENT_BEGIN "\n${NUSPEC_SUBELEMENT_INDENT_SIZE}<repository")
    set(NUSPEC_REPOSITORY_CONTENT_END " />")
    set(NUSPEC_REPOSITORY_CONTENT "${NUSPEC_REPOSITORY_CONTENT_BEGIN}")
    # Attributes of /package/metadata/repository
    if(NOT "${NUARG_REPOSITORY_TYPE}" STREQUAL "")
        string(APPEND NUSPEC_REPOSITORY_CONTENT " type=\"${NUARG_REPOSITORY_TYPE}\"")
    endif()
    if(NOT "${NUARG_REPOSITORY_URL}" STREQUAL "")
        string(APPEND NUSPEC_REPOSITORY_CONTENT " url=\"${NUARG_REPOSITORY_URL}\"")
    endif()
    if(NOT "${NUARG_REPOSITORY_BRANCH}" STREQUAL "")
        string(APPEND NUSPEC_REPOSITORY_CONTENT " branch=\"${NUARG_REPOSITORY_BRANCH}\"")
    endif()
    if(NOT "${NUARG_REPOSITORY_COMMIT}" STREQUAL "")
        string(APPEND NUSPEC_REPOSITORY_CONTENT " commit=\"${NUARG_REPOSITORY_COMMIT}\"")
    endif()
    # End /package/metadata/repository
    string(APPEND NUSPEC_REPOSITORY_CONTENT "${NUSPEC_REPOSITORY_CONTENT_END}")
    if(NOT "${NUSPEC_REPOSITORY_CONTENT}" STREQUAL "${NUSPEC_REPOSITORY_CONTENT_BEGIN}${NUSPEC_REPOSITORY_CONTENT_END}")
        string(APPEND NUSPEC_CONTENT "${NUSPEC_REPOSITORY_CONTENT}")
    endif()
    # Optional collection metadata subelements
    # Section: /package/metadata/dependencies -- add package dependencies that are marked as PUBLIC or INTERFACE
    # in previous nuget_add_dependencies() calls.
    nuget_internal_nuspec_add_dependencies("${NUSPEC_SUBELEMENT_INDENT_SIZE}" "${NUSPEC_CONTENT}" NUSPEC_CONTENT)
    # End /package/metadata
    string(APPEND NUSPEC_CONTENT "\n${NUSPEC_INDENT_SIZE}</metadata>")
    set(${OUT_NUSPEC_CONTENT} "${NUSPEC_CONTENT}" PARENT_SCOPE)
    set(${OUT_PACKAGE_ID} "${NUARG_PACKAGE}" PARENT_SCOPE)
endfunction()

## Internal. Section: /package/metadata/dependencies in .nuspec XML file.
## Automatically generated based on previous nuget_add_dependencies() calls.
## Only dependencies marked as PUBLIC or INTERFACE are added (ie. PRIVATE dependencies are omitted).
function(nuget_internal_nuspec_add_dependencies NUSPEC_INDENT_SIZE NUSPEC_CONTENT OUT_NUSPEC_CONTENT)
    # Begin /package/metadata/dependencies
    set(NUSPEC_DEPENDENCIES_CONTENT_BEGIN "\n${NUSPEC_INDENT_SIZE}<dependencies>")
    set(NUSPEC_DEPENDENCIES_CONTENT_END "\n${NUSPEC_INDENT_SIZE}</dependencies>")
    set(NUSPEC_DEPENDENCIES_CONTENT "${NUSPEC_DEPENDENCIES_CONTENT_BEGIN}")
    # For each dependency that should be in /package/metadata/dependencies
    nuget_get_dependencies(DEPENDENCIES)
    set(NUSPEC_SUBELEMENT_INDENT_SIZE "${NUSPEC_INDENT_SIZE}${NUGET_NUSPEC_INDENT_SIZE}")
    foreach(DEPENDENCY IN LISTS DEPENDENCIES)
        nuget_get_dependency_usage("${DEPENDENCY}" USAGE)
        if("${USAGE}" STREQUAL "PRIVATE")
            continue()
        endif()
        nuget_get_dependency_version("${DEPENDENCY}" VERSION)
        string(APPEND NUSPEC_DEPENDENCIES_CONTENT "\n${NUSPEC_SUBELEMENT_INDENT_SIZE}<dependency id=\"${DEPENDENCY}\" version=\"${VERSION}\" />")
    endforeach()
    # End /package/metadata/dependencies
    string(APPEND NUSPEC_DEPENDENCIES_CONTENT "${NUSPEC_DEPENDENCIES_CONTENT_END}")
    if(NOT "${NUSPEC_DEPENDENCIES_CONTENT}" STREQUAL "${NUSPEC_DEPENDENCIES_CONTENT_BEGIN}${NUSPEC_DEPENDENCIES_CONTENT_END}")
        string(APPEND NUSPEC_CONTENT "${NUSPEC_DEPENDENCIES_CONTENT}")
    endif()
    set(${OUT_NUSPEC_CONTENT} "${NUSPEC_CONTENT}" PARENT_SCOPE)
endfunction()

## Internal. Section: /package/files in .nuspec XML file (FILES as section identifier CMake argument).
function(nuget_internal_nuspec_process_files_args NUSPEC_INDENT_SIZE NUSPEC_CONTENT OUT_NUSPEC_CONTENT)
    # Begin /package/files
    set(NUSPEC_FILES_CONTENT_BEGIN "\n${NUSPEC_INDENT_SIZE}<files>")
    set(NUSPEC_FILES_CONTENT_END "\n${NUSPEC_INDENT_SIZE}</files>")
    set(NUSPEC_FILES_CONTENT "${NUSPEC_FILES_CONTENT_BEGIN}")
    set(ARGS_HEAD "")
    set(ARGS_TAIL ${ARGN})
    set(NUSPEC_SUBELEMENT_INDENT_SIZE "${NUSPEC_INDENT_SIZE}${NUGET_NUSPEC_INDENT_SIZE}")
    while(NOT "${ARGS_TAIL}" STREQUAL "")
        nuget_internal_helper_cut_arg_list(CMAKE_CONDITIONAL_SECTION "${ARGS_TAIL}" ARGS_HEAD ARGS_TAIL)
        list(LENGTH ARGS_HEAD ARGS_HEAD_LENGTH)
        if(ARGS_HEAD_LENGTH GREATER_EQUAL 2)
            list(GET ARGS_HEAD 0 MAYBE_CMAKE_INCLUDE_CONDITION_IDENTIFIER)
            if("${MAYBE_CMAKE_INCLUDE_CONDITION_IDENTIFIER}" STREQUAL "CMAKE_CONDITIONAL_SECTION")
                list(GET ARGS_HEAD 1 CMAKE_CONDITIONAL_SECTION)
                nuget_internal_helper_list_sublist("${ARGS_HEAD}" 2 -1 ARGS_HEAD)
            endif()
        endif()
        nuget_internal_nuspec_add_files_conditionally("${NUSPEC_SUBELEMENT_INDENT_SIZE}" "${NUSPEC_FILES_CONTENT}" NUSPEC_FILES_CONTENT
            "${CMAKE_CONDITIONAL_SECTION}" ${ARGS_HEAD}
        )
    endwhile()
    # End /package/files
    string(APPEND NUSPEC_FILES_CONTENT "${NUSPEC_FILES_CONTENT_END}")
    if("${NUSPEC_FILES_CONTENT}" STREQUAL "${NUSPEC_FILES_CONTENT_BEGIN}${NUSPEC_FILES_CONTENT_END}")
        message(FATAL_ERROR "Assembled expression for generating /package/files node(s) for .nuspec XML file(s) is empty.")
    endif()
    string(APPEND NUSPEC_CONTENT "${NUSPEC_FILES_CONTENT}")
    set(${OUT_NUSPEC_CONTENT} "${NUSPEC_CONTENT}" PARENT_SCOPE)
endfunction()

## Internal.
function(nuget_internal_nuspec_add_files_conditionally NUSPEC_INDENT_SIZE NUSPEC_CONTENT OUT_NUSPEC_CONTENT CMAKE_CONDITIONAL_SECTION)
    # Input: check for a CMAKE_CONDITIONAL_SECTION parameter pack
    if(NOT "${CMAKE_CONDITIONAL_SECTION}" STREQUAL "")
        string(APPEND NUSPEC_CONTENT "$<${CMAKE_CONDITIONAL_SECTION}:")
    endif()
    nuget_internal_helper_error_if_empty("${ARGN}" "Input expression for generating elements under /package/files node(s) for .nuspec XML file(s) is empty: no FILE_SRC and FILE_TARGET arguments were provided.")
    # Loop over parameter pack
    set(ARGS_HEAD "")
    set(ARGS_TAIL ${ARGN})
    while(NOT "${ARGS_TAIL}" STREQUAL "")
        nuget_internal_helper_cut_arg_list(FILE_SRC "${ARGS_TAIL}" ARGS_HEAD ARGS_TAIL)
        nuget_internal_nuspec_add_file_conditionally("${NUSPEC_INDENT_SIZE}" "${NUSPEC_CONTENT}" NUSPEC_CONTENT
            "${CMAKE_CONDITIONAL_SECTION}" ${ARGS_HEAD}
        )
    endwhile()
    # Close generator expression if this was a CMAKE_CONDITIONAL_SECTION parameter pack
    if(NOT "${CMAKE_CONDITIONAL_SECTION}" STREQUAL "")
        string(APPEND NUSPEC_CONTENT ">")
    endif()
    set(${OUT_NUSPEC_CONTENT} "${NUSPEC_CONTENT}" PARENT_SCOPE)
endfunction()

## Internal.
function(nuget_internal_nuspec_add_file_conditionally NUSPEC_INDENT_SIZE NUSPEC_CONTENT OUT_NUSPEC_CONTENT CMAKE_CONDITIONAL_SECTION)
    # Inputs
    # See https://docs.microsoft.com/en-us/nuget/reference/nuspec#file-element-attributes
    set(options "")
    set(oneValueArgs FILE_SRC FILE_TARGET)
    set(multiValueArgs FILE_EXCLUDE)
    cmake_parse_arguments(NUARG
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN}
    )
    nuget_internal_helper_error_if_unparsed_args(
        "${NUARG_UNPARSED_ARGUMENTS}"
        "${NUARG_KEYWORDS_MISSING_VALUES}"
    )
    nuget_internal_helper_error_if_empty("${NUARG_FILE_SRC}"
        "FILE_SRC must not be empty: it is a required attribute (src) of "
        "a .nuspec XML file's /package/files/file element."
    )
    # Actual functionality
    string(APPEND NUSPEC_CONTENT "\n${NUSPEC_INDENT_SIZE}<file")
    string(APPEND NUSPEC_CONTENT " src=\"${NUARG_FILE_SRC}\"")
    if(NOT "${NUARG_FILE_TARGET}" STREQUAL "")
        string(APPEND NUSPEC_CONTENT " target=\"${NUARG_FILE_TARGET}\"")
    endif()
    if(NOT "${NUARG_FILE_EXCLUDE}" STREQUAL "")
        string(APPEND NUSPEC_CONTENT " exclude=\"${NUARG_FILE_EXCLUDE}\"")
    endif()
    string(APPEND NUSPEC_CONTENT " /$<ANGLE-R>")
    set(${OUT_NUSPEC_CONTENT} "${NUSPEC_CONTENT}" PARENT_SCOPE)
endfunction()

## Internal. Section: CMake-specific (without special section identifier CMake argument).
## Write output .nuspec XML file(s) conditionally for provided configurations in CMAKE_CONFIGURATIONS intersected with
## the available configurations in the current build system this function is actually called from. No error is raised if
## a given configuration is not available -- the output file is simply not generated for that in the current build system.
## Not raising an error if a given configuration is unavailable makes it possible to reuse the same nuget_generate_nuspec_files()
## calls across different build systems without adjustments or writing additional code for generating the values of the
## CMAKE_CONFIGURATIONS argument.
function(nuget_internal_nuspec_generate_output NUSPEC_CONTENT PACKAGE_ID CMAKE_OUTPUT_DIR CMAKE_CONFIGURATIONS)
    # Inputs
    nuget_internal_helper_error_if_empty("${NUSPEC_CONTENT}" "NUSPEC_CONTENT to be written is empty: cannot generate .nuspec file's content.")
    nuget_internal_helper_error_if_empty("${PACKAGE_ID}" "PACKAGE_ID to be written is empty: cannot generate .nuspec filename.")
    if("${CMAKE_OUTPUT_DIR}" STREQUAL "")
        set(OUTPUT_FILE "${CMAKE_BINARY_DIR}/CMakeNuGetTools/nuspec/${PACKAGE_ID}.$<CONFIG>.nuspec")
    else()
        set(OUTPUT_FILE "${CMAKE_OUTPUT_DIR}/${PACKAGE_ID}.$<CONFIG>.nuspec")
    endif()
    # Actual functionality
    if("${CMAKE_CONFIGURATIONS}" STREQUAL "")
        file(GENERATE OUTPUT "${OUTPUT_FILE}" CONTENT "${NUSPEC_CONTENT}")
        message(STATUS "Written \"${OUTPUT_FILE}\" file(s).")
    else()
        set(CONDITIONS "$<OR:")
        foreach(CONFIGURATION IN LISTS CMAKE_CONFIGURATIONS)
            string(APPEND CONDITIONS "${CONDITIONS_SEPARATOR}$<CONFIG:${CONFIGURATION}>")
            set(CONDITIONS_SEPARATOR ",")
        endforeach()
        string(APPEND CONDITIONS ">")
        file(GENERATE OUTPUT "${OUTPUT_FILE}" CONTENT "${NUSPEC_CONTENT}" CONDITION "${CONDITIONS}")
        message(STATUS "Written \"${OUTPUT_FILE}\" file(s) for \"${CMAKE_CONFIGURATIONS}\" configuration(s).")
    endif()
endfunction()

## Internal.
function(nuget_internal_merge_second_nuspec_file_into_first FILEPATH_ACC FILEPATH_IN)
    set(FILES_NODE_BEGIN_STR "<files>")
    set(FILES_NODE_END_STR "</files>")
    string(LENGTH "${FILES_NODE_BEGIN_STR}" FILES_NODE_BEGIN_LEN)
    string(LENGTH "${FILES_NODE_END_STR}" FILES_NODE_END_LEN)
    # Inputs: FILEPATH_IN
    file(STRINGS "${FILEPATH_IN}" LINES_IN NEWLINE_CONSUME ENCODING UTF-8)
    string(FIND "${LINES_IN}" "${FILES_NODE_BEGIN_STR}" LINES_IN_FILES_NODE_BEGIN_POS)
    if(${LINES_IN_FILES_NODE_BEGIN_POS} EQUAL -1)
        message(FATAL_ERROR "Cannot merge: did not find the \"${FILES_NODE_BEGIN_STR}\" part of the .nuspec files node in \"${FILEPATH_IN}\".")
    endif()
    string(FIND "${LINES_IN}" "${FILES_NODE_END_STR}" LINES_IN_FILES_NODE_END_POS REVERSE)
    if(${LINES_IN_FILES_NODE_END_POS} EQUAL -1)
        message(FATAL_ERROR "Cannot merge: did not find the \"${FILES_NODE_END_STR}\" part of the .nuspec files node in \"${FILEPATH_IN}\".")
    endif()
    # Inputs: FILEPATH_ACC
    file(STRINGS "${FILEPATH_ACC}" LINES_ACC NEWLINE_CONSUME ENCODING UTF-8)
    string(FIND "${LINES_ACC}" "${FILES_NODE_BEGIN_STR}" LINES_ACC_FILES_NODE_BEGIN_POS)
    if(${LINES_ACC_FILES_NODE_BEGIN_POS} EQUAL -1)
        message(FATAL_ERROR "Cannot merge: did not find the \"${FILES_NODE_BEGIN_STR}\" part of the .nuspec files node in \"${FILEPATH_ACC}\".")
    endif()
    string(FIND "${LINES_ACC}" "${FILES_NODE_END_STR}" LINES_ACC_FILES_NODE_END_POS REVERSE)
    if(${LINES_ACC_FILES_NODE_END_POS} EQUAL -1)
        message(FATAL_ERROR "Cannot merge: did not find the \"${FILES_NODE_END_STR}\" part of the .nuspec files node in \"${FILEPATH_ACC}\".")
    endif()
    # Check: FILEPATH_ACC and FILEPATH_IN contents outside the files node should not differ
    string(SUBSTRING "${LINES_IN}" 0 ${LINES_IN_FILES_NODE_BEGIN_POS} LINES_IN_UNTIL_FILES_NODE_BEGIN)
    string(SUBSTRING "${LINES_ACC}" 0 ${LINES_ACC_FILES_NODE_BEGIN_POS} LINES_ACC_UNTIL_FILES_NODE_BEGIN)
    if(NOT "${LINES_IN_UNTIL_FILES_NODE_BEGIN}" STREQUAL "${LINES_ACC_UNTIL_FILES_NODE_BEGIN}")
        message(FATAL_ERROR "Cannot merge: file content before the .nuspec files node of \"${FILEPATH_ACC}\" and \"${FILEPATH_IN}\" differs.")
    endif()
    string(SUBSTRING "${LINES_IN}" ${LINES_IN_FILES_NODE_END_POS} -1 LINES_IN_AFTER_FILES_NODE_END)
    string(SUBSTRING "${LINES_ACC}" ${LINES_ACC_FILES_NODE_END_POS} -1 LINES_ACC_AFTER_FILES_NODE_END)
    if(NOT "${LINES_IN_AFTER_FILES_NODE_END}" STREQUAL "${LINES_ACC_AFTER_FILES_NODE_END}")
        message(FATAL_ERROR "Cannot merge: file content after the .nuspec files node of \"${FILEPATH_ACC}\" and \"${FILEPATH_IN}\" differs.")
    endif()
    # Create merged content
    # NOTE: no need to check for duplicate <file> element entries when merging FILEPATH_ACC and FILEPATH_IN: nuget pack does not
    # seem to complain when file elements with the same src and target attributes are present in the files node of a .nuspec file.
    string(SUBSTRING "${LINES_ACC}" 0 ${LINES_ACC_FILES_NODE_END_POS} NEW_LINES_ACC)
    string(APPEND NEW_LINES_ACC "${NUGET_NUSPEC_INDENT_SIZE}<!-- Below merged from: \"${FILEPATH_IN}\" -->")
    math(EXPR LINES_IN_AFTER_FILES_NODE_BEGIN_POS "${LINES_IN_FILES_NODE_BEGIN_POS} + ${FILES_NODE_BEGIN_LEN}")
    string(SUBSTRING "${LINES_IN}" ${LINES_IN_AFTER_FILES_NODE_BEGIN_POS} -1 LINES_IN_AFTER_FILES_NODE_BEGIN)
    string(APPEND NEW_LINES_ACC "${LINES_IN_AFTER_FILES_NODE_BEGIN}")
    # Output: overwrite FILEPATH_ACC with merged content
    file(WRITE "${FILEPATH_ACC}" "${NEW_LINES_ACC}")
endfunction()

## Internal.
function(nuget_internal_merge_n_nuspec_files FILEPATH_ACC)
    nuget_internal_helper_error_if_empty("${FILEPATH_ACC}" "No filepath to a .nuspec file provided as a basis for merge operation.")
    nuget_internal_helper_error_if_empty("${ARGN}" "No .nuspec filepaths provided for merge operation.")
    # Read first input file (FILEPATH_BASE)
    list(GET ARGN 0 FILEPATH_BASE)
    file(STRINGS "${FILEPATH_BASE}" LINES_BASE NEWLINE_CONSUME ENCODING UTF-8)
    # Inputs: FILEPATH_BASE
    set(FILES_NODE_BEGIN_STR "<files>")
    string(LENGTH "${FILES_NODE_BEGIN_STR}" FILES_NODE_BEGIN_LEN)
    string(FIND "${LINES_BASE}" "${FILES_NODE_BEGIN_STR}" LINES_BASE_FILES_NODE_BEGIN_POS)
    if(${LINES_BASE_FILES_NODE_BEGIN_POS} EQUAL -1)
        message(FATAL_ERROR "Cannot merge: did not find the \"${FILES_NODE_BEGIN_STR}\" part of the .nuspec files node in \"${FILEPATH_BASE}\".")
    endif()
    # Initialize content of FILEPATH_ACC with FILEPATH_BASE adding comment referring to original source file
    math(EXPR LINES_BASE_AFTER_FILES_NODE_BEGIN_POS "${LINES_BASE_FILES_NODE_BEGIN_POS} + ${FILES_NODE_BEGIN_LEN}")
    string(SUBSTRING "${LINES_BASE}" 0 ${LINES_BASE_AFTER_FILES_NODE_BEGIN_POS} NEW_LINES_BASE)
    string(APPEND NEW_LINES_BASE "\n${NUGET_NUSPEC_INDENT_SIZE}${NUGET_NUSPEC_INDENT_SIZE}<!-- Below loaded from: \"${FILEPATH_BASE}\" -->")
    string(SUBSTRING "${LINES_BASE}" ${LINES_BASE_AFTER_FILES_NODE_BEGIN_POS} -1 LINES_BASE_AFTER_FILES_NODE_BEGIN)
    string(APPEND NEW_LINES_BASE "${LINES_BASE_AFTER_FILES_NODE_BEGIN}")
    file(WRITE "${FILEPATH_ACC}" "${NEW_LINES_BASE}")
    # Merge rest of the input files into FILEPATH_ACC
    nuget_internal_helper_list_sublist("${ARGN}" 1 -1 FILEPATHS_IN_TAIL)
    foreach(FILEPATH_IN IN LISTS FILEPATHS_IN_TAIL)
        nuget_internal_merge_second_nuspec_file_into_first("${FILEPATH_ACC}" "${FILEPATH_IN}")
    endforeach()
endfunction()

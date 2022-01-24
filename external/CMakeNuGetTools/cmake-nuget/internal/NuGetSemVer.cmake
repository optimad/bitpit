## Include implementation
include("${CMAKE_CURRENT_LIST_DIR}/NuGetSemVer.git.cmake")

## Public interface.
function(nuget_git_get_semantic_version)
    set(options "")
    set(oneValueArgs TAG_PREFIX PRERELEASE_LABEL FULL CORE MAJOR MINOR PATCH PRERELEASE)
    set(multiValueArgs "")
    cmake_parse_arguments(NUARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGV})
    nuget_internal_helper_error_if_unparsed_args("${NUARG_UNPARSED_ARGUMENTS}" "${NUARG_KEYWORDS_MISSING_VALUES}")
    nuget_internal_git_get_semantic_version_with_prerelease_override(
        "${NUARG_TAG_PREFIX}"
        "${NUARG_PRERELEASE_LABEL}"
        MAJOR
        MINOR
        PATCH
        PRERELEASE
    )
    if(NOT "${NUARG_FULL}" STREQUAL "")
        if(NOT "${PRERELEASE}" STREQUAL "")
            set(${NUARG_FULL} "${MAJOR}.${MINOR}.${PATCH}-${PRERELEASE}" PARENT_SCOPE)
        else()
            set(${NUARG_FULL} "${MAJOR}.${MINOR}.${PATCH}" PARENT_SCOPE)
        endif()
    endif()
    if(NOT "${NUARG_CORE}" STREQUAL "")
        set(${NUARG_CORE} "${MAJOR}.${MINOR}.${PATCH}" PARENT_SCOPE)
    endif()
    if(NOT "${NUARG_MAJOR}" STREQUAL "")
        set(${NUARG_MAJOR} "${MAJOR}" PARENT_SCOPE)
    endif()
    if(NOT "${NUARG_MINOR}" STREQUAL "")
        set(${NUARG_MINOR} "${MINOR}" PARENT_SCOPE)
    endif()
    if(NOT "${NUARG_PATCH}" STREQUAL "")
        set(${NUARG_PATCH} "${PATCH}" PARENT_SCOPE)
    endif()
    if(NOT "${NUARG_PRERELEASE}" STREQUAL "")
        set(${NUARG_PRERELEASE} "${PRERELEASE}" PARENT_SCOPE)
    endif()
endfunction()

## Public interface.
function(nuget_git_get_mapped_semantic_version)
    set(options "")
    set(oneValueArgs TAG_PREFIX
        BRANCH FULL CORE MAJOR MINOR PATCH PRERELEASE
    )
    set(multiValueArgs BRANCH_NAME_REGEXES PRERELEASE_PREFIX_LABELS PRERELEASE_POSTFIX_FLAGS
        NO_PRELEASE_WHEN_ON_TAG_FLAGS COMMIT_COUNT_IN_PRELEASE_FLAGS COMMIT_COUNT_ONLY_ADDS_ONE_FLAGS)
    cmake_parse_arguments(NUARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGV})
    nuget_internal_helper_error_if_unparsed_args("${NUARG_UNPARSED_ARGUMENTS}" "${NUARG_KEYWORDS_MISSING_VALUES}")
    nuget_internal_git_get_semantic_version_applying_rules(
        "${NUARG_TAG_PREFIX}"
        "${NUARG_BRANCH_NAME_REGEXES}"
        "${NUARG_PRERELEASE_PREFIX_LABELS}"
        "${NUARG_PRERELEASE_POSTFIX_FLAGS}"
        "${NUARG_NO_PRELEASE_WHEN_ON_TAG_FLAGS}"
        "${NUARG_COMMIT_COUNT_IN_PRELEASE_FLAGS}"
        "${NUARG_COMMIT_COUNT_ONLY_ADDS_ONE_FLAGS}"
        BRANCH
        MAJOR
        MINOR
        PATCH
        PRERELEASE
    )
    if(NOT "${NUARG_BRANCH}" STREQUAL "")
        set(${NUARG_BRANCH} "${BRANCH}" PARENT_SCOPE)
    endif()
    if(NOT "${NUARG_FULL}" STREQUAL "")
        if(NOT "${PRERELEASE}" STREQUAL "")
            set(${NUARG_FULL} "${MAJOR}.${MINOR}.${PATCH}-${PRERELEASE}" PARENT_SCOPE)
        else()
            set(${NUARG_FULL} "${MAJOR}.${MINOR}.${PATCH}" PARENT_SCOPE)
        endif()
    endif()
    if(NOT "${NUARG_CORE}" STREQUAL "")
        set(${NUARG_CORE} "${MAJOR}.${MINOR}.${PATCH}" PARENT_SCOPE)
    endif()
    if(NOT "${NUARG_MAJOR}" STREQUAL "")
        set(${NUARG_MAJOR} "${MAJOR}" PARENT_SCOPE)
    endif()
    if(NOT "${NUARG_MINOR}" STREQUAL "")
        set(${NUARG_MINOR} "${MINOR}" PARENT_SCOPE)
    endif()
    if(NOT "${NUARG_PATCH}" STREQUAL "")
        set(${NUARG_PATCH} "${PATCH}" PARENT_SCOPE)
    endif()
    if(NOT "${NUARG_PRERELEASE}" STREQUAL "")
        set(${NUARG_PRERELEASE} "${PRERELEASE}" PARENT_SCOPE)
    endif()
endfunction()

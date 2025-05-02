
# Configure self containment for the specified target
#
# \param TARGET_NAME is the name of the target for which the self-containment will be configured
# \param TARGET_INSTALL_DIR is the directory where the target will be installed, the directory
# can be an absolute path or a path relative to the CMake install prefix
# \param THIRD_PARTY_PATTERNS is a semi-color separated list of patterns that will be used to
# identify the third party libraries
# \param THIRD_PARTY_INSTALL_DIR is the directory where the third-party libraries will be
# installed, the directory can be an absolute path or a path relative to the CMake install
# prefix
function(configure_self_containment TARGET_NAME TARGET_INSTALL_DIR THIRD_PARTY_PATTERNS THIRD_PARTY_INSTALL_DIR)
    get_self_containment_additional_search_dirs(ADDITIONAL_SEARCH_DIRS)

    set(COLLECT_THIRD_PARTY_LIBRARIES_TARGET ${TARGET_NAME}_self_containment_third_party_collection)

    string(REPLACE ";" "$<SEMICOLON>" THIRD_PARTY_PATTERNS "${THIRD_PARTY_PATTERNS}")

    if(NOT IS_ABSOLUTE "${THIRD_PARTY_INSTALL_DIR}")
        get_filename_component(ABSOLUTE_THIRD_PARTY_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${THIRD_PARTY_INSTALL_DIR}" ABSOLUTE)
    else()
        set(ABSOLUTE_THIRD_PARTY_INSTALL_DIR "${THIRD_PARTY_INSTALL_DIR}")
    endif()

    get_target_property(TARGET_TYPE ${TARGET_NAME} TYPE)
    if(TARGET_TYPE STREQUAL "EXECUTABLE")
        set(LIBRARY_FILES "")
        set(EXECUTABLE_FILES "$<TARGET_FILE:${TARGET_NAME}>")
    else()
        set(LIBRARY_FILES "$<TARGET_FILE:${TARGET_NAME}>")
        set(EXECUTABLE_FILES "")
    endif()

    if (${WIN32})
        set(EXCLUDED_DEPENDENCIES "$<JOIN:$<TARGET_RUNTIME_DLLS:${TARGET_NAME}>,$<SEMICOLON>>")
    else()
        set(EXCLUDED_DEPENDENCIES "")
    endif()

    add_custom_target(${COLLECT_THIRD_PARTY_LIBRARIES_TARGET} COMMAND
        ${CMAKE_COMMAND} -E echo "Collecting third-party libraries for target ${TARGET_NAME}" &&
        ${CMAKE_COMMAND} -E echo "Third-party libraries will be copied to ${ABSOLUTE_THIRD_PARTY_INSTALL_DIR}" &&
        ${CMAKE_COMMAND}
        -D executable_list="${EXECUTABLE_FILES}"
        -D library_list="${LIBRARY_FILES}"
        -D include_filename_patterns="${THIRD_PARTY_PATTERNS}"
        -D exclude_filepath_list="${EXCLUDED_DEPENDENCIES}"
        -D copy_to_directory="${ABSOLUTE_THIRD_PARTY_INSTALL_DIR}"
        -D additional_search_directories="${ADDITIONAL_SEARCH_DIRS}"
        -P ${PROJECT_SOURCE_DIR}/cmake/CopyDependencies.cmake
    )

    install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} --build . --config $<CONFIG> --target ${COLLECT_THIRD_PARTY_LIBRARIES_TARGET})")

    set_run_path(${TARGET_NAME} "${TARGET_INSTALL_DIR}" "${CMAKE_INSTALL_LIBDIR};${THIRD_PARTY_INSTALL_DIR}")
endfunction()

# Get additional search directories in which to search for third-party libraries.
#
# \param ADDITIONAL_SEARCH_DIRS on output will contain the search directories
function(get_self_containment_additional_search_dirs ADDITIONAL_SEARCH_DIRS)
    if(WIN32 AND MKL_FOUND)
        set(${ADDITIONAL_SEARCH_DIRS} "${OMP_DLL_DIR}" PARENT_SCOPE)
    elseif(NOT "$ENV{LD_LIBRARY_PATH}" STREQUAL "")
        string(REPLACE ":" "$<SEMICOLON>" ENV_SEARCH_DIRS $ENV{LD_LIBRARY_PATH})
        set(${ADDITIONAL_SEARCH_DIRS} ${ENV_SEARCH_DIRS} PARENT_SCOPE)
    else()
        set(${ADDITIONAL_SEARCH_DIRS} "." PARENT_SCOPE)
    endif()
endfunction()

# Set the run path (rpath) option for the specified target
#
# \param TARGET_NAME is the name of the target for which the run path will be set
# \param ORIGIN_DIR is the directory where the target will be installed, the directory can be
# an absolute path or a path relative to the CMake install prefix
# \param LIBRARY_DIRS are the directories that will be added to the run path, directories can
# be an absolute paths or a paths relative to the CMake install prefix
function(set_run_path TARGET_NAME ORIGIN_DIR LIBRARY_DIRS)
    if(NOT IS_ABSOLUTE "${ORIGIN_DIR}")
        get_filename_component(ABSOLUTE_ORIGIN_DIR "${CMAKE_INSTALL_PREFIX}/${ORIGIN_DIR}" ABSOLUTE)
    else()
        set(ABSOLUTE_ORIGIN_DIR "${ORIGIN_DIR}")
    endif()

    set(INSTALL_RPATH "$ORIGIN")
    foreach (LIBRARY_DIR IN LISTS LIBRARY_DIRS)
        if(NOT IS_ABSOLUTE "${LIBRARY_DIR}")
            get_filename_component(ABSOLUTE_LIBRARY_DIR "${CMAKE_INSTALL_PREFIX}/${LIBRARY_DIR}" ABSOLUTE)
        else()
            set(ABSOLUTE_LIBRARY_DIR "${LIBRARY_DIR}")
        endif()

        file(RELATIVE_PATH RELATIVE_LIBRARY_DIR ${ABSOLUTE_ORIGIN_DIR} ${ABSOLUTE_LIBRARY_DIR})
        if(NOT "${RELATIVE_LIBRARY_DIR}" STREQUAL "")
            set(INSTALL_RPATH "${INSTALL_RPATH}:$ORIGIN/${RELATIVE_LIBRARY_DIR}")
        endif()
    endforeach()

    set_target_properties(${TARGET_NAME} PROPERTIES SKIP_BUILD_RPATH FALSE)
    set_target_properties(${TARGET_NAME} PROPERTIES BUILD_RPATH_USE_ORIGIN TRUE)
    set_target_properties(${TARGET_NAME} PROPERTIES INSTALL_RPATH "${INSTALL_RPATH}")
    set_target_properties(${TARGET_NAME} PROPERTIES INSTALL_RPATH_USE_LINK_PATH FALSE)
endfunction()

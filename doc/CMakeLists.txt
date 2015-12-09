# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Add a target to generate API documentation with Doxygen
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" OFF)

IF(BUILD_DOCUMENTATION)
  FIND_PACKAGE(Doxygen)
  IF(NOT DOXYGEN_FOUND)
    MESSAGE(FATAL_ERROR "Doxygen is needed to build the documentation.")
  ENDIF()

  option(DOC_EXTRACT_PRIVATE "If set, all private members are included in the documentation" OFF)
  IF (DOC_EXTRACT_PRIVATE)
    SET(DOXY_EXTRACT_PRIVATE "YES")
  ELSE ()
    SET(DOXY_EXTRACT_PRIVATE "NO")
  ENDIF ()

  SET( doxyfile_in          ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in     )
  SET( doxyfile             ${PROJECT_BINARY_DIR}/Doxyfile              )
  SET( doxy_html_index_file ${CMAKE_CURRENT_BINARY_DIR}/html/index.html )
  SET( doxy_output_root     ${CMAKE_CURRENT_BINARY_DIR}                 ) # Pasted into Doxyfile.in
  SET( doxy_input           ${PROJECT_SOURCE_DIR}/src                   ) # Pasted into Doxyfile.in

  if (NOT "${DOC_EXTRACT_PRIVATE_INTERNAL}" STREQUAL "${DOC_EXTRACT_PRIVATE}")
    if (EXISTS ${doxy_html_index_file})
      FILE (REMOVE ${doxy_html_index_file})
    endif ()
    set(DOC_EXTRACT_PRIVATE_INTERNAL ${DOC_EXTRACT_PRIVATE} CACHE INTERNAL "This is the value of the last time DOC_EXTRACT_PRIVATE was set successfully." FORCE)
  endif ()

  CONFIGURE_FILE( ${doxyfile_in} ${doxyfile} @ONLY )

  ADD_CUSTOM_COMMAND( OUTPUT ${doxy_html_index_file}
                      COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
                      # The following should be ${doxyfile} only but it
                      # will break the dependency.
                      # The optimal solution would be creating a
                      # custom_command for ${doxyfile} generation
                      # but I still have to figure out how...
                      MAIN_DEPENDENCY ${doxyfile} ${doxyfile_in}
                      DEPENDS ${doxy_extra_files} ${doxyfile_in}
                              ${COMMON_HEADERS}
                              ${BASE_SOURCES} ${BASE_HEADERS}
                              ${IO_SOURCES} ${IO_HEADERS}
                              ${SA_SOURCES} ${SA_HEADERS}
                              ${EXAMPLES_SOURCES}
                      COMMENT "Generating HTML documentation")

  ADD_CUSTOM_TARGET( doc ALL DEPENDS ${doxy_html_index_file} )
  add_custom_target(clean-doc COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  INSTALL( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION doc/${CMAKE_PROJECT_NAME}-${BITP_MESH_VERSION} )
ENDIF()

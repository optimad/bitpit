# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Add a target to generate API documentation with Doxygen
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" OFF)

IF(BUILD_DOCUMENTATION)
  FIND_PACKAGE(Doxygen)
  IF(NOT DOXYGEN_FOUND)
    MESSAGE(FATAL_ERROR "Doxygen is needed to build the documentation.")
  ENDIF()

  SET( doxyfile_in          ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in     )
  SET( doxyfile             ${PROJECT_BINARY_DIR}/Doxyfile              )
  SET( doxy_html_index_file ${CMAKE_CURRENT_BINARY_DIR}/html/index.html )
  SET( doxy_output_root     ${CMAKE_CURRENT_BINARY_DIR}                 ) # Pasted into Doxyfile.in
  SET( doxy_input           ${PROJECT_SOURCE_DIR}/src                   ) # Pasted into Doxyfile.in

  CONFIGURE_FILE( ${doxyfile_in} ${doxyfile} @ONLY )

  ADD_CUSTOM_COMMAND( OUTPUT ${doxy_html_index_file}
                      COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
                      # The following should be ${doxyfile} only but it
                      # will break the dependency.
                      # The optimal solution would be creating a
                      # custom_command for ${doxyfile} generation
                      # but I still have to figure out how...
                      MAIN_DEPENDENCY ${doxyfile} ${doxyfile_in}
                      DEPENDS ${PATCHMAN_SOURCES} ${PATCHMAN_HEADERS}
                      COMMENT "Generating HTML documentation")

  ADD_CUSTOM_TARGET( doc ALL DEPENDS ${doxy_html_index_file} )

  INSTALL( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION doc/${CMAKE_PROJECT_NAME} )
ENDIF()
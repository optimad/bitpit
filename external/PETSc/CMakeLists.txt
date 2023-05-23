
#------------------------------------------------------------------------------------#
# CMake targets
#------------------------------------------------------------------------------------#
add_custom_target(clean-external-petsc COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(FILES "cmake/FindPETSc.cmake" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")
install(FILES "cmake/ResolveCompilerPaths.cmake" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")
install(FILES "cmake/FindPackageMultipass.cmake" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")
install(FILES "cmake/CorrectWindowsPaths.cmake" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")

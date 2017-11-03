# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------------#
# CMake targets
#------------------------------------------------------------------------------------#
add_custom_target(clean-external-petsc COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(FILES "cmake/FindPETSc.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${BITPIT_INSTALL_CMAKE_DIR}")
install(FILES "cmake/ResolveCompilerPaths.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${BITPIT_INSTALL_CMAKE_DIR}")
install(FILES "cmake/FindPackageMultipass.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${BITPIT_INSTALL_CMAKE_DIR}")
install(FILES "cmake/CorrectWindowsPaths.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${BITPIT_INSTALL_CMAKE_DIR}")

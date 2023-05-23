
#------------------------------------------------------------------------------------#
# CMake targets
#------------------------------------------------------------------------------------#
add_custom_target(clean-external-lapacke COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(FILES "cmake/FindLAPACKE.cmake" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")


#------------------------------------------------------------------------------------#
# CMake targets
#------------------------------------------------------------------------------------#
add_custom_target(clean-external-metis COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(FILES "cmake/FindMETIS.cmake" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")

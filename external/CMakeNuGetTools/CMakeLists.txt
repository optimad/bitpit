# Specify the version being used as well as the language
cmake_minimum_required(VERSION 3.13)

#------------------------------------------------------------------------------------#
# CMake targets
#------------------------------------------------------------------------------------#
add_custom_target(clean-external-cmakenugettools COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(DIRECTORY "cmake-nuget" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")
install(FILES "bitpitNuGet.cmake" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")
# Install script for directory: /home/haysam/Dev/BitP_Mesh/src/patchman

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/haysam/bitpit")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/BitP_Mesh" TYPE FILE FILES
    "/home/haysam/Dev/BitP_Mesh/src/patchman/patchman_version.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/patchman/cell.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/patchman/reference.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/patchman/BitP_Mesh_PATCHMAN.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/patchman/patch_cartesian.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/patchman/element.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/patchman/patchman.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/patchman/output_manager.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/patchman/node.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/patchman/interface.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/patchman/collapsedArray2D.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/patchman/patch_octree.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/patchman/piercedVector.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/patchman/patch.hpp"
    )
endif()


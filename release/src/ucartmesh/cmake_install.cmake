# Install script for directory: /home/haysam/Dev/BitP_Mesh/src/ucartmesh

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
    "/home/haysam/Dev/BitP_Mesh/src/ucartmesh/BitP_Mesh_UCARTMESH.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/ucartmesh/Class_UCartMesh.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/ucartmesh/UCartMesh.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/ucartmesh/Class_UCartMesh.tpp"
    "/home/haysam/Dev/BitP_Mesh/src/ucartmesh/UCartMesh.tpp"
    )
endif()


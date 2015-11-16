# Install script for directory: /home/haysam/Dev/BitP_Mesh/src/pablo

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
    "/home/haysam/Dev/BitP_Mesh/src/pablo/mpi_datatype_conversion.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Data_Comm_Interface.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Intersection.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Map.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Octant.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/ioFunct.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Data_LB_Interface.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Array.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/PABLO_version.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Local_Tree.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Comm_Buffer.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/BitP_Mesh_PABLO.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Log.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Para_Tree.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Global.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/inlinedFunct.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/logFunct.hpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Data_LB_Interface.tpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Local_Tree_3D.tpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Global_2D.tpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Para_Tree_3D.tpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Map.tpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Data_Comm_Interface.tpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Comm_Buffer.tpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Para_Tree_2D.tpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Octant_3D.tpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Intersection_2D.tpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Octant_2D.tpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Intersection_3D.tpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Local_Tree_2D.tpp"
    "/home/haysam/Dev/BitP_Mesh/src/pablo/Class_Global_3D.tpp"
    )
endif()


# UseBITPIT.cmake
# -----------
#
# This file sets up include directories, link directories, and
# compiler settings for a project to use BITPIT.  It should not be
# included directly, but rather through the BITPIT_USE_FILE setting
# obtained from BITPITConfig.cmake.

if(BITPIT_USE_FILE_INCLUDED)
  return()
endif()
set(BITPIT_USE_FILE_INCLUDED 1)

# Update CMAKE_MODULE_PATH so includes work.
list(APPEND CMAKE_MODULE_PATH ${BITPIT_CMAKE_DIR})

# Setp up variables needed by external dependencies
list(FIND BITPIT_EXTERNAL_DEPENDENCIES "PETSc" _PETSc_index)
if (${_PETSc_index} GREATER -1)
    if (NOT PETSC_DIR AND DEFINED ENV{PETSC_DIR})
        set(DEFAULT_PETSC_DIR "$ENV{PETSC_DIR}")
    else()
        set(DEFAULT_PETSC_DIR "")
    endif()
    set(PETSC_DIR "${DEFAULT_PETSC_DIR}" CACHE PATH "Installation directory of PETSC library")

    if (NOT PETSC_ARCH AND DEFINED ENV{PETSC_ARCH})
        set(DEFAULT_PETSC_ARCH "$ENV{PETSC_ARCH}")
    else()
        set(DEFAULT_PETSC_ARCH "")
    endif()
    set(PETSC_ARCH "${DEFAULT_PETSC_ARCH}" CACHE STRING "Build architecture")
endif()
unset(_PETSc_index)

# Add compiler flags needed to use BITPIT.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${BITPIT_REQUIRED_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BITPIT_REQUIRED_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${BITPIT_REQUIRED_EXE_LINKER_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${BITPIT_REQUIRED_SHARED_LINKER_FLAGS}")
set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${BITPIT_REQUIRED_MODULE_LINKER_FLAGS}")

# Add preprocessor definitions needed to use BITPIT.
set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${BITPIT_DEFINITIONS})

# Add include directories needed to use BITPIT.
include_directories(${BITPIT_INCLUDE_DIRS})

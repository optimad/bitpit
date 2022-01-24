# Specify the version being used as well as the language
cmake_minimum_required(VERSION 3.10)

# CBLAS
add_subdirectory(CBLAS)

# LAPACKE
add_subdirectory(LAPACKE)

# PETSc
add_subdirectory(PETSc)

# METIS
add_subdirectory(METIS)
#CMakeNuGetTools
if(MSVC)
    add_subdirectory(CMakeNuGetTools)
endif()

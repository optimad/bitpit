# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# CBLAS
add_subdirectory(CBLAS)

# LAPACKE
add_subdirectory(LAPACKE)

# PETSc
add_subdirectory(PETSc)

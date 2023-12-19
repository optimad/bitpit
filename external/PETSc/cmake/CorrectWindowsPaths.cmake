# CorrectWindowsPaths - this module defines one macro
#
# CONVERT_CYGWIN_PATH( PATH )
#  This uses the command cygpath (provided by cygwin or msys2) to convert
#  unix-style paths into paths useable by cmake on windows
#
# BEWARE WIN32 enviroments will need the variable PETSC_CYGPATH_EXECUTABLE 
# (containing the path to cygpath.exe command coming with Cygwin or Msys2)  
# to enable path conversion from PETSc Unix paths to win compliant paths. 

macro (CONVERT_CYGWIN_PATH _path)
  if (PETSC_CYGPATH_EXECUTABLE)
    EXECUTE_PROCESS(COMMAND ${PETSC_CYGPATH_EXECUTABLE} -m ${${_path}}
      OUTPUT_VARIABLE ${_path})
      string (STRIP ${${_path}} ${_path})
  endif (PETSC_CYGPATH_EXECUTABLE)
endmacro (CONVERT_CYGWIN_PATH)


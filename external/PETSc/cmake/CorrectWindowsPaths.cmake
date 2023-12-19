# CorrectWindowsPaths - this module defines one macro
#
# CONVERT_CYGWIN_PATH( PATH )
#  This uses the command cygpath (provided by cygwin) to convert
#  unix-style paths into paths useable by cmake on windows
#
# The variable PETSC_CYGPATH_EXECUTABLE can be used to specify the path
# of cygpath.

macro (CONVERT_CYGWIN_PATH _path)
  if (WIN32)
    set(CYGPATH_EXECUTABLE "cygpath.exe")
    if (PETSC_CYGPATH_DIR)
      set(CYGPATH_EXECUTABLE "${PETSC_CYGPATH_DIR}/${CYGPATH_EXECUTABLE}")
    endif (PETSC_CYGPATH_DIR)

    EXECUTE_PROCESS(COMMAND ${CYGPATH_EXECUTABLE} -m ${${_path}}
      OUTPUT_VARIABLE ${_path})
    string (STRIP ${${_path}} ${_path})
  endif (WIN32)
endmacro (CONVERT_CYGWIN_PATH)


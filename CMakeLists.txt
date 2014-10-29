#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here
project(OPTIMAD_UCartMesh)


SET(TargetLib ucartmesh)
SET(TargetTest test_ucartmesh)

SET(DepLibs opt_base opt_IO  )

SET(BLDIR build/lib/ )
SET(IDIR include )

EXEC_PROGRAM(pwd ../../ OUTPUT_VARIABLE OptLibsDir)

ADD_DEFINITIONS(-std=c++0x)

include_directories( ${IDIR} )
FOREACH( lib ${DepLibs} )

    set(INAME "${OptLibsDir}/${lib}/${IDIR}")
    set(LNAME "${OptLibsDir}/${lib}/${BLDIR}")

    message( ${INAME} )
    include_directories( ${INAME} )
    link_directories( ${LNAME} )

ENDFOREACH( lib ${DepLibs})

ADD_SUBDIRECTORY(lib)
ADD_SUBDIRECTORY(test_lib )

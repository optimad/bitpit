SET(TargetLib opt_LA)
SET(TargetTest test_LA)

SET(DepLibs opt_base opt_IO )

SET(BLDIR build/lib/ )
SET(IDIR include )

EXEC_PROGRAM(pwd ../src/base OUTPUT_VARIABLE BaseLibDir)
EXEC_PROGRAM(pwd ../src/IO OUTPUT_VARIABLE IOLibDir)
EXEC_PROGRAM(pwd ../../ OUTPUT_VARIABLE OptLibsDir)

ADD_DEFINITIONS(-std=c++0x)

include_directories( ${BaseLibDir} )
include_directories( ${IOLibDir} )
include_directories( ${IDIR} )
FOREACH( lib ${DepLibs} )

    set(INAME "${OptLibsDir}/${lib}/${IDIR}")
    set(LNAME "${OptLibsDir}/${lib}/${BLDIR}")

    message( ${INAME} )
    include_directories( ${INAME} )
    link_directories( ${LNAME} )

ENDFOREACH( lib ${DepLibs})

ADD_SUBDIRECTORY(test_lib )

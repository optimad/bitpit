SET(TargetLib opt_IO)
SET(TargetTest test_IO)

SET(DepLibs opt_base )

SET(BLDIR build/lib/ )

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message( "building debug")
    SET(BLDIR debug/lib/ )
#    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -g -O0")
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

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

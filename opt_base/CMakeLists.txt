SET(TargetLib opt_base)
SET(TargetTest test_ob)

SET(DepLibs)

SET(BLDIR build/lib/ )

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message( "building debug")
    SET(BLDIR debug/lib/ )
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

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

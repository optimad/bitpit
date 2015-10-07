#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here

include_directories( ${IOLibDir} )

file(GLOB test_src "*.cpp")
add_executable( ${TargetTest} ${test_src})

target_link_libraries(${TargetTest} ${BASE_LIBRARY})
TARGET_LINK_LIBRARIES(${TargetTest} ${IO_LIBRARY} )

FOREACH( lib ${DepLibs} )

    set(LNAME "${OptLibsDir}/${lib}/${BLDIR}/lib${lib}.a")

    if(EXISTS  ${LNAME})
        TARGET_LINK_LIBRARIES(${TargetTest} ${lib} )
    endif(EXISTS  ${LNAME})

ENDFOREACH( lib ${DepLibs})

#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here

file(GLOB lib_src "*.cpp")

add_library( ${TargetLib} STATIC ${lib_src})

FOREACH( lib ${DepLibs} )

    set(LNAME "${OptLibsDir}/${lib}/${BLDIR}/lib${lib}.a")

    if(EXISTS  ${LNAME})
        TARGET_LINK_LIBRARIES(${TargetLib} ${lib} )
    endif(EXISTS  ${LNAME})

ENDFOREACH( lib ${DepLibs})

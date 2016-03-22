set (ARMATUS_MAIN Armatus.cpp ArmatusUtil.cpp ArmatusParams.cpp ArmatusDAG.cpp IntervalScheduling.cpp)

include_directories( 
${Boost_INCLUDE_DIR}
${TOP_SOURCE_DIR}/include
)

link_directories(${Boost_LIBRARY_DIR})

add_executable(armatus ${ARMATUS_MAIN})
target_link_libraries(armatus ${Boost_LIBRARIES} z)

install (TARGETS armatus DESTINATION bin)

set (ARMATUS_MAIN Armatus.cpp )

include_directories( 
${TOP_SOURCE_DIR}/include
)

link_directories()

add_executable(armatus ${ARMATUS_MAIN})
target_link_libraries(armatus)

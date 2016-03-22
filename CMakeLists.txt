cmake_minimum_required(VERSION 2.8)

project(ARMATUS)

## Prefer static to dynamic libraries
## We have to build boost on TravisCI each time for this, so avoid it for now
#SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

set(CMAKE_CXX_FLAGS "-O3 -g -std=c++11 -Wall")

set (USING_GCC FALSE)
set (USING_CLANG FALSE)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    set (USING_GCC TRUE)
    set (USING_CLANG FALSE)
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set (USING_CLANG TRUE)
    set (USING_GCC FALSE)
endif ()

if (USING_GCC)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
endif ()


set(TOP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_VERBOSE_MAKEFILE true)

find_package(Boost COMPONENTS program_options iostreams graph system)

add_subdirectory(src)


# build a CPack driven installer package
# include (InstallRequiredSystemLibraries)
# set (CPACK_RESOURCE_FILE_LICENSE  
#      "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
# set (CPACK_PACKAGE_VERSION_MAJOR "${Tutorial_VERSION_MAJOR}")
# set (CPACK_PACKAGE_VERSION_MINOR "${Tutorial_VERSION_MINOR}")
# set (CPACK_GENERATOR "TGZ")
# include (CPack)

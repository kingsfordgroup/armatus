cmake_minimum_required(VERSION 2.8)

project(ARMATUS)

set(CMAKE_CXX_FLAGS "-O3 -g -std=c++11 -Wall")

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

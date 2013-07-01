cmake_minimum_required(VERSION 2.8)

project(ARMATUS)

set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ -Wall")

set(TOP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_VERBOSE_MAKEFILE true)

find_package(Boost COMPONENTS program_options iostreams graph system)

add_subdirectory(src)

cmake_minimum_required(VERSION 3.10)

project(inq VERSION 0.1)
enable_language(Fortran)

configure_file(inq_config.h.in inq_config.h)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_COMPILER mpic++)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(FortranCInterface)


FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")

#find_package(Boost REQUIRED COMPONENTS boost_serialization)

# Required for inq_config.h
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(external_libs)
add_subdirectory(src)






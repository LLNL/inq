cmake_minimum_required(VERSION 3.10)
include(GNUInstallDirs)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

project(inq VERSION 0.1)


if(ENABLE_CUDA)
  enable_language(CUDA)
  set(SRC_EXT cu)
else()
  set(SRC_EXT cpp)
endif()

enable_language(CXX)
enable_language(Fortran)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_COMPILER mpic++)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")

find_package(BLAS REQUIRED)
link_libraries(${BLAS_LIBRARIES})

find_package(LAPACK REQUIRED)
link_libraries(${LAPACK_LIBRARIES})

find_package(Libxc REQUIRED)
link_libraries(${Libxc_LIBRARIES})

find_package(Boost REQUIRED COMPONENTS serialization)
link_libraries(${Boost_LIBRARIES})

find_package(FFTW REQUIRED COMPONENTS DOUBLE_LIB)
link_libraries(${FFTW_LIBRARIES})

# Required for inq_config.h
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_BINARY_DIR}/external_libs/pseudopod)
configure_file(inq_config.h.in inq_config.h)

install(DIRECTORY share/ DESTINATION share/inq)

add_subdirectory(external_libs)
add_subdirectory(src)
add_subdirectory(tests)

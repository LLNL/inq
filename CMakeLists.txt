cmake_minimum_required(VERSION 3.11)
include(GNUInstallDirs)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

project(inq VERSION 0.1)

set(CMAKE_VERBOSE_MAKEFILE ON)

enable_testing()
include(CTest)

if(ENABLE_CUDA)
  set(SRC_EXT cu)
  enable_language(CUDA)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe \"--diag_suppress=implicit_return_from_non_void_function\"") #nvcc
  find_package(CUDA)
  set(INQ_CUDA_LIBRARIES ${CUDA_CUFFT_LIBRARIES} ${CUDA_cusolver_LIBRARY} ${CUDA_CUBLAS_LIBRARIES})
  link_libraries(${INQ_CUDA_LIBRARIES})
else()
  set(SRC_EXT cpp)
endif()

link_libraries(m)

add_subdirectory(external_libs)

enable_language(CXX)
enable_language(Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")

#FFTW has to go before blas to avoid unscrupulous (i.e. MKL) blas implementations that include FFTW and don't implement it properly
find_package(FFTW REQUIRED COMPONENTS DOUBLE_LIB)
include_directories(${FFTW_INCLUDE_DIRS})
link_libraries(${FFTW_LIBRARIES})

find_package(Heffte)
#get_target_property(HEFFTE_INCLUDE_DIRS Heffte::Heffte INTERFACE_INCLUDE_DIRECTORIES)
if(Heffte_FOUND)
  link_libraries(Heffte::Heffte)
endif()

find_package(BLAS REQUIRED)
link_libraries(${BLAS_LIBRARIES})

find_package(LAPACK REQUIRED)
link_libraries(${LAPACK_LIBRARIES})

set(Libxc_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/external_libs/libxc/src/ ${PROJECT_BINARY_DIR}/external_libs/libxc ${PROJECT_BINARY_DIR}/external_libs/libxc/gen_funcidx)
set(Libxc_LIBRARIES ${PROJECT_BINARY_DIR}/external_libs/libxc/libxc.a)
include_directories(${Libxc_INCLUDE_DIRS})
link_libraries(${Libxc_LIBRARIES})
link_libraries(xc)

find_package(Boost REQUIRED COMPONENTS serialization)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

# Required for inq_config.h
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_BINARY_DIR}/external_libs/pseudopod)
configure_file(inq_config.h.in inq_config.h)

#Generate inc++
set(INCXX scripts/inc++)
configure_file(${INCXX}.bash.in ${INCXX} @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${INCXX} DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(DIRECTORY share/ DESTINATION share/inq)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(examples)


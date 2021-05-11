cmake_minimum_required(VERSION 3.11) # for reference: ubuntu 20.04 has cmake 3.16
include(GNUInstallDirs)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

project(inq VERSION 0.1)

set(CMAKE_VERBOSE_MAKEFILE ON)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

if(POLICY CMP0104)
  cmake_policy(SET CMP0104 NEW)
endif()

enable_testing()
include(CTest)

if(ENABLE_CUDA)
  set(SRC_EXT cu)
  enable_language(CUDA)
  find_package(CUDA 11.0.0 REQUIRED) # deprecated in 3.17 replaced by FindCudaToolkit
  cuda_select_nvcc_arch_flags(ARCH_FLAGS ${CMAKE_CUDA_ARCHITECTURES}) # needs cmake 3.7, may autodetect at configure time
  list(APPEND CUDA_NVCC_FLAGS ${ARCH_FLAGS})
  string(REPLACE "-gencode;" "--generate-code=" ARCH_FLAGS "${ARCH_FLAGS}")
  set(CMAKE_CUDA_FLAGS "--forward-unknown-to-host-linker ${ARCH_FLAGS} -std=c++17 --expt-relaxed-constexpr --extended-lambda --Werror=cross-execution-space-call -Xcudafe \"--diag_suppress=implicit_return_from_non_void_function\"")
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

find_package(MPI REQUIRED)
link_libraries(MPI::MPI_CXX) # needs cmake 3.9
message(STATUS "MPI Executable for running programs:" ${MPIEXEC_EXECUTABLE})
message(STATUS "MPI number of processors detected on the host system: " ${MPIEXEC_MAX_NUMPROCS})

#FFTW has to go before blas to avoid unscrupulous (i.e. MKL) blas implementations that include FFTW and don't implement it properly
find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
include_directories(PkgConfig::FFTW)
link_libraries(PkgConfig::FFTW)

if(ENABLE_HEFFTE)
  find_package(Heffte)
  #get_target_property(HEFFTE_INCLUDE_DIRS Heffte::Heffte INTERFACE_INCLUDE_DIRECTORIES)
  link_libraries(Heffte::Heffte)
endif()

find_package(BLAS REQUIRED)
link_libraries(${BLAS_LIBRARIES})

set(BLA_VENDOR Generic) # ESSL doesn't implement dgelss_, now needed in the code
find_package(LAPACK REQUIRED)
link_libraries(${LAPACK_LIBRARIES})

set(Libxc_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/external_libs/libxc/src/ ${PROJECT_BINARY_DIR}/external_libs/libxc ${PROJECT_BINARY_DIR}/external_libs/libxc/gen_funcidx)
set(Libxc_LIBRARIES ${PROJECT_BINARY_DIR}/external_libs/libxc/libxc.a)
include_directories(${Libxc_INCLUDE_DIRS})
link_libraries(${Libxc_LIBRARIES})
link_libraries(xc)

find_package(Boost REQUIRED COMPONENTS serialization)
message(STATUS "Boost actual version found " ${Boost_VERSION})
if(Boost_VERSION LESS_EQUAL 105300) # may be higher than 1.53.0
	message(STATUS "Applying option BOOST_NO_AUTO_PTR to patch for old Boost")
	add_compile_options(-DBOOST_NO_AUTO_PTR)
endif()
include_directories(Boost::serialization) #${Boost_INCLUDE_DIRS})
link_libraries(Boost::serialization)

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


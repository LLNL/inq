# -*- indent-tabs-mode: t -*-
cmake_minimum_required(VERSION 3.18) # for reference: ubuntu 20.04 has cmake 3.16
message("CMake version ${CMAKE_VERSION}")
include(GNUInstallDirs)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

project(inq VERSION 0.1)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(FETCHCONTENT_UPDATES_DISCONNECTED ON CACHE BOOL "")

message(STATUS "downloading Catch2 v3p3")
Include(FetchContent)
FetchContent_Declare(
	Catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG        7cf2f88e50f0d1de324489c31db0314188423b6d
	)
FetchContent_MakeAvailable(Catch2)
include_directories(${catch2_SOURCE_DIR}/src)

enable_testing()
include(CTest)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(ENABLE_NCCL FALSE)
if(ENABLE_CUDA)
  set(SRC_EXT cu)
  enable_language(CUDA)
  find_package(CUDA 11.0.0 REQUIRED) # deprecated in 3.17 replaced by FindCudaToolkit

#  nvcc 11 doesn't work with gcc 10, use option `--pass-thru -DCMAKE_CUDA_HOST_COMPILER=g++-9`") # https://gitlab.com/npneq/inq/-/issues/107

  # temporary workaround to detect CUDA arch, for between the deprecation of FindCUDA in CMake 3.10 and the availability of FindCudaToolkit in CMake 3.17 # https://stackoverflow.com/a/68223399/225186
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
	  message("-DCMAKE_CUDA_ARCHITECTURES is not set, automatic architecture detection...")
      include(FindCUDA/select_compute_arch)
      CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU_CCS_1)
      string(STRIP "${INSTALLED_GPU_CCS_1}" INSTALLED_GPU_CCS_2)
      string(REPLACE " " ";" INSTALLED_GPU_CCS_3 "${INSTALLED_GPU_CCS_2}")
      string(REPLACE "." "" CUDA_ARCH_LIST "${INSTALLED_GPU_CCS_3}")
      SET(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_LIST})
      set_property(GLOBAL PROPERTY CUDA_ARCHITECTURES "${CUDA_ARCH_LIST}")
      set_property(GLOBAL PROPERTY CUDA_STANDARD 17)  # with cmake 3.13 this is not enough, -std=c++17 is needed below
	  message("setting CUDA_ARCHITECTURES global property to ${CUDA_ARCH_LIST}")
  else()
	  message("-DCMAKE_CUDA_ARCHITECTURES is set...")
	  set(CUDA_ARCH_LIST "${CMAKE_CUDA_ARCHITECTURES}")
	  message("setting CUDA_ARCH_LIST to ${CMAKE_CUDA_ARCHITECTURES}")
      set_property(GLOBAL PROPERTY CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
      set_property(GLOBAL PROPERTY CUDA_STANDARD 17)  # with cmake 3.13 this is not enough, -std=c++17 is needed below
	  message("setting CUDA_ARCHITECTURES global property to ${CMAKE_CUDA_ARCHITECTURES}")
	  set(ARCH_FLAGS "--generate-code=arch=compute_${CMAKE_CUDA_ARCHITECTURES},code=[compute_${CMAKE_CUDA_ARCHITECTURES},sm_${CMAKE_CUDA_ARCHITECTURES}]")
  endif()

  message("... at the end ARCH_FLAGS was set to ${ARCH_FLAGS}")

  string(APPEND CMAKE_CUDA_FLAGS " --forward-unknown-to-host-linker ${ARCH_FLAGS} -std=c++17 --expt-relaxed-constexpr --extended-lambda --Werror=cross-execution-space-call -Xcudafe \"--diag_suppress=implicit_return_from_non_void_function\"")

  if(${CUDA_CUFFT_LIBRARIES} STREQUAL "CUDA_cufft_LIBRARY-NOTFOUND")
    include_directories("${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/include")
    link_directories   ("${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/lib64")
    link_libraries("-lcufft")
  else()
    link_libraries(${CUDA_CUFFT_LIBRARIES})
  endif()
  if(${CUDA_cusolver_LIBRARY} STREQUAL "CUDA_cusolver_LIBRARY-NOTFOUND")
    include_directories("${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/include")
    link_directories   ("${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/lib64")
    link_libraries("-lcusolver")
  else()
    link_libraries(${CUDA_cusolver_LIBRARY})
  endif()
  if(${CUDA_cublas_LIBRARY} STREQUAL "CUDA_cublas_LIBRARY-NOTFOUND")
    include_directories("${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/include")
    link_directories   ("${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/lib64")
    link_libraries("-lcublas")
  else()
    link_libraries(${CUDA_cublas_LIBRARY})
  endif()
	find_package(NCCL)

	set(ENABLE_NCCL ${NCCL_FOUND})
	if(${ENABLE_NCCL})
		include_directories(${NCCL_INCLUDE_DIRS})
		link_libraries(${NCCL_LIBRARIES})
	endif()
	
else()
  set(SRC_EXT cpp)
endif()

link_libraries(m)

add_subdirectory(external_libs)

enable_language(Fortran)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")

find_package(MPI REQUIRED COMPONENTS CXX)
include_directories(${MPI_CXX_INCLUDE_PATH})
link_libraries(MPI::MPI_CXX) # needs cmake 3.9
message(STATUS "MPI Executable for running programs:" ${MPIEXEC_EXECUTABLE})
message(STATUS "MPI number of processors detected on the host system: " ${MPIEXEC_MAX_NUMPROCS})
include(CMakePrintHelpers)
cmake_print_properties(TARGETS MPI::MPI_CXX PROPERTIES INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES)

#Check for MPI_Isendrecv_replace that sometimes MPI doesn't have
set(CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_CXX)
include(CheckCXXSymbolExists)
check_cxx_symbol_exists(MPI_Isendrecv_replace mpi.h HAVE_MPI_ISENDRECV_REPLACE)

#FFTW has to go before blas to avoid unscrupulous (i.e. MKL) blas implementations that include FFTW and don't implement it properly
find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
include_directories(PkgConfig::FFTW)
link_libraries(PkgConfig::FFTW)
cmake_print_properties(TARGETS PkgConfig::FFTW PROPERTIES INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES)

if(ENABLE_HEFFTE)
  find_package(Heffte REQUIRED)
  #get_target_property(HEFFTE_INCLUDE_DIRS Heffte::Heffte INTERFACE_INCLUDE_DIRECTORIES)
  link_libraries(Heffte::Heffte)
endif()

find_package(BLAS REQUIRED)
link_libraries(${BLAS_LIBRARIES})

set(BLA_VENDOR Generic) # ESSL doesn't implement dgelss_, now needed in the code
find_package(LAPACK REQUIRED)
link_libraries(${LAPACK_LIBRARIES})

#LIBXC
set(Libxc_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/external_libs/libxc/src/ ${PROJECT_BINARY_DIR}/external_libs/libxc ${PROJECT_BINARY_DIR}/external_libs/libxc/gen_funcidx)
set(Libxc_LIBRARIES ${PROJECT_BINARY_DIR}/external_libs/libxc/libxc.a)
include_directories(${Libxc_INCLUDE_DIRS})
link_libraries(${Libxc_LIBRARIES})
link_libraries(xc)

# Multi
set(multi_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/external_libs/multi/include)
include_directories(${PROJECT_SOURCE_DIR}/external_libs/multi/include)

########################
# SPGLIB
########################

message(STATUS "Downloading spglib")

set(CMAKE_BUILD_TYPE_BACKUP ${MAKE_BUILD_TYPE})
set(CMAKE_BUILD_TYPE Release)

FetchContent_Declare(
  spglib
	GIT_REPOSITORY https://gitlab.com/npneq/spglib.git
	GIT_TAG        f04dadd3175a60279f86780630045d3b3f64b8aa
	)
FetchContent_MakeAvailable(spglib)

set(CMAKE_BUILD_TYPE ${MAKE_BUILD_TYPE_BACKUP})

set(spglib_INCLUDE_DIRS ${spglib_SOURCE_DIR}/src/)
set(spglib_LIBRARIES ${spglib_BINARY_DIR}/libsymspg.a)
include_directories(${spglib_INCLUDE_DIRS})
link_libraries(symspg_static)

########################
# BOOST
########################

find_package(Boost REQUIRED COMPONENTS serialization filesystem system)
message(STATUS "Boost actual version found " ${Boost_VERSION})
if(Boost_VERSION LESS_EQUAL 105300) # may be higher than 1.53.0
	message(STATUS "Applying option BOOST_NO_AUTO_PTR to patch for old Boost")
	add_compile_options(-DBOOST_NO_AUTO_PTR)
endif()
include_directories(Boost::serialization Boost::filesystem Boost::system) #${Boost_INCLUDE_DIRS})
link_libraries(Boost::serialization Boost::filesystem Boost::system)

# Required for inq_config.h
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_BINARY_DIR}/external_libs/pseudopod)
configure_file(inq_config.h.in inq_config.h)

#Generate inc++
set(INCXX scripts/inc++)
if(ENABLE_CUDA)
set(INCPP_COMPILER ${CMAKE_CUDA_COMPILER})
set(INCPP_FLAGS ${CMAKE_CUDA_FLAGS})
else()
set(INCPP_COMPILER ${CMAKE_CXX_COMPILER})
set(INCPP_FLAGS ${CMAKE_CXX_FLAGS})
endif()
configure_file(${INCXX}.bash.in ${INCXX} @ONLY)
configure_file(scripts/run_test.sh.in scripts/run_test.sh @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${INCXX} DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(DIRECTORY share/ DESTINATION share/inq)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(speed_tests)
add_subdirectory(examples)
add_subdirectory(benchmarks)



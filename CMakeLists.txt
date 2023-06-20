# -*- indent-tabs-mode: t -*-
# Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.21)
message("CMake version ${CMAKE_VERSION}")
include(GNUInstallDirs)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

project(inq VERSION 0.1)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(FETCHCONTENT_UPDATES_DISCONNECTED ON CACHE BOOL "")
set(FETCHCONTENT_QUIET FALSE)
Include(FetchContent)

add_library(inq INTERFACE)
target_include_directories(inq
    INTERFACE
     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
     $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
     $<INSTALL_INTERFACE:src>
)

target_compile_options(inq INTERFACE
  $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:
  --expt-relaxed-constexpr --extended-lambda --Werror=cross-execution-space-call -Xcudafe=--display_error_number -Xcudafe=--diag_error=incompatible_assignment_operands -Xcudafe=--diag_error=returning_ptr_to_local_variable -Xcudafe=--diag_error=subscript_out_of_range -Xcudafe=--diag_error=used_before_set -Xcudafe=--diag_error=undefined_preproc_id -Xcudafe=--diag_error=implicit_func_decl -Xcudafe=--diag_error=implicit_return_from_non_void_function -Xcudafe=--diag_error=missing_type_specifier
  >
)

enable_testing()
include(CTest)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(ENABLE_CUDA)
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

    target_include_directories(inq INTERFACE "${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/include")
    target_link_directories   (inq INTERFACE "${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/lib64")
    target_link_libraries     (inq INTERFACE "-lcufft")
  else()
    link_libraries(${CUDA_CUFFT_LIBRARIES})
    target_link_libraries(inq INTERFACE ${CUDA_CUFFT_LIBRARIES})
  endif()
  if(${CUDA_cusolver_LIBRARY} STREQUAL "CUDA_cusolver_LIBRARY-NOTFOUND")
    include_directories("${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/include")
    link_directories   ("${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/lib64")
    link_libraries("-lcusolver")

    target_include_directories(inq INTERFACE "${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/include")
    target_link_directories   (inq INTERFACE "${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/lib64")
    target_link_libraries     (inq INTERFACE "-lcusolver")
  else()
    link_libraries(${CUDA_cusolver_LIBRARY})
    target_link_libraries(inq INTERFACE ${CUDA_cusolver_LIBRARY})
  endif()
  if(${CUDA_cublas_LIBRARY} STREQUAL "CUDA_cublas_LIBRARY-NOTFOUND")
    include_directories("${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/include")
    link_directories   ("${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/lib64")
    link_libraries("-lcublas")

    target_include_directories(inq INTERFACE "${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/include")
    target_link_directories   (inq INTERFACE "${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/lib64")
    target_link_libraries     (inq INTERFACE "-lcublas")
  else()
    link_libraries(${CUDA_cublas_LIBRARY})
    target_link_libraries(inq INTERFACE ${CUDA_cublas_LIBRARY})
  endif()

	if(ENABLE_NCCL)
		find_package(NCCL REQUIRED)
		include_directories(${NCCL_INCLUDE_DIRS})
		link_libraries(${NCCL_LIBRARIES})

		target_include_directories(inq INTERFACE ${NCCL_INCLUDE_DIRS})
		target_link_libraries(inq INTERFACE ${NCCL_LIBRARIES})
  endif()
endif()

link_libraries(m)

enable_language(Fortran)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")

find_package(MPI REQUIRED COMPONENTS CXX)
include_directories(${MPI_CXX_INCLUDE_PATH})
link_libraries(MPI::MPI_CXX) # needs cmake 3.9

target_link_libraries(inq INTERFACE MPI::MPI_CXX)

message(STATUS "MPI Executable for running programs:" ${MPIEXEC_EXECUTABLE})
message(STATUS "MPI number of processors detected on the host system: " ${MPIEXEC_MAX_NUMPROCS})
include(CMakePrintHelpers)
cmake_print_properties(TARGETS MPI::MPI_CXX PROPERTIES INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES)

#Check for MPI_Isendrecv_replace that sometimes MPI doesn't have
set(CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_CXX)
include(CheckCXXSymbolExists)
check_cxx_symbol_exists(MPI_Isendrecv_replace mpi.h HAVE_MPI_ISENDRECV_REPLACE)

#FFTW has to go before blas to avoid unscrupulous (i.e. MKL) implementations that include FFTW and don't implement the full FFTW interface
find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW fftw3 IMPORTED_TARGET)
if(FFTW_FOUND)
  link_libraries(PkgConfig::FFTW)
  #include_directories(PkgConfig::FFTW)  # commented since it always generates a spurious -I...PkgConfig::FFTW option
  cmake_print_properties(TARGETS PkgConfig::FFTW PROPERTIES INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES)
  target_link_libraries(inq INTERFACE PkgConfig::FFTW)
else()  # some systems do not have fftw3 registered in pkg-config, in this case it falls back to good old -lfftw3 which can fail later
  if(NOT DEFINED FFTW_LIBRARIES)
    message(WARNING "FFTW not detected by pkg-config, using fftw3 or define -DFFTW_LIBRARIES (e.g. =-lfftw3) and -DFFTW_INCLUDE_DIRS (e.g. =/usr/include)")
    link_libraries(fftw3)
    include_directories(${FFTW_INCLUDE_DIRS})  # commented since it always generates a spurious -I...PkgConfig::FFTW option
    target_link_libraries(inq INTERFACE fftw3)
    target_include_directories(inq INTERFACE ${FFTW_INCLUDE_DIRS})
  else()
    link_libraries(${FFTW_LIBRARIES})
    include_directories(${FFTW_INCLUDE_DIRS})  # commented since it always generates a spurious -I...PkgConfig::FFTW option
    target_link_libraries(inq INTERFACE ${FFTW_LIBRARIES})
    target_include_directories(inq INTERFACE ${FFTW_INCLUDE_DIRS})
  endif()
  check_cxx_symbol_exists(fftw_plan_dft_1d fftw3.h FFTW3_CHECKED)
  if(FFTW_FOUND)
  else()
    message(FATAL_ERROR "FFTW not found")
  endif()
endif()

if(ENABLE_HEFFTE)
  find_package(Heffte REQUIRED)
  #get_target_property(HEFFTE_INCLUDE_DIRS Heffte::Heffte INTERFACE_INCLUDE_DIRECTORIES)
  link_libraries(Heffte::Heffte)
  target_link_libraries(inq INTERFACE Heffte::Heffte)
endif()

find_package(BLAS REQUIRED)
link_libraries(${BLAS_LIBRARIES})

target_link_libraries(inq INTERFACE ${BLAS_LIBRARIES})

# ESSL doesn't implement dgelss_, needed in the code, may need to override set(BLA_VENDOR Generic)
find_package(LAPACK REQUIRED)
link_libraries(${LAPACK_LIBRARIES})

target_link_libraries(inq INTERFACE ${LAPACK_LIBRARIES})

#GPU::RUN
add_subdirectory(${PROJECT_SOURCE_DIR}/external_libs/gpurun)
set(gpurun_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/external_libs/gpurun/include)
target_include_directories(inq INTERFACE ${PROJECT_SOURCE_DIR}/external_libs/gpurun/include)

#PSEUDOPOD
add_subdirectory(external_libs/pseudopod)
include_directories(${PROJECT_BINARY_DIR}/external_libs/pseudopod)
target_include_directories(inq INTERFACE ${PROJECT_BINARY_DIR}/external_libs/pseudopod)

#LIBXC
message("\nCONFIGURING LIBXC\n")
set(DISABLE_FXC ON)
set(BUILD_TESTING OFF)
add_subdirectory(external_libs/libxc EXCLUDE_FROM_ALL)
message("\nDONE CONFIGURING LIBXC\n")
set(Libxc_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/external_libs/libxc/src/ ${PROJECT_BINARY_DIR}/external_libs/libxc ${PROJECT_BINARY_DIR}/external_libs/libxc/gen_funcidx)
include_directories(${Libxc_INCLUDE_DIRS})

target_include_directories(inq INTERFACE ${Libxc_INCLUDE_DIRS})
target_link_libraries(inq INTERFACE ${Libxc_LIBRARIES} xc)

# Multi
set(multi_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/external_libs/multi/include)
include_directories(${PROJECT_SOURCE_DIR}/external_libs/multi/include)

add_subdirectory(${PROJECT_SOURCE_DIR}/external_libs/multi)
target_link_libraries(inq INTERFACE multi)

# MPI3
set(mpi3_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/external_libs/mpi3/include)
include_directories(${PROJECT_SOURCE_DIR}/external_libs/mpi3/include)

add_subdirectory(${PROJECT_SOURCE_DIR}/external_libs/mpi3)
target_link_libraries(inq INTERFACE mpi3)

########################
# SPDLOG
########################

message(STATUS "Downloading spdlog")

FetchContent_Declare(
  spdlog
	GIT_REPOSITORY https://github.com/gabime/spdlog.git
	GIT_TAG        v1.11.0
  GIT_PROGRESS   TRUE
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(spdlog)

set(spdlog_INCLUDE_DIRS ${spdlog_SOURCE_DIR}/include/)
include_directories(${spdlog_INCLUDE_DIRS})
target_include_directories(inq INTERFACE ${spdlog_INCLUDE_DIRS})
target_link_libraries(inq INTERFACE spdlog)

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
  GIT_PROGRESS   TRUE
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(spglib)

set_property(DIRECTORY ${spglib_SOURCE_DIR} PROPERTY EXCLUDE_FROM_ALL YES)

set(CMAKE_BUILD_TYPE ${MAKE_BUILD_TYPE_BACKUP})

set(spglib_INCLUDE_DIRS ${spglib_SOURCE_DIR}/src/)
include_directories(${spglib_INCLUDE_DIRS})
link_libraries(symspg_static)

target_include_directories(inq
  INTERFACE
  ${spglib_INCLUDE_DIRS}
)
target_link_libraries(inq INTERFACE symspg_static)
target_link_libraries(inq
   INTERFACE
   ${spglib_LIBRARIES}
)

########################
# GEMMI
########################

target_include_directories(inq INTERFACE ${PROJECT_SOURCE_DIR}/external_libs/gemmi/include)
include_directories(${PROJECT_SOURCE_DIR}/external_libs/gemmi/include)

########################
# BOOST
########################

find_package(Boost REQUIRED COMPONENTS serialization filesystem system)
message(STATUS "Boost actual version found " ${Boost_VERSION})
if(Boost_VERSION LESS_EQUAL "1.53.0")  # 105300
	message(STATUS "Applying option BOOST_NO_AUTO_PTR to patch for old Boost")
	add_compile_options(-DBOOST_NO_AUTO_PTR)
endif()
include_directories(Boost::serialization Boost::filesystem Boost::system) #${Boost_INCLUDE_DIRS})
link_libraries(Boost::serialization Boost::filesystem Boost::system)

target_include_directories(inq INTERFACE Boost::serialization Boost::filesystem Boost::system) #${Boost_INCLUDE_DIRS})
target_link_libraries(inq INTERFACE Boost::serialization Boost::filesystem Boost::system)

#######################
# GENERATE inq_config.h
######################

include_directories(${PROJECT_BINARY_DIR})

target_include_directories(inq INTERFACE ${PROJECT_BINARY_DIR})
target_include_directories(inq INTERFACE ${PROJECT_SOURCE_DIR}/external_libs)

configure_file(inq_config.h.in inq_config.h)
configure_file(scripts/run_test.sh.in scripts/run_test.sh @ONLY)

install(DIRECTORY share/ DESTINATION share/inq)

# this makes this library CMake's FetchContent friendly https://www.foonathan.net/2022/06/cmake-fetchcontent/
if(NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  return()
endif()

message(STATUS "Downloading Catch2 v3p3")
FetchContent_Declare(
	Catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG        7cf2f88e50f0d1de324489c31db0314188423b6d
  GIT_PROGRESS   TRUE
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(Catch2)
include_directories(${catch2_SOURCE_DIR}/src)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(speed_tests)
add_subdirectory(examples)
add_subdirectory(benchmarks)

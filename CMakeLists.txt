cmake_minimum_required(VERSION 3.11) # for reference: ubuntu 20.04 has cmake 3.16
message("CMake version ${CMAKE_VERSION}")
include(GNUInstallDirs)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

project(inq VERSION 0.1)

set(CMAKE_VERBOSE_MAKEFILE ON)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

#if(POLICY CMP0104)
#  cmake_policy(SET CMP0104 NEW)
#endif()

enable_testing()
include(CTest)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(ENABLE_CUDA)
  set(SRC_EXT cu)
  enable_language(CUDA)
  find_package(CUDA 11.0.0 REQUIRED) # deprecated in 3.17 replaced by FindCudaToolkit
#  if(NOT DEFINED ENV{CUDA_ARCH_OVERRIDE})
#    cuda_select_nvcc_arch_flags(ARCH_FLAGS ${CMAKE_CUDA_ARCHITECTURES}) # needs cmake 3.7, may autodetect at configure time, in local system reboot after updating drivers
#  endif()
#  list(APPEND CUDA_NVCC_FLAGS ${ARCH_FLAGS})
#  string(REPLACE "-gencode;" "--generate-code=" ARCH_FLAGS "${ARCH_FLAGS}")

#  nvcc 11 doesn't work with gcc 10, use option `--pass-thru -DCMAKE_CUDA_HOST_COMPILER=g++-9`") # https://gitlab.com/npneq/inq/-/issues/107

  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    # temporary workaround to detect CUDA arch, for between the deprecation of FindCUDA in CMake 3.10 and the availability of FindCudaToolkit in CMake 3.17 # https://stackoverflow.com/a/68223399/225186
    include(FindCUDA/select_compute_arch)
    CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU_CCS_1)
    string(STRIP "${INSTALLED_GPU_CCS_1}" INSTALLED_GPU_CCS_2)
    string(REPLACE " " ";" INSTALLED_GPU_CCS_3 "${INSTALLED_GPU_CCS_2}")
    string(REPLACE "." "" CUDA_ARCH_LIST "${INSTALLED_GPU_CCS_3}")
    SET(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_LIST})
    set_property(GLOBAL PROPERTY CUDA_ARCHITECTURES "${CUDA_ARCH_LIST}")
    set_property(GLOBAL PROPERTY CUDA_STANDARD 17)
  endif()

  string(APPEND CMAKE_CUDA_FLAGS " --forward-unknown-to-host-linker ${ARCH_FLAGS} --expt-relaxed-constexpr --extended-lambda --Werror=cross-execution-space-call -Xcudafe \"--diag_suppress=implicit_return_from_non_void_function\"")
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
#  set(INQ_CUDA_LIBRARIES ${CUDA_CUFFT_LIBRARIES} ${CUDA_cusolver_LIBRARY} ${CUDA_CUBLAS_LIBRARIES})
#  link_libraries(${INQ_CUDA_LIBRARIES})
else()
  set(SRC_EXT cpp)
endif()

link_libraries(m)

add_subdirectory(external_libs)

enable_language(Fortran)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
link_libraries(MPI::MPI_CXX) # needs cmake 3.9
message(STATUS "MPI Executable for running programs:" ${MPIEXEC_EXECUTABLE})
message(STATUS "MPI number of processors detected on the host system: " ${MPIEXEC_MAX_NUMPROCS})
include(CMakePrintHelpers)
cmake_print_properties(TARGETS MPI::MPI_CXX PROPERTIES INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES)

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

#SPGLIB

message("\nCONFIGURING SPGLIB\n")
set(CMAKE_BUILD_TYPE_BACKUP ${MAKE_BUILD_TYPE})
set(CMAKE_BUILD_TYPE Release)
add_subdirectory(external_libs/spglib)
set(CMAKE_BUILD_TYPE ${MAKE_BUILD_TYPE_BACKUP})

message("\nDONE CONFIGURING SPGLIB\n")

set(spglib_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/external_libs/spglib/src/)
set(spglib_LIBRARIES ${PROJECT_BINARY_DIR}/external_libs/spglib/libsymspg.a)
include_directories(${spglib_INCLUDE_DIRS})
link_libraries(symspg_static)


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

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${INCXX} DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(DIRECTORY share/ DESTINATION share/inq)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(examples)


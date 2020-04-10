#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT(inq, alpha, xavier@llnl.gov)
AM_INIT_AUTOMAKE([subdir-objects tar-ustar])
AC_CONFIG_SRCDIR([src/hamiltonian/ks_hamiltonian.hpp])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

# enable cuda first, since this is required for the compiler selection
AC_ARG_ENABLE(cuda, AS_HELP_STRING([--enable-cuda], [Cuda support (experimental)]))

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX(mpic++)
AX_CXX_COMPILE_STDCXX([14], [], [])

#we may need to set the linker for the cuda case, so pass that to automake
if test x"$CXXLD" == x ; then
   CXXLD=$CXX
fi
AC_SUBST(CXXLD)

AC_PROG_CC(cc)
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

AC_SYS_LARGEFILE

# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_FSEEKO
AC_FUNC_STRTOD
AC_CHECK_FUNCS([gettimeofday memset mkdir pow sqrt strspn uname])

if test x"${enable_cuda}" == x"yes" ; then
  AC_MSG_NOTICE([
  ================================================================================
    CHECKING FOR CUDA
  ================================================================================])

  AC_DEFINE(HAVE_CUDA, 1, [defined if cuda support is enabled])
  ACX_CUDA
fi

AC_MSG_NOTICE([
================================================================================
  CHECKING FOR FFT LIBRARIES
================================================================================])

found_fft=no

#FFTW3
if test x$found_fft != xyes ; then
  ACX_FFTW3
fi

if test x$acx_fftw3_ok == xyes ; then
  found_fft=yes
  fft_lib=fftw3
  CFLAGS="$CFLAGS $CFLAGS_FFTW3"
  CXXFLAGS="$CXXFLAGS $CFLAGS_FFTW3"
  LIBS="$LIBS $LIBS_FFTW3"
fi

AC_MSG_NOTICE([
================================================================================
  CHECKING FOR FORTRAN LIBRARIES
================================================================================])
dnl We need to check for fortran, to link with blas and other libraries
AC_PROG_FC()
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS
		
AC_MSG_NOTICE([
================================================================================
  CHECKING FOR REQUIRED LINEAR ALGEBRA LIBRARIES (BLAS AND LAPACK)
================================================================================])
dnl check for BLAS
ACX_BLAS([], AC_MSG_ERROR([could not find required BLAS library]))

dnl check for LAPACK
ACX_LAPACK([], AC_MSG_ERROR([could not find required LAPACK library]))

AC_MSG_NOTICE([
================================================================================
  CHECKING FOR REQUIRED LIBXC
================================================================================])
AC_LANG_PUSH([C++])

ACX_LIBXC

AC_LANG_POP

AC_MSG_NOTICE([
================================================================================
  CHECKING FOR MPI
================================================================================])
AC_LANG_PUSH([C++])

AC_CHECK_HEADER([mpi.h], [], [AC_MSG_ERROR([could not find MPI])])

AC_LANG_POP

AC_MSG_NOTICE([
================================================================================
  CONFIGURATION DONE, GENERATING MAKEFILES
================================================================================])

AC_CONFIG_SUBDIRS([external_libs/pseudopod])

AC_CONFIG_FILES([Makefile
                 external_libs/Makefile
		             external_libs/catch2/Makefile
                 external_libs/rapidxml/Makefile
                 external_libs/tinyformat/Makefile
                 external_libs/pcg-cpp/Makefile
                 share/Makefile
                 share/unit_tests_data/Makefile
                 src/Makefile])

AC_OUTPUT

AC_MSG_NOTICE([

================================================================================
  CONFIGURATION SUMMARY
================================================================================

C compiler       :  $CC $CFLAGS
C++ compiler     :  $CXX $CXXFLAGS
Fortran compiler :  $FC $FCFLAGS
FFT library      :  $fft_lib
BLAS             :  $LIBS_BLAS
LAPACK           :  $LIBS_LAPACK

])




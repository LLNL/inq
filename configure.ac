#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT(inq, alpha, xavier@llnl.gov)
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_SRCDIR([src/pseudo/element.hpp])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX(mpic++)
AC_PROG_CC(mpicc)
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

AC_MSG_NOTICE([
================================================================================
  CHECKING FOR FORTRAN LIBRARIES
================================================================================])
dnl We need to check for fortran, to link with blas and other libraries
AC_PROG_FC(mpif90)
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS

AC_MSG_NOTICE([
================================================================================
  CONFIGURATION DONE, GENERATING MAKEFILES
================================================================================])
AC_CONFIG_FILES([Makefile
                 external_libs/Makefile
		 external_libs/catch2/Makefile
 		 external_libs/rapidxml/Makefile
		 share/Makefile
		 share/pseudopotentials/Makefile
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
])




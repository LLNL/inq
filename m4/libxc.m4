## Copyright (C) 2010-2015 M. Marques, X. Andrade, D. Strubbe, M. Oliveira
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
##

AC_DEFUN([ACX_LIBXC], [
acx_libxc_ok=no
acx_libxc_v3=no
acx_libxc_v4=no

dnl Check if the library was given in the command line
dnl if not, use environment variables or defaults
AC_ARG_WITH(libxc-prefix, [AS_HELP_STRING([--with-libxc-prefix=DIR], [Directory where libxc was installed.])])

# Set FCFLAGS_LIBXC only if not set from environment
if test x"$CXXFLAGS_LIBXC" = x; then
  case $with_libxc_prefix in
    "") CXXFLAGS_LIBXC="-I/usr/include" ;;
    *)  CXXFLAGS_LIBXC="-I$with_libxc_prefix/include" ;;
  esac
fi

AC_ARG_WITH(libxc-include, [AS_HELP_STRING([--with-libxc-include=DIR], [Directory where libxc headers were installed.])])
case $with_libxc_include in
  "") ;;
  *)  FCFLAGS_LIBXC="-I$with_libxc_include" ;;
esac

dnl Backup LIBS and FCFLAGS
acx_libxc_save_LIBS="$LIBS"
acx_libxc_save_CXXFLAGS="$CXXFLAGS"

CXXFLAGS="$CXXFLAGS_LIBXC $acx_libxc_save_CXXFLAGS"

AC_CHECK_LIB([xc], [xc_lda_exc])

acx_libxc_ok=$ac_cv_lib_xc_xc_lda_exc

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_libxc_ok" = xyes; then
  AC_DEFINE(HAVE_LIBXC, 1, [Defined if you have the LIBXC library.])
else
  AC_MSG_ERROR([Could not find required libxc library ( >= v 2.0.0).])
fi

AC_SUBST(CXXFLAGS_LIBXC)
AC_SUBST(LIBS_LIBXC)
CXXFLAGS="$acx_libxc_save_CXXFLAGS"
LIBS="$acx_libxc_save_LIBS"
])dnl ACX_LIBXC

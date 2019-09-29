## Copyright (C) 2016 X. Andrade
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

AC_DEFUN([ACX_CUDA],
[

  acx_cuda_ok=yes

  dnl BACKUP LIBS AND CFLAGS
  acx_cuda_save_LIBS="$LIBS"
  acx_cuda_save_LDFLAGS="$LDFLAGS"
  acx_cuda_save_CFLAGS="$CFLAGS"

  dnl Check if the library was given in the command line
  AC_ARG_WITH(cuda-prefix, [AS_HELP_STRING([--with-cuda-prefix=DIR], [Directory where cuda was installed.])])

  # Set CFLAGS_CUDA only if not set from environment
  if test x"$CFLAGS_CUDA" = x; then
    case $with_cuda_prefix in
      "") CFLAGS_CUDA="-I/usr/include" ;;
      *)  CFLAGS_CUDA="-I$with_cuda_prefix/include" ;;
    esac
  fi

  CFLAGS="$CFLAGS_CUDA $acx_cuda_save_CFLAGS"

  if test ! -z "$with_cuda_prefix"; then
    LDFLAGS="-L$with_cuda_prefix/lib64"
  else
    LDFLAGS=""
  fi

  LIBS=""

  AC_CHECK_LIB(cuda, cuInit, [], [acx_cuda_ok=no], [$acx_cuda_save_LIBS])
  AC_CHECK_LIB(nvrtc, nvrtcCreateProgram, [], [acx_cuda_ok=no], [$acx_cuda_save_LIBS])
  AC_CHECK_LIB(cublas, cublasCreate_v2, [], [acx_cuda_ok=no], [$acx_cuda_save_LIBS])
  AC_CHECK_LIB(cufft, cufftPlan3d, [], [acx_cuda_ok=no], [$acx_cuda_save_LIBS])

  LIBS_CUDA="$LDFLAGS $LIBS"

  AC_MSG_CHECKING([for whether we have all the required cuda libraries])

  if test x$acx_cuda_ok != xyes; then

    AC_MSG_RESULT([no])
    AC_MSG_ERROR([Could not find the cuda library])

  else
      AC_MSG_RESULT([yes ($CFLAGS_CUDA $LIBS_CUDA)])
  fi

  AC_SUBST(CFLAGS_CUDA)
  AC_SUBST(LIBS_CUDA)

  CFLAGS="$acx_cuda_save_CFLAGS"
  LDFLAGS="$acx_cuda_save_LDFLAGS"
  LIBS="$acx_cuda_save_LIBS"

])

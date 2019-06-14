# Configure path for the GNU Scientific Library
# Christopher R. Gabriel <cgabriel@linux.it>, April 2000
# extensively rewritten, D Strubbe 2012

AC_DEFUN([AX_PATH_GSL],
[
  AC_ARG_WITH(gsl-prefix,[AS_HELP_STRING([--with-gsl-prefix=DIR], [Prefix where GSL is installed])],
    gsl_prefix="$withval", gsl_prefix="")
  AC_ARG_WITH(gsl-exec-prefix,[AS_HELP_STRING([--with-gsl-exec-prefix=DIR], [Exec prefix where GSL is installed])],
    gsl_exec_prefix="$withval", gsl_exec_prefix="")

  if test "x${GSL_CONFIG+set}" != xset ; then
    if test "x$gsl_prefix" != x ; then
      GSL_CONFIG="$gsl_prefix/bin/gsl-config"
    fi
    if test "x$gsl_exec_prefix" != x ; then
      GSL_CONFIG="$gsl_exec_prefix/bin/gsl-config"
    fi
  fi

  AC_PATH_PROG(GSL_CONFIG, gsl-config, no)
  min_gsl_major_version=1
  min_gsl_minor_version=9

  if test "$GSL_CONFIG" = "no" ; then
    no_gsl=yes
    AC_MSG_ERROR([The gsl-config script installed by GSL could not be found.
    If GSL was installed in PREFIX, make sure PREFIX/bin is in
    your path, or set the GSL_CONFIG environment variable to the
    full path to gsl-config.])
  else
    AC_MSG_CHECKING(for GSL - version >= $min_gsl_major_version.$min_gsl_minor_version)
    no_gsl=""
    GSL_CFLAGS=`$GSL_CONFIG --cflags`
    GSL_LIBS=`$GSL_CONFIG --libs`

    gsl_major_version=`$GSL_CONFIG --version | \
      sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${gsl_major_version}" = "x" ; then
      gsl_major_version=0
    fi

    gsl_minor_version=`$GSL_CONFIG --version | \
      sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${gsl_minor_version}" = "x" ; then
      gsl_minor_version=0
    fi

    gsl_micro_version=`$GSL_CONFIG --version | \
      sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${gsl_micro_version}" = "x" ; then
      gsl_micro_version=0
    fi

    AC_MSG_RESULT($gsl_major_version.$gsl_minor_version.$gsl_micro_version)
    if test $gsl_major_version -lt $min_gsl_major_version -o \
            $gsl_major_version -eq $min_gsl_major_version -a $gsl_minor_version -lt $min_gsl_minor_version; then
      AC_MSG_ERROR([A newer version of GSL is required.])
    fi

    AC_MSG_CHECKING(whether GSL can be linked)
    ac_save_CFLAGS="$CFLAGS"
    ac_save_LIBS="$LIBS"
    CFLAGS="$CFLAGS $GSL_CFLAGS"
    LIBS="$LIBS $GSL_LIBS"

    AC_LINK_IFELSE([AC_LANG_PROGRAM([
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
],[
  gsl_spline x;
  gsl_asinh(1.0);
])],
      AC_MSG_RESULT(yes),[
      AC_MSG_RESULT(no)
      no_gsl=yes
      AC_MSG_ERROR([The test program failed to link. See the file config.log for the
      exact error that occured. This usually means GSL was incorrectly installed
      or that you have moved GSL since it was installed. In the latter case, you
      can edit the gsl-config script: $GSL_CONFIG.])])
  fi

  CFLAGS="$ac_save_CFLAGS"
  LIBS="$ac_save_LIBS"

  AC_SUBST(GSL_CFLAGS)
  AC_SUBST(GSL_LIBS)
])

AU_ALIAS([AM_PATH_GSL], [AX_PATH_GSL])

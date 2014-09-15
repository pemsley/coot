# SYNOPSIS
#
#   SINGLE_FFTW2
#
# DESCRIPTION
#
#   Test for the single-precision version of FFTW2 libraries (fftw and rfftw).
#   FFTW version 2 can be compiled in either single or double precision
#   (default is double, option --enable-float is for single).
#   The single-precision version can be compiled with or without prefix
#   (option --enable-type-prefix adds 's' to both headers and libraries).
#   This macro first checks for prefixed version and then for unprefixed,
#   checking if single precision is used.
#
#   If the library is found FFTW2_LIBS is set (with AC_SUBST).
#   Otherwise scripts stops with an error message.
#   FFTW2_PREFIX_S is defined (AC_DEFINE) if prefix is used. Used it for
#   including headers:
#     #ifdef FFTW2_PREFIX_S
#     # include <srfftw.h>
#     #else
#     # include <rfftw.h>
#     #endif
#
#   Example:
#     AC_SEARCH_LIBS(cos, m, , AC_MSG_ERROR([math library not found.]))
#     SINGLE_FFTW2
#     LIBS="$FFTW2_LIBS $LIBS"
#
# LICENSE
#
#   Public Domain
#

AC_DEFUN([SINGLE_FFTW2],
[
saved_LIBS="$LIBS"
AC_LANG_PUSH(C++)

AC_MSG_CHECKING([for prefixed single-precision FFTW2 (sfftw.h)])
FFTW2_LIBS="-lsrfftw -lsfftw"
LIBS="$FFTW2_LIBS $saved_LIBS"
AC_TRY_LINK([#include <sfftw.h>],
            [float a; fftw_real *p = &a; return *fftw_version], 
            have_fftw=yes, have_fftw=no)
AC_MSG_RESULT($have_fftw)
if test $have_fftw = yes; then
  AC_DEFINE(FFTW2_PREFIX_S, 1, [Define if FFTW2 is prefixed.])
else
  AC_MSG_CHECKING([for not prefixed single-precision FFTW2 (fftw.h)])
  FFTW2_LIBS="-lrfftw -lfftw"
  LIBS="$FFTW2_LIBS $saved_LIBS"
  AC_TRY_LINK([#include <fftw.h>],
              [float a; fftw_real *p = &a; return *fftw_version], 
              [AC_MSG_RESULT(yes)],
              [AC_MSG_RESULT(no)
               AC_MSG_ERROR([single-precision FFTW 2 library not found.])])
fi

AC_LANG_POP(C++)
LIBS="$saved_LIBS"
AC_SUBST(FFTW2_LIBS)
])

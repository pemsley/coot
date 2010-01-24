# simple libpng configure based on GTK+ configure

dnl AM_PATH_GTK_2_0([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND [, MODULES]]]])
dnl Test for GTK+, and define GTK_CFLAGS and GTK_LIBS, if gthread is specified in MODULES, 
dnl pass to pkg-config
dnl
AC_DEFUN([AM_PATH_LIBPNG],[
LIBPNG_CFLAGS="-DUSE_LIBPNG"
LIBPNG_LIBS="-lpng"

# BL: workaround needed for new MinGW
ac_cv_build_alias=${ac_cv_build_alias:=$build_alias}

case $ac_cv_build_alias in 

   *-mingw*)

dnl 
dnl Get the cflags and libraries from pkg-config (only for windows)
dnl
  AC_MSG_CHECKING([for LIBPNG])

     AC_PATH_PROG(PKG_CONFIG, pkg-config, no)
     pkg_config_args=libpng
     LIBPNG_CFLAGS="-DUSE_LIBPNG "`$PKG_CONFIG $pkg_config_args --cflags`
     # need to use ../png.h since we need zlib.h too
     LIBPNG_CFLAGS=$(echo $LIBPNG_CFLAGS | sed 's/\/libpng12//')
     LIBPNG_LIBS=`$PKG_CONFIG $pkg_config_args --libs`

      ac_save_CFLAGS="$CFLAGS"
      ac_save_LIBS="$LIBS"
      CFLAGS="$CFLAGS $LIBPNG_CFLAGS"
      LIBS="$LIBPNG_LIBS $LIBS"
    dnl
    dnl Now check if the installed LIBPNG is sufficiently new. (Also sanity
    dnl checks the results of pkg-config to some extent)
    dnl actually only 'check' if its there
    dnl
      rm -f conf.libpngtest
      AC_TRY_RUN([
#include <png.h>
#include <stdio.h>
#include <stdlib.h>

int 
main ()
{
  int major, minor, micro;
  char *tmp_version;

  system ("touch conf.libpngtest");

  return 1;
}
])
     rm -f conf.libpngtest
     AC_MSG_RESULT([yes])
     CFLAGS="$ac_save_CFLAGS"
     LIBS="$ac_save_LIBS"
     break;;
  esac

  AC_SUBST(LIBPNG_CFLAGS)
  AC_SUBST(LIBPNG_LIBS)
])

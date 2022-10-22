AC_DEFUN([AM_PATH_GLOB],
[
AC_PROVIDE([AM_PATH_GLOB])

AC_ARG_WITH(glob-prefix, [  --with-glob-prefix=PFX Prefix where glob can been found],
 glob_prefix="$withval",
 glob_prefix="")

saved_LIBS="$LIBS"
saved_CFLAGS="$CFLAGS"

if test x$glob_prefix != x; then
 ac_GLOB_CFLAGS="-I$glob_prefix/include"
 ac_GLOB_LIBS="-L$glob_prefix/lib -lglob"
else
 ac_GLOB_CFLAGS=""
 ac_GLOB_LIBS="-lglob"
fi

AC_MSG_CHECKING([for glob (in windows)])

       LIBS="$save_LIBS $ac_GLOB_LIBS"
       CFLAGS="$save_CXXFLAGS $ac_GLOB_CFLAGS"
       #
       AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <glob.h> ]], [[ ]])],[have_generic_glob=yes],[have_generic_glob=no])
       AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <glob.h> ]], [[ ]])],[have_glob=yes],[have_glob=no])
       AC_MSG_RESULT($have_glob)

if test x$have_glob = xyes; then

  GLOB_CFLAGS="$ac_GLOB_CFLAGS"
  GLOB_LIBS="$ac_GLOB_LIBS"

else

  echo Oops couldnt find glob.h compilation will fail!

fi

LIBS="$saved_LIBS"
CFLAGS="$saved_CFLAGS"

AC_SUBST(GLOB_CFLAGS)
AC_SUBST(GLOB_LIBS)

])

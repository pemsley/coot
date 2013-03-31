 
# AM_PATH_GTKCANVAS([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
AC_DEFUN([AM_PATH_GTKCANVAS],
[
AC_PROVIDE([AM_PATH_GTKCANVAS])

AC_ARG_WITH(gtkcanvas-prefix, [  --with-gtkcanvas-prefix=PFX Prefix where GTKCANVAS has been built],
 gtkcanvas_prefix="$withval",
 gtkcanvas_prefix="")


saved_LIBS="$LIBS"
saved_CFLAGS="$CFLAGS"

if test x$gtkcanvas_prefix != x; then


#   # first test for imlib, setting IMLIB_CFLAGS and IMLIB_LIBS
# case $ac_cv_build_alias in 
#
#  *-mingw*)
#    break;;
#  *)
#    AM_PATH_IMLIB
#    ;;
# esac

  # echo debug: here we have IMLIB_CFLAGS: $IMLIB_CFLAGS and IMLIB_LIBS $IMLIB_LIBS   

	# the majority of cases, we will try to configure with:
	# --gtkcanvas-prefix=/some/thing
	#

 GTKCANVAS_CFLAGS="-DHAVE_GTK_CANVAS -I$gtkcanvas_prefix/include -I$gtkcanvas_prefix/include/gtk-canvas -I$gtkcanvas_prefix/include/gnome-1.0 $IMLIB_CFLAGS"
 #
 # Similarly for gtkcanvas, the uninstalled library position is simply in
 # $gtkcanvas_prefix, but the installed is in the standard prefixed subdirectory.
 #
 # SGI compiler CC (CXX=CC) needs -lm to link maths library, but 
 # GCC c++ does not.
 #
 GTKCANVAS_LDOPTS="-L$gtkcanvas_prefix/lib -lgtk-canvas $IMLIB_LIBS -lart_lgpl -lm"
else
 # the compiler looks in the "standard" places for GTKCANVAS.  
 GTKCANVAS_CFLAGS="-DHAVE_GTK_CANVAS"
 GTKCANVAS_LDOPTS="-lgtk-canvas -lart_lgpl"
fi

AC_MSG_CHECKING([for GtkCanvas])

	LIBS="$save_LIBS $GTKCANVAS_LDOPTS $GTK_LIBS"
	CFLAGS="$save_CFLAGS $GTKCANVAS_CFLAGS $GTK_CFLAGS"
	#
	# AC_TRY_LINK uses the c compiler (set by AC_LANG), so we will
	# temporarily reassign $CC to the c++ compiler.
 	#
	AC_TRY_LINK([#include "gtk-canvas.h"] ,[ GtkCanvas *a;  ], have_gtkcanvas=yes, have_gtkcanvas=no)
	AC_MSG_RESULT($have_gtkcanvas)

if test x$have_gtkcanvas = xyes; then

 LIBS="$saved_LIBS"
 CFLAGS="$saved_CFLAGS"
 GTKCANVAS_CFLAGS="$GTKCANVAS_CFLAGS"
 GTKCANVAS_LIBS="$GTKCANVAS_LDOPTS"
ifelse([$1], , :, [$1])

else

 # we want to fail if a prefix was given but gtkcanvas was not found
 # have_gtkcanvas=no here
 # 
 if test x$gtkcanvas_prefix != x ; then
	# fail
	AC_MSG_ERROR([--with-gtkcanvas-prefix was specified, but no gtkcanvas was found])
 fi

 LIBS="$saved_LIBS"
 CFLAGS="$saved_CFLAGS"
 GTKCANVAS_LIBS=""
 GTKCANVAS_CFLAGS=""
 ifelse([$2], , :, [$2])

fi

AC_SUBST(GTKCANVAS_CFLAGS)
AC_SUBST(GTKCANVAS_LIBS)

])

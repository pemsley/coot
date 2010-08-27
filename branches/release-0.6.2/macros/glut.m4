# glut.m4
# 
# Copyright 2003 The University of York
# Author: Paul Emsley
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

# AM_PATH_GLUT([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
AC_DEFUN([AM_PATH_GLUT],
[
AC_REQUIRE([AM_PATH_GTK])
AC_PROVIDE([AM_PATH_GLUT])


AC_ARG_WITH(
	glut-prefix, [  --with-glut-prefix=PFX for GLUT Installation],
	glut_prefix="$withval",
	glut_prefix="")
	

AC_MSG_CHECKING([for GLUT])

saved_LIBS="$LIBS"
saved_CFLAGS="$CFLAGS"

if test x$glut_prefix != x ; then
        # the prefix was specified....
	# SGIs need -lX11, -lXmu and -lm to link glut apps, so we use
	# the GL and GTK libs to find them... sigh.
	GLUT_CFLAGS="-I$glut_prefix/include $GL_CFLAGS"
	# vanilla UNIX:
	GLUT_LDOPTS="-L$glut_prefix/$acl_libdirstem -lglut -L/usr/X11R6/$acl_libdirstem -lXmu $GL_LIBS"

        # but there is no libXmu on Darwin... argh!

        # special case:
	# BL: workaround needed for new MinGW
	ac_cv_build_alias=${ac_cv_build_alias:=$build_alias}
	
	case $ac_cv_build_alias in

	powerpc-apple-darwin*)
	   GLUT_LDOPTS="-L/usr/X11R6/$acl_libdirstem -L$glut_prefix/$acl_libdirstem -lglut $GL_LIBS"
	   GLUT_CFLAGS="-I/usr/X11R6/include $GLUT_CFLAGS"
	   break;;

	i686-pc-linux-gnu*)
	   if test -e /etc/fedora-release ; then
	      GLUT_LDOPTS="-L$glut_prefix/$acl_libdirstem -lglut $GL_LIBS -L/usr/X11R6/$acl_libdirstem -lXi -lXmu"
	   fi
	   break;;

        i686-apple-darwin*) 
	   # Mactel, same as G4
	   GLUT_LDOPTS="-L/usr/X11R6/$acl_libdirstem -L$glut_prefix/$acl_libdirstem -lglut $GL_LIBS"
	   GLUT_CFLAGS="-I/usr/X11R6/include $GLUT_CFLAGS"
	   break;;

	i686-pc-linux-gnu*)
	   if test -e /etc/fedora-release ; then
	      GLUT_LDOPTS="-L$glut_prefix/$acl_libdirstem -lglut $GL_LIBS -L/usr/X11R6/lib -lXi -lXmu"
	   fi
	   break;;
	
	i686-pc-cygwin*)
 	   # This is linking for freeglut.  with conventional glut, we
 	   # don't pass an argument and it links OK - it seems.
	   GLUT_LDOPTS="-L$glut_prefix/lib -lglut -L/usr/X11R6/lib $GL_LIBS -lXext -lX11 -lm"
	   break;;

	*-mingw*)
	   # BL says:: this is for MinGW Windows
	   # for now we try just with GL_LIBS may have to be GTK_LIBS
	   # we have freeglut on windows now, so test for this		
	   if test -e $glut_prefix/GL/include/freeglut.h ; then
	      GLUT_LDOPTS="-L$glut_prefix/lib -lfreeglut $GL_LIBS"
	   else
	      GLUT_LDOPTS="-L$glut_prefix/lib -lglut32 $GL_LIBS"
	   fi
	   break;;

	mips-sgi-irix6*)
	   GLUT_LDOPTS="-L$glut_prefix/lib -lglut $GL_LIBS -lXmu -lX11 -lm"
	   break;;

       x86_64-*-linux-gnu)
	   # e.g. our 64-bit Suse 9.x machine (hkl101) [no Xmu]
  	   GLUT_LDOPTS="-L$glut_prefix/lib -lglut -L/usr/X11R6/lib $GL_LIBS"
	   break;;

	esac

        # as above, but on Ubuntu 64-bit, we dont have Xmu (and we dont
        # have ac_cv_build_alias either!)
	#
        case "$ac_cv_build" in 
           x86_64-unknown-linux-gnu)
  	   GLUT_LDOPTS="-L$glut_prefix/lib -lglut -L/usr/X11R6/lib $GL_LIBS"
	   break;;
        esac


dnl Recall that on cygwin we need 
dnl GLUT_LDOPTS="-L$glut_prefix/lib" -lglut -L/usr/X11R6/lib $GL_LIBS -lXxf86vm -lXext -lX11" 
dnl it seems to me that -L/usr/X11R6/lib should be part of GL_LIBS, it is not clear to 
dnl me why it isn't. 

dnl	# now tinker for the special case additions:
dnl Actually, I no longer think that this is a special case
dnl redhat 9 needs it too... 
dnl	case $ac_cv_build_alias in

dnl	mips*)
dnl      	   GLUT_LDOPTS="$GLUT_LDOPTS -lXmu $GTK_LIBS"
dnl	   break ;;
dnl	esac

else

	GLUT_CFLAGS=""
	GLUT_LDOPTS="-lglut $GL_LIBS -lXmu $GTK_LIBS"

        # special case:
        case $ac_cv_build_alias in
	powerpc-apple-darwin*)
	   GLUT_LDOPTS="-L$glut_prefix/lib -lglut $GL_LIBS"
	   break;;

	i686-pc-linux-gnu*)
	   if test -e /etc/fedora-release ; then
	      GLUT_LDOPTS="-lglut $GL_LIBS -L/usr/X11R6/lib -lXi -lXmu"
	   fi
	   break;;
	
	*-mingw*)
	   # BL says:: this is for MinGW Windows
	   # for now we try just with GL_LIBS may have to be GTK_LIBS
	   GLUT_LDOPTS="-lglut32 $GL_LIBS $GTK_LIBS"
	   break;;

	mips-sgi-irix6*)
	   GLUT_LDOPTS="-lglut $GL_LIBS -lXmu -lX11 -lm"
	   break;

	esac

fi


	LIBS="$saved_LIBS $GLUT_LDOPTS"
	CFLAGS="$saved_CFLAGS $GLUT_CFLAGS"

	AC_TRY_LINK([#include <GL/glut.h>] , [glutWireCube(1);], have_glut=yes, have_glut=no)
	
	AC_MSG_RESULT($have_glut)

if test x$have_glut = xyes; then

	LIBS="$saved_LIBS"
	CFLAGS="$saved_CFLAGS"
	GLUT_CFLAGS="$GLUT_CFLAGS"
	GLUT_LIBS="$GLUT_LDOPTS"
	ifelse([$1], , :, [$1])

else
	LIBS="$saved_LIBS"
 	CFLAGS="$saved_CFLAGS"
	GLUT_LIBS=""
	GLUT_CFLAGS=""
	ifelse([$2], , :, [$2])
fi

AC_SUBST(GLUT_CFLAGS)
AC_SUBST(GLUT_LIBS)

])

# libglade.m4
# 
# Copyright 2007 The University of Oxford
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
 

AC_DEFUN([AM_PATH_LIBGLADE],
[
AC_PROVIDE([AM_PATH_LIBGLADE])

AC_ARG_WITH(libglade-prefix,
[  --with-libglade-prefix=LIBGLADE_PREFIX  Location of libglade],
LIBGLADE_PREFIX="$withval")

# echo xxxxxxxxxxxxxxxxxxxxxxxxxx LIBGLADE_PREFIX is $LIBGLADE_PREFIX

if TRUE = TRUE ; then

   AC_MSG_CHECKING(libglade)

   saved_LIBS="$LIBS"
   saved_CXXFLAGS="$CXXFLAGS"
   saved_CFLAGS="$CFLAGS"


   # Note:  I tried $LIBGLADE_PREFIX/lib/libglade-2.0.la, like Ezra 
   # says is the preferred method now, but libguile test in configure
   # fails to link if I do that:
   # /home/paule/glade-3/lib/libglade-2.0.la: file not recognized: File format not recognized
   # collect2: ld returned 1 exit status
   # 
   # so, go back to old way of -Lxx -lxx

   # 
   if test x$LIBGLADE_PREFIX != x ; then
      ac_LIBGLADE_LIBS="-L$LIBGLADE_PREFIX/lib -lglade-2.0"
      ac_LIBGLADE_CFLAGS="-I$LIBGLADE_PREFIX/include/libglade-2.0 -DUSE_LIBGLADE"

# Comment this out until I think what I was trying to do, currently it is 
# nonsense in the usual case (where --with-libglade is not specificied)
#
   else
      ac_LIBGLADE_LIBS=-lglade-2.0
      ac_LIBGLADE_CFLAGS="-I/usr/include/libglade-2.0"

   fi
   
   LIBS="$LIBS $ac_LIBGLADE_LIBS $pkg_cv_GTK_LIBS"
   CFLAGS="$CFLAGS $ac_LIBGLADE_CFLAGS $pkg_cv_GTK_CFLAGS -DUSE_LIBGLADE"
   AC_TRY_LINK([#include <glade/glade.h>], [GladeXML *xml = glade_xml_new("x",NULL,NULL);], have_libglade=yes, have_libglade=no)

   CXXFLAGS="$saved_CXXFLAGS"
   CFLAGS="$saved_CFLAGS"
   LIBS="$saved_LIBS"
   AC_MSG_RESULT($have_libglade)

   if test "$have_libglade" = yes ; then
      LIBGLADE_LIBS=$ac_LIBGLADE_LIBS
      LIBGLADE_CFLAGS=$ac_LIBGLADE_CFLAGS
   fi

fi

AC_SUBST(LIBGLADE_CFLAGS)
AC_SUBST(LIBGLADE_LIBS)

])

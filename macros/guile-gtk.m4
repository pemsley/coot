# guile-gtk.m4
# 
# Copyright 2006 The University of York
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


# The AM_PATH_GUILE is for where the libs are to be found.  This is for
# whether we should compile with Guile-Gtk, giving us compilation
# defines.

AC_DEFUN([AM_WITH_GUILE_GTK],
[AC_PROVIDE([AM_WITH_GUILE_GTK])

AC_ARG_WITH(guile-gtk, [  --with-guile-gtk build with guile-gtk to allow main gui extensions],
	with_guile_gtk="$withval",
	with_guile_gtk="")

if test x$with_guile_gtk != x ; then
   USE_GUILE_GTK=-DUSE_GUILE_GTK
else
   USE_GUILE_GTK=
fi

AC_SUBST(USE_GUILE_GTK)
])


# Note, use AM_WITH_GUILE_GTK first, this function depends on USE_GUILE_GTK
# 
# AM_PATH_GUILE_GTK([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
AC_DEFUN([AM_PATH_GUILE_GTK],
[
# AC_REQUIRE([AM_PATH_GTK])
AC_PROVIDE([AM_PATH_GUILE_GTK])


AC_ARG_WITH(
	guile-gtk-prefix, [  --with-guile-gtk-prefix=PRFX for guile-gtk Installation],
	guile_gtk_prefix="$withval",
	guile_gtk_prefix="")

AC_MSG_CHECKING([for guile-gtk])


saved_LIBS="$LIBS"
saved_CFLAGS="$CFLAGS"

if test x$guile_gtk_prefix != x; then
   # the prefix was specified....
   # now, were were GTK2 or GTK1? 
   # we have the variable coot_gtk2 in configure.in, which is either TRUE or FALSE
   # Let's use that here.  Perhaps there is a better (set by gtk macro?) variable?
   #
   # Why do I use the .la file here?
   # Let's revert to standard lib usage.
   # GUILE_GTK_LIBS=$guile_gtk_prefix/lib/libguilegtk-2.0.la
   #
   # mac os x doesn't have that file when guile-gtk is installed.
   GUILE_GTK_LIBS="-L$guile_gtk_prefix/lib -lguilegtk-2.0"
   # do we need to set this in fact?
   GUILE_GTK_CFLAGS="-I$guile_gtk_prefix/include"

else 
   if test -n "$USE_GUILE_GTK" ; then 
      # not sure if this code gets executed (much).
      GUILE_GTK_LIBS="-lguilegtk-2.0"
   fi
fi

AC_MSG_RESULT([$guile_gtk_prefix])
AC_SUBST(GUILE_GTK_LIBS)
AC_SUBST(GUILE_GTK_CFLAGS)

])

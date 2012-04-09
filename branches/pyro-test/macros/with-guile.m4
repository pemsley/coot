# with-guile.m4
# 
# Copyright 2004 The University of York
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



AC_DEFUN([AM_WITH_GUILE], 
[AC_PROVIDE([AM_USE_GUILE])

# guile, not with-guile
AC_ARG_WITH(guile, [  --with-guile=PFX Prefix where GUILE has been installed],
 with_guile_prefix="$withval",
 with_guile_prefix="")


saved_LIBS="$LIBS"
saved_CFLAGS="$CFLAGS"

if test x$with_guile = xyes; then
   echo Congratulations, you are using Guile
   coot_guile=true
else 
   echo Not using guile
   coot_guile=false
fi

if test x$coot_guile = xtrue ; then
   COOT_USE_GUILE=-DUSE_GUILE
else 
   COOT_USE_GUILE=""
fi

AC_SUBST(COOT_USE_GUILE)
])


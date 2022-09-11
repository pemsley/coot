# mmdb.m4
# 
# Copyright 2002 The University of York
# Copyright 2014 by Medical Research Council
# Author: Paul Emsley
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope qthat it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA


AC_DEFUN([AM_PATH_MMDB2],
[
AC_PROVIDE([AM_PATH_MMDB2])


AC_MSG_CHECKING([for mmdb2])

if ${PKG_CONFIG} mmdb2 ; then 
   MMDB_CXXFLAGS="$($PKG_CONFIG --cflags mmdb2)"
   MMDB_LIBS="$($PKG_CONFIG --libs mmdb2)"
   coot_found_mmdb2=yes
else
   coot_found_mmdb2=no
fi

AC_MSG_RESULT($coot_found_mmdb2)

if test $coot_found_mmdb2 = no ; then
   AC_MSG_FAILURE([mmdb2 not found])
else
   save_CXXFLAGS="$CXXFLAGS"
   CXXFLAGS="$CXXFLAGS $MMDB_CXXFLAGS"
   # do we have a verison of mmdb2 that has distances in links?
   AC_LANG_PUSH(C++)
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <mmdb2/mmdb_manager.h>]], [[mmdb::Link *l; mmdb::realtype ll=l->dist ]])],[have_link_distance=yes],[have_link_distance=no])
   if test $have_link_distance = yes ; then
      MMDB_CXXFLAGS="$MMDB_CXXFLAGS -DMMDB_HAS_LINK_DISTANCE"
   fi
   CXXFLAGS=$save_CXXFLAGS
   AC_LANG_POP(C++)
fi

AC_SUBST(MMDB_CXXFLAGS)
AC_SUBST(MMDB_LIBS)

])


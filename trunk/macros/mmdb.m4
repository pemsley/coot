# mmdb.m4
# 
# Copyright 2002 The University of York
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

# Note that mmdb as distributed by Eugene does not install and does
# conform to normal unix standard notions of the behaviour of include 
# and libraries.
#
# So we will depend on you doing something by hand, and that is to link
# (or copy) mmdb.a to libmmdb.a (and that it is in the top directory of 
# mmdb).
#

 
# AM_PATH_MMDB([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
# ----------------------------------------------------------
# set up for MMDB
#
AC_DEFUN([AM_PATH_MMDB],
[
AC_PROVIDE([AM_PATH_MMDB])

# you may want to change to mmdb (not mmdb-prefix) in future to be 
# consistent with clipper.
AC_ARG_WITH(mmdb-prefix,
  AC_HELP_STRING( [--with-mmdb-prefix=PFX], [Prefix where MMDB has been installed] ),
  [if test "x${withval}" != xyes ; then
   mmdb_prefix=$withval
   else
   mmdb_prefix=""    
   fi],
  [mmdb_prefix=""] )

saved_LIBS="$LIBS"
saved_CXXFLAGS="$CXXFLAGS"

if test x$mmdb_prefix != x; then
	# very likely the majority of cases, we will try to configure with:
	# --with-mmdb=/some/thing
	#

   if test -r "$mmdb_prefix/include/mmdb2/mmdb_manager.h"; then
      ac_MMDB_CXXFLAGS="-I$mmdb_prefix/include"
   fi
   #
   # SGI compiler CC (CXX=CC) needs -lm to link maths library, but 
   # GCC c++ does not.
   #
   # not libdirstem because mmdb is in lib no matter what.
   #  ac_MMDB_LDOPTS="-L$mmdb_prefix/$acl_libdirstem -L$mmdb_prefix -lmmdb -lm"
   ac_MMDB_LIBS="-L$mmdb_prefix/lib -lmmdb2 -lm"
else
   # the compiler looks in the "standard" places for MMDB.  In real life,
   # it would be quite unlikely that MMDB would be installed in /usr/include, 
   # /usr/lib etc. so this code will not usually find the right dependencies.
   ac_MMDB_CXXFLAGS=""
   ac_MMDB_LIBS="-lmmdb2 -lm"
fi

AC_MSG_CHECKING([for MMDB])

   LIBS="$save_LIBS $ac_MMDB_LIBS"
   CXXFLAGS="$save_CXXFLAGS $ac_MMDB_CXXFLAGS"
   #
   # AC_TRY_LINK uses the c compiler (set by AC_LANG), so we will
   # temporarily reassign $CC to the c++ compiler.
   #
   AC_LANG_PUSH(C++)
   AC_TRY_LINK([#include    <mmdb2/mmdb_manager.h>] ,[ mmdb::Manager a;  ], have_mmdb=yes, have_mmdb=no)
   AC_LANG_POP(C++)  # the language we have just quit
   AC_MSG_RESULT($have_mmdb)

if test x$have_mmdb = xyes; then

  LIBS="$saved_LIBS"
  HASH_FLAG=-DHAVE_MMDB_IGNORE_HASH
  CISPEP_FLAG=-DHAVE_MMDB_WITH_CISPEP
  CXXFLAGS="$saved_CXXFLAGS"
  MMDB_CXXFLAGS="$ac_MMDB_CXXFLAGS $HASH_FLAG $CISPEP_FLAG"
  MMDB_LIBS="$ac_MMDB_LIBS"

else 
  
  AC_MSG_FAILURE([ --with-mmdb-prefix specified but include files not found])

fi

AC_SUBST(MMDB_CXXFLAGS)
AC_SUBST(MMDB_LIBS)

])

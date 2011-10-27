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

 # should ideally be MMDB_CXXFLAGS="-I$MMDB_prefix/include", and the like
 # when MMDB and dependencies get installed - we infact, include both
 # directories.
 #  
ac_mmdb_dirs='
.
include
include/mmdb
lib
src
lib/src
lib/src/mmdb'
for ac_dir in $ac_mmdb_dirs; do
   if test -r "$mmdb_prefix/$ac_dir/mmdb_manager.h"; then
      ac_MMDB_CXXFLAGS="-I$mmdb_prefix/$ac_dir"
      break
   fi
done
 #
 # SGI compiler CC (CXX=CC) needs -lm to link maths library, but 
 # GCC c++ does not.
 #
 ac_MMDB_LDOPTS="-L$mmdb_prefix/$acl_libdirstem -L$mmdb_prefix -lmmdb -lm"
else
 # the compiler looks in the "standard" places for MMDB.  In real life,
 # it would be quite unlikely that MMDB would be installed in /usr/include, 
 # /usr/lib etc. so this code will not usually find the right dependencies.
 ac_MMDB_CXXFLAGS=""
 ac_MMDB_LDOPTS="-lmmdb -lm"
fi

AC_MSG_CHECKING([for MMDB])

	LIBS="$save_LIBS $ac_MMDB_LDOPTS"
	CXXFLAGS="$save_CXXFLAGS $ac_MMDB_CXXFLAGS"
	#
	# AC_TRY_LINK uses the c compiler (set by AC_LANG), so we will
	# temporarily reassign $CC to the c++ compiler.
 	#
	AC_LANG_PUSH(C++)
	AC_TRY_LINK([#include "mmdb_manager.h"] ,[ CMMDBManager a;  ], have_generic_mmdb=yes, have_generic_mmdb=no)
        AC_TRY_COMPILE([#include "mmdb_manager.h"] ,[ CAtom at; const char *name = "test"; at.SetAtomName(name); CMMDBManager *m; cpstr sg = m->GetSpaceGroup(); ], have_mmdb=yes, have_mmdb=no)
	# version 1.10 CISPEP code currently not used because of licence problem
	AC_TRY_COMPILE([#include "mmdb_manager.h"] ,[ CMMDBManager *m; m->SetFlag(MMDBF_IgnoreHash)  ], have_mmdb_ignore_hash=yes, have_mmdb_ignore_hash=no)	
	AC_TRY_COMPILE([#include "mmdb_manager.h"] ,[ CMMDBManager *m; m->GetModel(1)->GetNumberOfCisPeps(); ], have_mmdb_with_cispep=yes, have_mmdb_with_cispep=no)	
	AC_TRY_COMPILE([#include "mmdb_manager.h"] ,[ CLinkR c; ], have_mmdb_with_linkr=yes, have_mmdb_with_linkr=no)	
	AC_LANG_POP(C++)  # the language we have just quit
	AC_MSG_RESULT($have_mmdb)

if test x$have_mmdb = xyes; then

 LIBS="$saved_LIBS"
 if test $have_mmdb_ignore_hash = yes ; then
    HASH_FLAG=-DHAVE_MMDB_IGNORE_HASH
 else
    HASH_FLAG=
 fi
 if test "$have_mmdb_with_cispep" = "yes" ; then
    CISPEP_FLAG=-DHAVE_MMDB_WITH_CISPEP
 fi
 if test "$have_mmdb_with_linkr" = "no" ; then
    CLINKR_FLAG=-DMMDB_WITHOUT_LINKR
 fi
 CXXFLAGS="$saved_CXXFLAGS"
 MMDB_CXXFLAGS="$ac_MMDB_CXXFLAGS $HASH_FLAG $CISPEP_FLAG $CLINKR_FLAG"
 MMDB_LIBS="$ac_MMDB_LDOPTS"
ifelse([$1], , :, [$1])

else

 # no (CAtom checked) mmdb:

 if test x$have_generic_mmdb = xyes ; then
   echo Opps - You have an mmdb library, but it is out of date.
   echo You need version 1.12
 fi

 LIBS="$saved_LIBS"
 CXXFLAGS="$saved_CXXFLAGS"
 MMDB_LIBS=""
 MMDB_CXXFLAGS=""
 ifelse([$2], , :, [$2])

fi

AC_SUBST(MMDB_CXXFLAGS)
AC_SUBST(MMDB_LIBS)

])

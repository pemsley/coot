# mmdb-ssm.m4
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


AC_DEFUN([AM_WITH_CCP4SRS],
[AC_PROVIDE([AM_USE_CCP4SRS])


AC_ARG_WITH(ccp4srs-prefix, 
	AC_HELP_STRING( [--with-ccp4srs-prefix=PFX], [Prefix where CCP4SRS has been installed] ),
	[ with_ccp4srs_prefix="$withval" ],
 with_ccp4srs_prefix="")

AC_MSG_CHECKING([for CCP4SRS])

if test x$with_ccp4srs_prefix != x; then

   CCP4SRS_CXXFLAGS="-DHAVE_CCP4SRS"
   CCP4SRS_LIBS="-L$with_ccp4srs_prefix/$acl_libdirstem -lccp4srs"

ac_mmdb_dirs='
.
include
'

   for ac_dir in $ac_mmdb_dirs; do
      if test -r "$with_ccp4srs_prefix/$ac_dir/ccp4srs_types.h"; then
         ac_CCP4SRS_CXXFLAGS="-I$with_ccp4srs_prefix/$ac_dir"
         break
      fi
   done

  CCP4SRS_CXXFLAGS="$CCP4SRS_CXXFLAGS $ac_CCP4SRS_CXXFLAGS"
  
else 

   CCP4SRS_CXXFLAGS=""
   CCP4SRS_LIBS=""
   with_ccp4srs_prefix=no 

fi

AC_MSG_RESULT([$with_ccp4srs_prefix])

AC_SUBST(CCP4SRS_CXXFLAGS)
AC_SUBST(CCP4SRS_LIBS)
])

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


AC_DEFUN([AM_WITH_MMDBSSM],
[AC_PROVIDE([AM_USE_MMDBSSM])


AC_ARG_WITH(ssmlib-prefix, 
	AC_HELP_STRING( [--with-ssmlib-prefix=PFX], [Prefix where SSMLib has been installed] ),
	[ with_ssmlib_prefix="$withval" ],
          with_ssmlib_prefix=)

AC_MSG_CHECKING([for ssm library])

if test x$with_ssmlib_prefix != x; then

   if test -r "$with_ssmlib_prefix/include/ssm/ssm_superpose.h"; then

      MMDBSSM_CXXFLAGS="-DHAVE_SSMLIB"
      # AFAICS, ssm does not put its libs in $acl_libdirstem
      # MMDBSSM_LIBS="-L$with_ssmlib_prefix/$acl_libdirstem -lssm"
      MMDBSSM_LIBS="-L$with_ssmlib_prefix/lib -lssm"
      ac_MMDBSSM_CXXFLAGS="-I$with_ssmlib_prefix/include"

      # saved_LIBS="$LIBS"
      # saved_CPPFLAGS="$CPPFLAGS"
      # CPPFLAGS="$CPPFLAGS $ac_MMDBSSM_CPPFLAGS"
      # AC_CHECK_LIB([ssm], [ssm::DefineStructure])
      # LIBS="$saved_LIBS"
      # CPPFLAGS="$saved_CPPFLAGS"

   else 
      AC_MSG_FAILURE([ --with-ssmlib-prefix specified but include files not found])
   fi
   MMDBSSM_CXXFLAGS="$MMDBSSM_CXXFLAGS $ac_MMDBSSM_CXXFLAGS"
  
else 

   MMDBSSM_CXXFLAGS=""
   MMDBSSM_LIBS=""
   with_ssmlib_prefix=no 

fi

AC_MSG_RESULT([$with_ssmlib_prefix])

AC_SUBST(MMDBSSM_CXXFLAGS)
AC_SUBST(MMDBSSM_LIBS)
])

# mmdb-ssm.m4
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
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA


AC_DEFUN([AM_WITH_SSM],
[
AC_PROVIDE([AM_WITH_SSM])


AC_MSG_CHECKING([for ssm])

if ${PKG_CONFIG} ssm ; then 
   # dodgy hack for old bash
   LIBSSM_CXXFLAGS_TMP="$($PKG_CONFIG --cflags ssm)"
   LIBSSM_CXXFLAGS="$LIBSSM_CXXFLAGS_TMP -DHAVE_SSMLIB"
   LIBSSM_LIBS="$($PKG_CONFIG --libs ssm)"
   coot_found_ssm=yes
else
   coot_found_ssm=no
fi

AC_MSG_RESULT($coot_found_ssm)

if test $coot_found_ssm = no ; then
   AC_MSG_FAILURE([ssm not found])
fi

AC_SUBST(LIBSSM_CXXFLAGS)
AC_SUBST(LIBSSM_LIBS)

])


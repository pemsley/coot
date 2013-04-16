# coot-windows.m4
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


AC_DEFUN([AM_MINGW_WINDOWS],
[

AC_MSG_CHECKING([if this is MINGW on Windows])

 COOT_WINDOWS_CFLAGS=""
 have_windows_mingw=no

 # BL: workaround needed for new MinGW
 ac_cv_build_alias=${ac_cv_build_alias:=$build_alias}

 case $ac_cv_build_alias in 

  *-mingw*)
    COOT_WINDOWS_CFLAGS="-DWINDOWS_MINGW -DUSE_GNOME_CANVAS"
    have_windows_mingw=yes
    ;;
 esac

AC_MSG_RESULT([$have_windows_mingw])
AC_SUBST(COOT_WINDOWS_CFLAGS)

])


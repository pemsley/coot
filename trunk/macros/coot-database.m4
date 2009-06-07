# coot-database.m4
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

AC_DEFUN([AM_WITH_MYSQL_DATABASE], 
[AC_PROVIDE([AM_USE_DATABASE])

AC_ARG_WITH(database, [  --with-database  link with MYSQL database?],
   with_coot_database=true,
   with_coot_database="")

if test x$with_coot_database != x; then
   coot_database=true
   MYSQL_LIBS='-L/usr/lib/mysql -lmysqlclient'
   MYSQL_CFLAGS='-DUSE_MYSQL_DATABASE'
   # include files in /usr/include, so we don't need to specify that
   echo Linking with MYSQL database
else 
   echo Not using database
   coot_database=false
fi

# echo ............. in coot-database: MYSQL_LIBS=$MYSQL_LIBS
AC_SUBST(MYSQL_LIBS)
AC_SUBST(MYSQL_CFLAGS)
])


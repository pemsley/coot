# with-python.m4
# 
# Copyright 2004 The University of York
# Author: Paul Emsley
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
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


# AM_PATH_MMDB([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])

AC_DEFUN([AM_WITH_PYTHON], 
[AC_PROVIDE([AM_USE_PYTHON])

AC_MSG_CHECKING([for Python])


AC_ARG_WITH(python, [  --with-python=PFX Prefix where PYTHON has been installed],
 with_python_prefix="$withval",
 with_python_prefix="")


saved_LIBS="$LIBS"
saved_CFLAGS="$CFLAGS"

if test x$with_python != x; then

   #  ;^)
   #	

   # First check to see if python exists
   #
   python -c ''
   if test $? != 0 ; then
      echo python not found. 
      # This is an error because it was explicitly asked for
      # AC_MSG_ERROR([--with-python was specified, but no python was found])
	echo error
   fi
   
   # BL says: we need to have something special for MinGW since Python 
   # files are called different! E.g. no dots in version name
   # and put -lglob there too, dunno where else to put it
   
   # for MinGW
   if test $have_windows_mingw = yes ; then 
     py_cmd='import sys, os; drive, path = os.path.splitdrive(sys.prefix); drive = "/" + drive.replace(":",""); path = path.replace("\\", "/"); print drive + path + "/include"'
     PYTHON_CFLAGS="-DUSE_PYTHON -I`python -c "$py_cmd"`"
     py_cmd='import sys, os; drive, path = os.path.splitdrive(sys.prefix); drive = "/" + drive.replace(":",""); path = path.replace("\\", "/"); print drive + path + "/libs -lpython" + sys.version[[0]]+sys.version[[2]]'
     PYTHON_LIBS_PRE="-L`python -c "$py_cmd"`"
   else 
     # normal execution proceeds..
     PYTHON_CFLAGS="-DUSE_PYTHON -I`python -c 'import sys; print sys.prefix + "/include/python" + sys.version[[:3]]'`"
   # PYTHON_LIBS="-L/h/paule/build/lib/python2.2/config -lpython2.2 -lutil"
     PYTHON_LIBS_PRE="-L`python -c 'import sys; print sys.prefix + "/lib/python" + sys.version[[:3]] + "/config"'` -lpython`python -c 'import sys; print sys.version[[:3]]'`"
     break;;
   fi
	
   # now we have to deal with the -lutil issue.  On GNU/Linux, we need
   # it, MacOS X we don't (it doesn't exist and the compliation fails)
   # Cygwin is like MacOS X.
   # 
   # Recall that for RedHat, uname returns "Linux"
   #
   # consider using AM_CONDITIONAL, see the automake manual.
   # consider using config.guess output for python util check
   # config.guess returns e.g. i686-pc-linux-gnu, powerpc-apple-darwin7.0.0
   #
   UTIL_LIB="-lutil"
   py_uname=`uname`
   case	"$py_uname" in
	Darwin)
		UTIL_LIB=""
	;; 
	Cygwin*|CYGWIN*|cygwin*)
		UTIL_LIB=""
	;;
	# BL says:: same as for cygwin in mingw
	# dunno where else to put the glob at the moment
	MINGW*|Mingw*|mingw*)
		UTIL_LIB="-lglob"
	;;
   esac	

   PYTHON_LIBS="$PYTHON_LIBS_PRE $UTIL_LIB"
   AC_MSG_RESULT(yes)
   echo Cool, using Python
   coot_python=true	

   # BL says:: we shall check for PyGTK as well
   # let's try
   PKG_CHECK_MODULES(PYGTK, pygtk-2.0, have_pygtk=true, have_pygtk=false)
   AC_SUBST(PYGTK_CFLAGS)
   AC_SUBST(PYGTK_LIBS)

   if $have_pygtk; then
	AC_MSG_RESULT(yes)
   	echo Good we have pygtk
   	PYTHON_CFLAGS="$PYTHON_CFLAGS -DUSE_PYGTK"
   else
	AC_MSG_RESULT(no)
   	echo we dont have pygtk-2
   fi

else 

   AC_MSG_RESULT(no)
   echo Not using python
   PYTHON_CFLAGS=""
   PYTHON_LIBS=""
   # COOT_WRAP_PYTHON_CONVERT="cp"
   coot_python=false

fi

AC_SUBST(PYTHON_CFLAGS)
AC_SUBST(PYTHON_LIBS)

AM_CONDITIONAL(COOT_USE_PYTHON, test x$coot_python = xtrue)
])

dnl Stuart McNicholas writes:
dnl Yes, but make sure python is in your path.

dnl echo "import sys;print sys.prefix + '/include/python' + sys.version[:3]" | 
dnl python

dnl and 
dnl echo "import sys;print '-L' + sys.prefix + '/lib/python' + sys.version[:3] + 
dnl '/config' + ' -lpython' + sys.version[:3]" | python


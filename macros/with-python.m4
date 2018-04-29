# with-python.m4
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

   # First, add Daniel's indirection, python may not be python, it may
   # be $PYTHON:
   #
   if test x$PYTHON = x ; then
      PYTHON=python
   fi
   # Similar for python-config
   #
   if test x$PYTHON_CONFIG = x ; then
      PYTHON_CONFIG=python-config
   fi

   # Check to see if python exists:
   $PYTHON -c ''
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
      py_cmd='import sys, os; drive, path = os.path.splitdrive(sys.prefix); drive = "/" + drive.replace(":",""); path = path.replace("\\", "/"); print(drive + path + "/include")'
      PYTHON_CFLAGS="-DUSE_PYTHON -I`$PYTHON -c "$py_cmd"`"
      py_cmd='import sys, os; drive, path = os.path.splitdrive(sys.prefix); drive = "/" + drive.replace(":",""); path = path.replace("\\", "/"); print(drive + path + "/libs -lpython" + sys.version[[0]]+sys.version[[2]])'
      PYTHON_LIBS_PRE="-L`$PYTHON -c "$py_cmd"`"
   else 
      # normal execution proceeds..
      PYTHON_CFLAGS="-DUSE_PYTHON `$PYTHON_CONFIG --includes`"
      # PYTHON_LIBS="-L/h/paule/build/lib/python2.2/config -lpython2.2 -lutil"
      config_dir=`$PYTHON -c "import sys; print(sys.prefix + '/$acl_libdirstem/python' + sys.version[[:3]] + '/config')"`
      # echo  ======== config_dir: $config_dir
      PYTHON_LIBS_PRE="`$PYTHON_CONFIG --ldflags`"
   
      # extra hacking so that -ldl appears after -lpython2.x (needed
      # for correct linking on some systems)
      #
      for lib in $PYTHON_LIBS_PRE ;
      do
         case $lib in

           -lpython*)
              if test "$lib_hit" = "-ldl" ; then
                 lib_hit=python
              fi
              ;;

           -ldl)
              if test "$lib_hit" = python ; then
                 # do nothing, libs don't need fixing
                 :
              else 
                 lib_hit=-ldl
              fi
              ;;
         esac
     done

     if test "$lib_hit" = python ; then
        PYTHON_LIBS_PRE="$PYTHON_LIBS_PRE -ldl"
     fi

     # echo PYTHON_LIBS_PRE is $PYTHON_LIBS_PRE
     
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
   # BL: workaround needed for new MinGW
   ac_cv_build_alias=${ac_cv_build_alias:=$build_alias}

   case	$ac_cv_build_alias in
	*apple-darwin*)
		UTIL_LIB=""
	;; 
	Cygwin*|CYGWIN*|*cygwin*)
		UTIL_LIB=""
	;;
	# BL says:: same as for cygwin in mingw
	# dunno where else to put the glob at the moment
	MINGW*|Mingw*|*mingw*)
		#UTIL_LIB="-lglob"
		UTIL_LIB=""
	;;
   esac	

   # 20120516-PE
   # this is needed on Ubuntu - (or should that be gcc-4.6+?).  Not
   # needed on Ubuntu 10.04 (4.4.3) or (Fed 16 (4.6.3) without
   # --disable-shared)
   #
   # Now that lidia can dynamically add python modules, we need -ldl for Fed 14 (gcc 4.5.1)
   # 
   #
   coot_gcc_version=$(gcc -dumpversion)
   case "$coot_gcc_version" in 
      *.*.*)
       gcc_major_minor=${coot_gcc_version%.*}
        ;;  
      *.*)
        gcc_major_minor=$coot_gcc_version
        ;;
   esac


   save_CPPFLAGS="$CPPFLAGS"
   CPPFLAGS="$CPPFLAGS $PYTHON_CFLAGS"
   AC_CHECK_HEADER(Python.h, found_python_include_file=true, found_python_include_file=false)
   CPPFLAGS="$save_CPPFLAGS"
	
   # if found_python_include_file=false, then we can't compile with python or pygtk.

   if test "$found_python_include_file" = true ; then 

      PYTHON_LIBS="$PYTHON_LIBS_PRE $UTIL_LIB"

dnl    this results in a random yes on its own line.  Commenting it out.
dnl    echo Bmessage
dnl    AC_MSG_RESULT(yes)
dnl    echo done Bmessage

      echo Cool, using Python
      coot_python=true	

      # Check for PyGTK as well
      # let's try if --with-pygtk was given on the command line

      AC_ARG_WITH(pygtk, [  --with-pygtk=PFX Prefix where pygtk has been installed],
          with_pygtk_prefix="$withval",
          with_pygtk_prefix="no")

      if test "$with_pygtk_prefix" != "no" ; then

          PKG_CHECK_MODULES(PYGTK, pygtk-2.0, have_pygtk=true, have_pygtk=false)
          AC_SUBST(PYGTK_CFLAGS)
          AC_SUBST(PYGTK_LIBS)

          if $have_pygtk ; then
   	      echo Good we have pygtk
   	      PYTHON_CFLAGS="$PYTHON_CFLAGS -DUSE_PYGTK"
          else
	      AC_MSG_FAILURE([pygtk specified but not found])
              echo we dont have pygtk-2.0
          fi
      fi

   else 
      PYTHON_CFLAGS=""
      PYTHON_LIBS=""
      # COOT_WRAP_PYTHON_CONVERT="cp"
      coot_python=false
   fi

else 

dnl    AC_MSG_RESULT(no) yeuch.
   echo Not using python
   PYTHON_CFLAGS=""
   PYTHON_LIBS=""
   # COOT_WRAP_PYTHON_CONVERT is not used anywhere.  Perhaps it should be?
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


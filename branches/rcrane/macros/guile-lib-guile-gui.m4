
AC_DEFUN([AM_GUILE_LIB],
[
AC_MSG_CHECKING([for Guile-Lib])
if test -z "$ac_cv_path_GUILE" ; then 
   have_guile_lib=not_installed
else 
   $ac_cv_path_GUILE -c '(use-modules (sxml simple))'
   if test "$?" = 0  ; then 
      have_guile_lib=yes
   else 
      have_guile_lib=no
   fi
fi
AC_MSG_RESULT([$have_guile_lib])
if test "$have_guile_lib" = no  ; then 
   echo Must install guile-lib. http://home.gna.org/guile-lib/download/
   exit 2
fi
])


AC_DEFUN([AM_GUILE_GUI],
[
AC_MSG_CHECKING([for guile-gui])
if test -z "$ac_cv_path_GUILE" ; then 
   have_guile_gui=not_installed
else 
   $ac_cv_path_GUILE -c '(use-modules (gui paren-match))'
   if test "$?" = 0  ; then 
      have_guile_gui=yes
   else 
      have_guile_gui=no
   fi
fi
AC_MSG_RESULT([$have_guile_gui])
if test "$have_guile_gui" = no  ; then 
   echo Must install guile-gui for guile.
   exit 2
fi
])


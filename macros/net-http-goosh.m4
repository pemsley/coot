
AC_DEFUN([AM_GUILE_NET_HTTP],
[

dnl If Guile is configured for and found then ac_cv_path_GUILE is set to something sensible
dnl If that is the case, then we check for net-http, and gui (separately).
dnl If that is NOT the case (e.g. we configured only for python) then
dnl    we do not exit if net-http, or gui were not found.

AC_MSG_CHECKING([for net-hhtp])
if test -z "$ac_cv_path_GUILE" ; then 
   have_net_http=not_installed
else 
   $ac_cv_path_GUILE -c '(use-modules (oop goops) (oop goops describe) (net http))'
   if test "$?" = 0  ; then 
      have_net_http=yes
   else 
     have_net_http=no
   fi
fi
AC_MSG_RESULT([$have_net_http])
if test "$have_net_http" = no  ; then 
   echo Must install net-http for guile.
   exit 2
fi
])



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


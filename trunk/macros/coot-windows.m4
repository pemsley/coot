

AC_DEFUN([AM_MINGW_WINDOWS],
[

AC_MSG_CHECKING([if this is MINGW on Windows])

 COOT_WINDOWS_CFLAGS=""
 have_windows_mingw=no

 case $ac_cv_build_alias in 

  *-mingw*)
    COOT_WINDOWS_CFLAGS="-DWINDOWS_MINGW -DUSE_GNOME_CANVAS"
    have_windows_mingw=yes
    break;;
 esac

AC_MSG_RESULT([$have_windows_mingw])
AC_SUBST(COOT_WINDOWS_CFLAGS)

])


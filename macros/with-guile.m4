


AC_DEFUN([AM_WITH_GUILE], 
[AC_PROVIDE([AM_USE_GUILE])

# guile, not with-guile
AC_ARG_WITH(guile, [  --with-guile=PFX Prefix where GUILE has been installed],
 with_guile_prefix="$withval",
 with_guile_prefix="")


saved_LIBS="$LIBS"
saved_CFLAGS="$CFLAGS"

if test x$with_guile != x; then
   echo Congratulations, you are using Guile
   coot_guile=true
else 
   echo Not using guile
   coot_guile=false
fi

AM_CONDITIONAL(COOT_USE_GUILE, test x$coot_guile = xtrue)
])


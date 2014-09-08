
AC_DEFUN([AM_PATH_SQLITE3],
[
AC_ARG_WITH(sqlite3,[  --with-sqlite3 (optional)],
            sqlite3="$withval", sqlite3="")

    AC_MSG_CHECKING([for SQLite3])

    if ${PKG_CONFIG} sqlite3 ; then 
       SQLITE3_CFLAGS="-DUSE_SQLITE3 $($PKG_CONFIG --cflags sqlite3)"
       SQLITE3_LIBS="$($PKG_CONFIG --libs sqlite3)"
       coot_use_sqlite3=yes
    else
       SQLITE3_CFLAGS=""
       SQLITE3_LIBS=""
       coot_use_sqlite3=no
    fi

   AC_MSG_RESULT($coot_use_sqlite3)

   AC_SUBST(SQLITE3_CFLAGS)
   AC_SUBST(SQLITE3_LIBS)

])



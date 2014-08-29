
AC_DEFUN([AM_PATH_SQLITE3],
[
AC_ARG_WITH(sqlite3-prefix,[  --with-sqlite3-prefix=PFX   Prefix where SQLITE3 is installed (optional)],
            sqlite3_prefix="$withval", sqlite3_prefix="")

    AC_MSG_CHECKING([for SQLite3])

    if test -z "${PKG_CONFIG}"; then
      SQLITE3_CFLAGS=""
      SQLITE3_LIBS=""
      coot_use_sqlite3=no
    else
      SQLITE3_CFLAGS="-DUSE_SQLITE3 $($PKG_CONFIG --cflags sqlite3)"
      SQLITE3_LIBS=$($PKG_CONFIG --libs sqlite3)
      coot_use_sqlite3=yes
    fi

   AC_MSG_RESULT($coot_use_sqlite3)

   AC_SUBST(SQLITE3_CFLAGS)
   AC_SUBST(SQLITE3_LIBS)

])



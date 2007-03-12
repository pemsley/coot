
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


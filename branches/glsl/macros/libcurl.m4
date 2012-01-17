
dnl Copright 2002 by Daniel Stenberg 

dnl MY_CURL
dnl -------
dnl set my_cv_curl_vers to the version of libcurl or NONE
dnl if libcurl is not found or is too old

AC_DEFUN([COOT_CURL],[
 AC_CACHE_VAL(coot_cv_curl_vers,[
 coot_cv_curl_vers=NONE
 check="7.9.7"
 check_hex="070907"

 AC_MSG_CHECKING([for curl >= $check])

 if eval curl-config --version 2>/dev/null >/dev/null; then
   ver=`curl-config --version | sed -e "s/libcurl //g"`
   hex_ver=`curl-config --vernum | tr 'a-f' 'A-F'`
   ok=`echo "ibase=16; if($hex_ver>=$check_hex) $hex_ver else 0" | bc`

   if test x$ok != x0; then
     coot_cv_curl_vers="$ver"

     AC_MSG_RESULT([$coot_cv_curl_vers])
   else
     AC_MSG_RESULT(FAILED)
     AC_MSG_WARN([$ver is too old. Need version $check or higher.])
   fi
 else
   AC_MSG_RESULT(FAILED)
   AC_MSG_WARN([curl-config was not found])
 fi
	
])


if test $coot_cv_curl_vers != "NONE" ; then
   CURL_CFLAGS=$(curl-config --cflags)
   CURL_LIBS=$(curl-config --libs)
   USE_LIBCURL=-DUSE_LIBCURL
   AC_SUBST(CURL_LIBS)
   AC_SUBST(CURL_CFLAGS)
   AC_SUBST(USE_LIBCURL)
fi
])



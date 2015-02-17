

AC_INIT(configure.ac)

AM_INIT_AUTOMAKE(srs-inspect, 0.0)

AC_PROG_CXX

AC_MSG_CHECKING([which libtool initialization strategy to adopt])
AC_MSG_RESULT([m4_ifset([LT_INIT],[LT-INIT],[AC-PROG-LIBTOOL])])
m4_ifset([LT_INIT],
[LT_INIT],
[AC_PROG_LIBTOOL])

AM_WITH_CCP4SRS


AC_OUTPUT([
Makefile
])


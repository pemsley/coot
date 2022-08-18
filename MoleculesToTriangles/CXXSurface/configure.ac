# Process this file with autoconf to produce a configure script.
AC_INIT(libCXXSurface, 3.0, martin.noble at ncl.ac.uk)
# Prepare for automake
AM_INIT_AUTOMAKE
# Prepare for libtool
LT_INIT
# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_CONFIG_MACRO_DIRS([m4])
# Specify langauge C++
AC_LANG(C++)
AC_PROG_CXX
AC_LANG_PUSH([C++])

AC_COMPILE_IFELSE([AC_LANG_SOURCE(
  [[template <typename T>
    struct check
    {
      static_assert(sizeof(int) <= sizeof(T), "not big enough");
    };
    typedef check<check<bool>> right_angle_brackets;

    int a;
    decltype(a) b;

    typedef check<int> check_type;
    check_type c;
    check_type&& cr = static_cast<check_type&&>(c);]])],,
  AC_MSG_FAILURE(['$CXX $CXXFLAGS' does not accept ISO C++11]))

# Check for OpenMP
AC_OPENMP

AC_SUBST(MMDB_INCLUDE_DIR,/usr/local/include)
AC_SUBST(MMDB_LIB_DIR,/usr/local/lib)
AC_SUBST(CLIPPER_INCLUDE_DIR,/usr/local/include)
AC_SUBST(CLIPPER_LIB_DIR,/usr/local/lib)
AC_CONFIG_FILES(Makefile MoleculesToTrianglesCXXSurface.pc)
AC_OUTPUT

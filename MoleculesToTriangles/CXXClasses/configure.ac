# Process this file with autoconf to produce a configure script.
AC_INIT([MoleculesToTrianglesCXXClasses],[1.0],[martin.noble at ncl.ac.uk])
# Prepare for automake
AM_INIT_AUTOMAKE
# Check for availability of AR
AM_PROG_AR
# Prepare for libtool
LT_INIT
# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_CONFIG_MACRO_DIRS([m4])
AX_BOOST_BASE([1.48],,)
AX_BOOST_SYSTEM
AX_PYTHON_DEVEL
AX_BOOST_BASE
AX_BOOST_PYTHON
# Specify langauge C++
AC_LANG(C++)
AC_PROG_CXX
AC_LANG_PUSH([C++])
#Check for c++11 features
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
AX_CHECK_GL
AX_CHECK_GLU

PKG_CHECK_MODULES(MMDB, mmdb2 >= 2.0.1)
PKG_CHECK_MODULES(MoleculesToTrianglesCXXSurface, MoleculesToTrianglesCXXSurface >= 3)

# Check for OpenMP
AC_OPENMP

dnl :-)
dnl AC_SUBST(MMDB_INCLUDE_DIR,/usr/local/include)
dnl AC_SUBST(MMDB_LIB_DIR,/usr/local/lib)
dnl AC_SUBST(CLIPPER_INCLUDE_DIR,/usr/local/include)
dnl AC_SUBST(CLIPPER_LIB_DIR,/usr/local/lib)
dnl AC_SUBST(CXXSURFACE_INCLUDE_DIR,/usr/local/include)
dnl AC_SUBST(CXXSURFACE_LIB_DIR,/usr/local/lib)


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES(Makefile MoleculesToTrianglesCXXClasses.pc)
PKG_CHECK_MODULES([glfw3],[glfw3])
AC_OUTPUT

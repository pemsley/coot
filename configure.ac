dnl configure.in
dnl
dnl Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
dnl Copyright 2008 The University of Oxford
dnl Author: Paul Emsley
dnl
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 3 of the License, or (at
dnl your option) any later version.
dnl
dnl This program is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
dnl 02110-1301, USA
dnl
dnl This program is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
dnl implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
dnl
dnl Process this file with autoconf to produce a configure script.



AC_INIT([coot],[1.0.06-pre])
AC_CONFIG_SRCDIR(src/main.cc)
AM_INIT_AUTOMAKE

dnl Pick up the macros.
AM_ACLOCAL_INCLUDE(macros)
AC_CONFIG_MACRO_DIR(macros)

AC_SEARCH_LIBS([strerror],[cposix])
AC_PROG_CC
m4_warn([obsolete],
[The preprocessor macro `STDC_HEADERS' is obsolete.
  Except in unusual embedded environments, you can safely include all
  ISO C90 headers unconditionally.])dnl
# Autoupdate added the next two lines to ensure that your configure
# script's behavior did not change.  They are probably safe to remove.
AC_PROG_EGREP


AC_PROG_CXX
# LT_INIT that's old style. Now we use LT_INIT.
# It is not yet clear to me how AC_LANG_SOURCE should be used
# for LT_INIT
LT_INIT

# AC_PROG_INTLTOOL, no, not today.  Mac compilation problems, no perl/ XML or something.
# BL says Lets do some testing and use old style if needed
AC_MSG_CHECKING([which libtool initialization strategy to adopt])
AC_MSG_RESULT([m4_ifset([LT_INIT],[LT-INIT],[AC-PROG-LIBTOOL])])
m4_ifset([LT_INIT],
[LT_INIT],
[LT_INIT])

GNOME_COMPILE_WARNINGS
GNOME_CXX_WARNINGS

m4_define([v_git_revision_count], m4_esyscmd([(git rev-list --count HEAD 2>/dev/null || echo 0) | tr -d '\n']))
m4_define([v_git_commit], m4_esyscmd([git rev-parse HEAD | tr -d '\n']))
#dnl m4_define([vxx_git_revision_count], m4_esyscmd([git rev-list --count HEAD]))
# old gits can't use count with HEAD.
# we need to make the src and pyrogen directories if they do not exist
if test ! -d src     ; then mkdir src     ; fi
if test ! -d pyrogen ; then mkdir pyrogen ; fi
echo 'extern "C" { int git_revision_count() { return v_git_revision_count ; } }' > src/git-revision-count.cc
echo "#include <string>" >> src/git-revision-count.cc
echo 'std::string git_commit() { return "v_git_commit"; }' >> src/git-revision-count.cc
echo 'def revision_count():
      return ' v_git_revision_count > pyrogen/coot_git.py

AC_CHECK_HEADERS([sys/stdtypes.h], [COOT_STDTYPES="-DHAVE_SYS_STDTYPES_H=1"], [COOT_STDTYPES="-DHAVE_SYS_STDTYPES_H=0"])

dnl AC_OPENMP
AX_OPENMP
dnl thanks John Bollinger
AC_SUBST([OPENMP_CFLAGS])
AC_SUBST([OPENMP_CXXFLAGS])

# Note to self: must make these compulsory - not optional
#

AX_CXX_COMPILE_STDCXX_11
AX_CXX_HAVE_THREAD

# AX_CXX_COMPILE_STDCXX_11 sets HAVE_CXX11=1 (hopefully)
if test $HAVE_CXX11 = 1 ; then
   CXXFLAGS="$CXXFLAGS -std=c++11"
else
   echo we can\'t compile coot then
fi

AC_LIB_PREPARE_MULTILIB

AC_SUBST(COOT_STDTYPES)

eval EXPORT_DYNAMIC_FLAG_SPEC=$export_dynamic_flag_spec
AC_SUBST(EXPORT_DYNAMIC_FLAG_SPEC)


dnl We need to put a PKG_CHECK_MODULES in early so that pkg-config is
dnl setup first (it seems to be a side-effect of using PKG_CHECK_MODULES),
dnl setting ac_cv_env_PKG_CONFIG and related. This is because the conditional
dnl  use of PKG_CHECK_MODULES for GPP4 causes problems when setting up pkg-config
dnl
dnl So, let's check GLIB here, should be harmless.
dnl
PKG_CHECK_MODULES([GLIB], [glib-2.0])


dnl test for WINDOWS_MINGW
AM_MINGW_WINDOWS
dnl test for LIBPNG (Windows only - maybe shouif here rather than in macro?)
if test "$have_windows_mingw" = yes; then
  AM_PATH_GLOB
  AC_CONFIG_FILES(pyrogen/pyrogen.bat)
fi

# This seems not to be returning the correct value
# AM_PATH_LIBPNG
# so try this:
PKG_CHECK_MODULES([LIBPNG], [libpng >= 1.2])
if test -z "$LIBPNG_LIBS" ; then
   :
else
   LIBPNG_CFLAGS="-DUSE_LIBPNG=1 $LIBPNG_CFLAGS"
fi

# and while we are thinking about images, let's check for Cairo also
PKG_CHECK_MODULES([CAIRO], [cairo >= 1.14])


AM_PATH_SQLITE3

dnl =========================== crystallographic libraries ================================
dnl Check for gpp4

# I had the syntax wrong - it was allways passing
AC_ARG_WITH([backward], AS_HELP_STRING([--with-backward Build with the backward header (better debugging info on GNU/Linux)]), with_backward=true, with_backward=false)
if test "x$with_backward" = xtrue ; then
   USE_BACKWARD=use_backward
fi
AM_CONDITIONAL(USE_BACKWARD, [test x$USE_BACKWARD = xuse_backward])    # tested in utils/Makefile.am

AC_ARG_WITH([libdw], AS_HELP_STRING([--with-libdw build backward with libdw support]), with_libdw=true, with_libdw=false)
if test x$with_libdw == xtrue ; then
   BACKWARD_HAS_DW=true
fi
AM_CONDITIONAL(BACKWARD_HAS_DW, [test x$with_libdw = xtrue])    # tested in utils/Makefile.am

dnl Set platform-specific compilation flags - derived from clipper's configure
test "${target:+set}" = set || target="$host_os"
case "$target" in
   *linux* ) ;;
   *irix* )
	if test "`basename $CXX`" = CC; then
	   case "$CXXFLAGS" in
	      *LANG::std* ) ;;
              * )
		  CXXFLAGS="$CXXFLAGS -LANG:std -ptused"
           esac
           # This patches libtool... This change exists in libtool 1.4d and
           # is required such that -LANG: is not filtered out...
           AC_MSG_NOTICE([Patching libtool for irix CC])
           sed "s:no/\*-\*-irix\*):no/*-*-irix* | /*-*-irix*):" libtool > libtool.tp
           mv libtool.tp libtool
           chmod +x libtool
	fi
    ;;
esac


dnl Check for MMDB.
AM_PATH_MMDB2(, [echo Error locating necessary mmdb; exit 1],,)

AM_WITH_SSM
AM_WITH_CCP4SRS

AM_SINGLE_FFTW2

dnl Check for Clipper.  we exit if clipper is not found.
#
AM_PATH_CLIPPER(, [echo Error locating Clipper - a CNS-enabled version of Clipper \(2.1-090520-ac or later\) is required now; exit 1],AC_MSG_ERROR([Cannot find Clipper]))

# 20210313-PE-merge-weirdness
# we set default coot_guile to true
coot_guile=true
AM_WITH_GUILE(, [echo optional guile\?],,)
if test x$coot_guile = xtrue ; then
   GUILE_FLAGS
   AM_GUILE_LIB
   # guile-gui depends on a gui interface being present.  If it is, force check for guile-gui
   if test x$with_gnome != x ; then
      # AM_GUILE_GNOME
      :
   fi
   GUILE_VERSION=$($ac_cv_path_GUILE_CONFIG --version  2>&1 | awk '{split($NF, arr, "[[.]]");  print arr[[1]] "." arr[[2]]}')
   echo checking GUILE_VERSION: $GUILE_VERSION
   AC_SUBST(GUILE_VERSION)
fi


dnl =========================== GSL ================================

dnl test for GSL
# min_gsl_version=1.3
AX_PATH_GSL([1.3],[COOT_GSL_FLAGS=-DHAVE_GSL],[COOT_GSL_FLAGS=""])
AC_SUBST(COOT_GSL_FLAGS)

#GLM
AC_LANG_PUSH(C++)
AC_CHECK_HEADER([glm/glm.hpp],, AC_MSG_ERROR(GLM headers are required.))
AC_LANG_POP


dnl =========================== GTK, GOOCANVAS ===========================

AC_ARG_WITH(gtk4, [  --gtk4 build coot with gtk4 ], with_gtk4="yes", with_gtk4="no")

if test $with_gtk4 = yes ; then
   PKG_CHECK_MODULES(GTK, gtk4 >= 4.1)
else
   PKG_CHECK_MODULES(GTK, gtk+-3.0 >= 3.18)
fi

# goocanvas for GTK4?
# PKG_CHECK_MODULES([GOOCANVAS], [goocanvas-2.0]) # not for the moment - while trying to compile vs GTK4


PKG_CHECK_MODULES(GLFW3, glfw3 >= 3.0)

# this needs a rethink
MAKE_COOTRC_CMD="sed -e 's/\#GTK2 //' "
AC_SUBST(MAKE_COOTRC_CMD)



AM_PATH_PYTHON([3.7])
AX_PYTHON_DEVEL

AX_BOOST_BASE(,,[AC_MSG_ERROR([No boost found])])

# Boost-Python is needed for RDKit interface.
# Actually, this should only be tested if enhanced-ligand-tools was passed - not all the time as
# we currently have. And in that case, it should be
# AX_BOOST_PYTHON(,,[AC_MSG_ERROR([No boost-python found])]) because it is required.
# FIXME Later.
#
AX_BOOST_PYTHON


if test -z "$BOOST_PYTHON_LIB" ; then
   echo BOOST_PYTHON_LIB has not been set. Did you pass --with-boost-libdir\?
   exit
fi

# So now BOOST_PYTHON_LIB is set (and that is typically "boost_python").
# How do we use that?  At the moment we want to only use
# -l$(BOOST_PYTHON_LIB) when we have RDKIT_LIBS set

if test -z "$RDKIT_LIBS" ; then
   :
else
   if test -z "$BOOST_PYTHON_LIB" ; then
      :
   else
      L_BOOST_PYTHON=-l$BOOST_PYTHON_LIB
      if test "$have_windows_mingw" = yes; then
         COOT_WINDOWS_CFLAGS+=" -DBOOST_SYSTEM_NO_DEPRECATED"
      fi
   fi
fi
AC_SUBST(L_BOOST_PYTHON)


# Goocanvas is no longer optional
# AM_PATH_GOOCANVAS
# if test "$GOOCANVAS_CFLAGS" != "" ; then
#   GOOCANVAS_CFLAGS="-DHAVE_GOOCANVAS $GOOCANVAS_CFLAGS"
# fi

# this can go away in 2020 - no one will try to use that old RDKit.
#
# was RDKIT used and if it was does it have Cairo support (if it does
# have cairo support, we can use it to generate molecular (chemical)
# depictions.
#
if test -n "$RDKIT_LIBS" ; then

   case "$RDKIT_LIBS" in

     *MolDraw2D*)
          RDKIT_CXXFLAGS="$RDKIT_CXXFLAGS -DRDKIT_HAS_CAIRO_SUPPORT"
          ;;
   esac
fi


AC_ARG_WITH(enhanced-ligand-tools, [  --with-enhanced-ligand-tools with optional ligand-related dependencies ],
 coot_enhanced_ligand_additions="true",
 coot_enhanced_ligand_additions="false")

# but we can't use enhanced-ligand-tools without python:
if $coot_enhanced_ligand_additions = true ; then
   if test -n "$COOT_USE_PYTHON_TRUE" ; then
      echo Configuration error: cant use enhanced-ligand-tools without python
      exit 1
   fi
fi

AC_ARG_WITH(rdkit-prefix, [location of the RDKit package], rdkit_prefix="$withval", rdkit_prefix="")

if test x$rdkit_prefix = x ; then
  :
else
   rdkit_install_dir="$withval"
   # better RDKit linking I hope
   RDKIT_CXXFLAGS="-I$rdkit_install_dir/include/rdkit -DRDKIT_HAS_CAIRO_SUPPORT"
   RDKIT_LIBS="-L$rdkit_install_dir/lib -lRDKitMolDraw2D -lRDKitForceFieldHelpers -lRDKitDescriptors -lRDKitForceField -lRDKitSubstructMatch -lRDKitOptimizer -lRDKitDistGeomHelpers -lRDKitDistGeometry -lRDKitChemReactions -lRDKitAlignment -lRDKitEigenSolvers -lRDKitDepictor -lRDKitcoordgen -lRDKitMolChemicalFeatures -lRDKitSubgraphs -lRDKitPartialCharges -lRDKitFileParsers -lRDKitRDGeometryLib -lRDKitGraphMol -lRDKitShapeHelpers -lRDKitFingerprints -lRDKitMolAlign -lRDKitMolTransforms -lRDKitChemTransforms -lRDKitSmilesParse -lRDKitFilterCatalog -lRDKitCatalogs -lRDKitSubgraphs -lRDKitPartialCharges -lRDKitRDGeneral -lRDKitDataStructs -lboost_serialization -l$BOOST_PYTHON_LIB $PYTHON_LIBS"
fi

if test x$coot_enhanced_ligand_additions = xtrue ; then
   ENHANCED_LIGAND_TOOLS_CPPFLAGS="-DMAKE_ENHANCED_LIGAND_TOOLS"
else
   ENHANCED_LIGAND_TOOLS_CPPFLAGS=""
fi

echo ENHANCED_LIGAND_TOOLS_CPPFLAGS $ENHANCED_LIGAND_TOOLS_CPPFLAGS

AC_ARG_WITH(lidia-only, [  --with-lidia-only build lidia only (not coot) ],
 lidia_only="yes",
 lidia_only="no")
AM_CONDITIONAL([COND_LIDIA_ONLY], [test "$lidia_only" = yes])

AM_CONDITIONAL(MAKE_ENHANCED_LIGAND_TOOLS, test x$coot_enhanced_ligand_additions = xtrue)
AC_SUBST(ENHANCED_LIGAND_TOOLS_CPPFLAGS)

AC_SUBST(RDKIT_CXXFLAGS)
AC_SUBST(RDKIT_LIBS)

AC_SUBST(MOLTRIANGLES_CPPFLAGS)
AC_SUBST(MOLTRIANGLES_LIBS)

AC_SUBST(COOT_BUILD_INFO_STRING)


# #defines for SWIG:
#
# SWIG_FLAGS is something of a misnomer.  It actually defines
# COOT_USE_GTK2_INTERFACE, which is needed to correctly declare the
# correct callbacks in both GTK1 and GTK2 regimes.
#
SWIG_FLAGS="-DCOOT_USE_GTK2_INTERFACE"
AC_SUBST(SWIG_FLAGS)


dnl check for guile modules www and gtk here.

dnl check for GSL
dnl AC_CHECK_LIB(m,main)
dnl AC_CHECK_LIB(gslcblas,main)
dnl AC_CHECK_LIB(gsl,main)

#GEMMI
AC_LANG_PUSH(C++)
AC_CHECK_HEADER([gemmi/to_json.hpp], AC_MSG_NOTICE(GEMMI headers were found), AC_MSG_NOTICE(GEMMI headers were not found),)
if test "$ac_cv_header_gemmi_to_json_hpp" = yes ; then
   GEMMI_CXXFLAGS=-DUSE_GEMMI
fi
AC_SUBST(GEMMI_CXXFLAGS)
AC_LANG_POP


AM_WITH_MYSQL_DATABASE

COOT_CURL

# 20090613 We need to be able to determine the build type, and embed
# that value in to coot binary.  We need to do that for the updating
# Coot.  THe build type will be set by COOT_SYS_BUILD_TYPE.  Which is
# for now one variable.  I may split it into os/os-ver/python/arch etc
# later.
AM_COOT_SYS_BUILD_TYPE


dnl Fix this another time.
dnl dnl Add the languages which your application supports here.
dnl ALL_LINGUAS="es es_ES"
dnl # AM_GNU_GETTEXT
dnl GETTEXT_PACKAGE=coot
dnl dnl we need to export ALL_LINGUAS here otherwise it is not expanded
dnl dnl for GTk1 on cycle at least, so that we would get:
dnl dnl No rule to make target `@ALL_LINGUAS@.po', needed by `@ALL_LINGUAS@ .gmo'.
dnl AC_SUBST(ALL_LINGUAS)
dnl AC_SUBST(GETTEXT_PACKAGE)
dnl AM_GLIB_GNU_GETTEXT

dnl dnl For now: just include the declarations for config.h.in
dnl AH_TEMPLATE([ENABLE_NLS], [Define if NLS support is used])
dnl AH_TEMPLATE([PACKAGE_LOCALE_DIR], [Directory of locales])


dnl Set PACKAGE_LOCALE_DIR in config.h.
# if test "x${prefix}" = "xNONE"; then
  # AC_DEFINE_UNQUOTED(PACKAGE_LOCALE_DIR, "${ac_default_prefix}/${DATADIRNAME}/locale")
# else
  # AC_DEFINE_UNQUOTED(PACKAGE_LOCALE_DIR, "${prefix}/${DATADIRNAME}/locale")
# fi

AC_CONFIG_FILES([
Makefile
macros/Makefile
pdb-templates/Makefile
scheme/Makefile
python/Makefile
compat/Makefile
utils/Makefile
coot-utils/Makefile
protein_db/Makefile
ccp4mg-utils/Makefile
coords/Makefile
db-main/Makefile
density-contour/Makefile
angles/Makefile
mini-mol/Makefile
skeleton/Makefile
ligand/Makefile
geometry/Makefile
ideal/Makefile
build/Makefile
analysis/Makefile
high-res/Makefile
cootaneer/Makefile
cootilus/Makefile
greg-tests/Makefile
python-tests/Makefile
lidia-core/Makefile
lbg/Makefile
doc/Makefile
src/Makefile
src/coot
python/coot_load_modules.py
data/Makefile
data/rama-coeffs/Makefile
data/metal/Makefile
data/cho-models/Makefile
data/cho-acedrg/Makefile
glTF/Makefile
pixmaps/Makefile
pixmaps/icons/Makefile
textures/Makefile
rama-data/Makefile
man/Makefile
pli/Makefile
fonts/Makefile
pyrogen/Makefile
pyrogen/pyrogen
MoleculesToTriangles/Makefile
MoleculesToTriangles/CXXSurface/Makefile
MoleculesToTriangles/CXXSurface/MoleculesToTrianglesCXXSurface.pc
MoleculesToTriangles/CXXClasses/Makefile
MoleculesToTriangles/CXXClasses/MoleculesToTrianglesCXXClasses.pc
auxil/Makefile
auxil/catcheck/Makefile
auxil/rama/Makefile
auxil/sec-str/Makefile
auxil/strand-analysis/Makefile
auxil/vectors/Makefile
rcrane/Makefile
])
AC_OUTPUT
# removed setup/Makefile
# po/Makefile.in
# lbg/Makefile

# GETTEXT related things: When AM_GNU_GETTEXT get enabled, these get enabled.
# po/Makefile
# intl/Makefile - sigh.  Fix this later, not now.


# or maybe:
#],[sed -e "/POTFILES =/r po/POTFILES" po/Makefile.in > po/Makefile])

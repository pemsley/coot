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


AC_INIT([coot],[1.1.x-pre])
AC_CONFIG_SRCDIR(src/main.cc)
AM_INIT_AUTOMAKE

dnl Pick up the macros.
AM_ACLOCAL_INCLUDE(macros)
AC_CONFIG_MACRO_DIR(macros)

AC_PROG_CXX
# AC_PROG_LIBTOOL? that's old style. Now we use LT_INIT.
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

# m4_define([v_git_revision_count], m4_esyscmd([(git rev-list --count HEAD 2>/dev/null || echo 0) | tr -d '\n']))
# m4_define([v_git_commit], m4_esyscmd([git rev-parse HEAD | tr -d '\n']))
# autoconf and aclocal can be run outside of a git repo, so give the git command some protection and fallback
m4_define([v_git_revision_count], m4_esyscmd([bash ./get-the-git-commit-count.sh | tr -d '\n']))
m4_define([v_git_commit], m4_esyscmd([bash ./get-the-git-commit.sh | tr -d '\n']))

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
   # C++20 is required by newest versions of RDKit (since March 2025)
   # Coot itself should compile with C++17
   CXXFLAGS="$CXXFLAGS -std=c++20"
else
   echo We can\'t compile Coot without C++17 \(C++20 ideally\)
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

uname_s=$(uname -s)
AM_CONDITIONAL([OS_DARWIN], [test x$uname_s = xDarwin ])

# this variable is used by the bin/coot script:
COOT_BINARY_NAME=coot-1
if test x$uname_s = xDarwin ; then
   COOT_BINARY_NAME=MacCoot
fi
AC_SUBST(COOT_BINARY_NAME)


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

   case "${host_os}" in
    *darwin* )
      PKG_CHECK_MODULES(LIBDW, libdwarf)
      ;;
    *linux* )
      PKG_CHECK_MODULES(LIBDW, libdw >= 0.20)
      ;;
   esac
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

# Let's not think about CCP4SRS. I could never get the beast to work.
# AM_WITH_CCP4SRS

AM_SINGLE_FFTW2

echo PKG_CONFIG_PATH is $PKG_CONFIG_PATH
echo debug clipper:
echo pkg-config --modversion clipper
pkg-config --modversion clipper
echo pkg-config --libs clipper
pkg-config --libs clipper

dnl Check for Clipper.  we exit if clipper is not found.
#
AM_PATH_CLIPPER(, [echo Error locating Clipper - a CNS-enabled version of Clipper \(2.1-090520-ac or later\) is required now; exit 1],AC_MSG_ERROR([Cannot find Clipper]))

dnl # 20210313-PE-merge-weirdness
dnl # we set default coot_guile to true
dnl coot_guile=true
dnl AM_WITH_GUILE(, [echo optional guile\?],,)
dnl if test x$coot_guile = xtrue ; then
dnl    GUILE_FLAGS
dnl    AM_GUILE_LIB
dnl    # guile-gui depends on a gui interface being present.  If it is, force check for guile-gui
dnl    if test x$with_gnome != x ; then
dnl       # AM_GUILE_GNOME
dnl       :
dnl    fi
dnl    GUILE_VERSION=$($ac_cv_path_GUILE_CONFIG --version  2>&1 | awk '{split($NF, arr, "[[.]]");  print arr[[1]] "." arr[[2]]}')
dnl    echo checking GUILE_VERSION: $GUILE_VERSION
dnl    AC_SUBST(GUILE_VERSION)
dnl fi

AC_ARG_WITH([guile], AS_HELP_STRING([--with-guile]), with_guile=true, with_guile=false)
if test "x$with_guile" = xtrue ; then
   GUILE_FLAGS
fi
AM_CONDITIONAL(USE_GUILE, [test x$USE_GUILE = xuse_guile])



if test "x$GUILE_CFLAGS" = x ; then
   :
else
   COOT_USE_GUILE="-DUSE_GUILE"
   AC_SUBST(COOT_USE_GUILE)
fi

dnl =========================== GSL ================================

PKG_CHECK_MODULES([gsl], [gsl])

save_libs="$LIBS"
LIBS="$gsl_LIBS"
save_cflags="$CFLAGS"
CFLAGS="$gsl_CFLAGS"
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
LIBS="$save_libs"
CFLAGS="$save_cflags"

# I don't understand how this test works
# AX_PATH_GSL(GSL_VERSION, [], [AC_MSG_ERROR(could not find required version of GSL)])

# AX_PATH_GSL(1.3, [], [AC_MSG_ERROR(could not find required version of GSL)])


dnl =========================== GLM ================================

# While GLM has no pkgconfig glm.pc file, we have to be a bit more hand-to-mouth
# about its configuration
# PKG_CHECK_MODULES(GLM, glm >= 0.9.9)

# Copied and adapted from how gemmi and rdkit are found below
AC_ARG_WITH([glm], AS_HELP_STRING([--with-glm location of the GLM package]), with_glm="$withval", with_glm=false)

AC_MSG_CHECKING([for GLM])
if test "x$with_glm" = xfalse ; then
   if test -x $prefix/include/glm ; then
      #echo OK, good we found glm in $prefix
      AC_MSG_RESULT(yes)
      glm_install_dir=$prefix
      AC_SUBST(GLM_CXXFLAGS, "-I$glm_install_dir/include")
      AC_SUBST(GLM_LIBS, "-L$glm_install_dir/lib -lglm")
   else
      AC_MSG_RESULT(no)
      echo Configuration error. Could not find GLM.
      exit 1
   fi
else
   AC_MSG_RESULT(yes)
   AC_SUBST(GLM_CXXFLAGS, "-I$withval/include")
   AC_SUBST(GLM_LIBS, "-L$withval/lib -lglm")
fi

dnl =========================== GTK ===========================

PKG_CHECK_MODULES(GTK, gtk4 >= 4.4)

PKG_CHECK_MODULES(EPOXY, epoxy >= 1.5)

dnl =========================== Python =========================

AM_PATH_PYTHON([3.7])
AX_PYTHON_DEVEL
# this is not the canonical way to use arguments 2 and 3, but it works in this case
# AX_PYTHON_MODULE([requests], [echo requests found], [echo requests not found])

dnl =========================== PyGObject =========================

PKG_CHECK_MODULES(PyGObject, pygobject-3.0 >= 3.36)

dnl =========================== Boost =========================

AX_BOOST_BASE(,,[AC_MSG_ERROR([No boost found])])

AX_BOOST_THREAD

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
# How do we use that?
if test -z "$BOOST_PYTHON_LIB" ; then
    :
else
    L_BOOST_PYTHON=-l$BOOST_PYTHON_LIB
    if test "$have_windows_mingw" = yes; then
         COOT_WINDOWS_CFLAGS+=" -DBOOST_SYSTEM_NO_DEPRECATED"
    fi
fi
AC_SUBST(L_BOOST_PYTHON)

dnl =========================== RDKit =========================

# this does nothing now - keep it for backwards compatibility
AC_ARG_WITH(enhanced-ligand-tools, [  --with-enhanced-ligand-tools with optional ligand-related dependencies ],
 coot_enhanced_ligand_additions="true",
 coot_enhanced_ligand_additions="false")

# but we can't use enhanced-ligand-tools without python:
if $coot_enhanced_ligand_additions = true ; then
   echo no longer need to specify --with-enhanced-ligand-tools
   if test -n "$COOT_USE_PYTHON_TRUE" ; then
      echo Configuration error: cant use enhanced-ligand-tools without python
      exit 1
   fi
fi


AC_ARG_WITH(rdkit-prefix, AS_HELP_STRING([--with-rdkit-prefix location of the RDKit package]),
            rdkit_prefix="$withval", rdkit_prefix=false)

AC_SUBST(EMBEDDED_PYTHON_LIBS, "$PYTHON_LIBS")
pl="$PYTHON_LIBS"

case "${host_os}" in
    *darwin* )
      pl="-undefined dynamic_lookup"
      # python libs for coot-bin executable (embedding Python) will need regular linking
      # Fix that latter.
      PYTHON_LIBS="-undefined dynamic_lookup"
      ;;
esac

if test x$rdkit_prefix = xfalse ; then
   echo checking for RDKit in $prefix
   if test -x $prefix/include/rdkit ; then
      echo OK, good we found RDKit in $prefix
      rdkit_install_dir=$prefix
      RDKIT_CXXFLAGS="-I$rdkit_install_dir/include/rdkit -DRDKIT_HAS_CAIRO_SUPPORT"
   else
      echo Configuration error. No, bad... we do not have necessary RDKIT
      exit 1
   fi
else
   rdkit_install_dir="$withval"
   # better RDKit linking I hope
   RDKIT_CXXFLAGS="-I$rdkit_install_dir/include/rdkit -DRDKIT_HAS_CAIRO_SUPPORT"
fi


# test here for presence of coordgen in the RDKit.
# To clarify, this is not related to how to configure Coot, it is related to how the RDKit was configured (and linked).
# We need to clarify that so that Coot can be linked correctly.
librdkitcoordgen=
AC_ARG_WITH([coordgen], AS_HELP_STRING([--with-coordgen Specify that the RDKit was compiled with coordgen]), with_rdkit_coordgen="$withval", with_rdkit_coordgen=false)
if test -n "$with_rdkit_coordgen"; then
    librdkitcoordgen=-lRDKitcoordgen
fi

# test here for presence of Inchi support in the RDKit.
# To clarify, this is not related to how to configure Coot, it is related to how the RDKit was configured (and linked).
# We need to clarify that so that Coot can be linked correctly.
librdkitinchi=
if test -f "$rdkit_install_dir/include/rdkit/GraphMol/inchi.h"; then
    echo Enabling RDKit Inchi support.
    # echo If this step causes linking problems, consider adding -linchi or -lRDKitInchi
    librdkitinchi="-lRDKitRDInchiLib"
else
   echo RDKit Inchi support disabled.
fi

RDKIT_LIBS="-L$rdkit_install_dir/lib -lRDKitMolDraw2D -lRDKitForceFieldHelpers -lRDKitDescriptors -lRDKitForceField -lRDKitSubstructMatch -lRDKitOptimizer -lRDKitDistGeomHelpers -lRDKitDistGeometry -lRDKitChemReactions -lRDKitAlignment -lRDKitEigenSolvers -lRDKitDepictor $librdkitcoordgen -lRDKitMolChemicalFeatures -lRDKitPartialCharges -lRDKitFileParsers -lRDKitRDGeometryLib -lRDKitGraphMol $librdkitinchi -lRDKitShapeHelpers -lRDKitFingerprints -lRDKitMolAlign -lRDKitMolTransforms -lRDKitChemTransforms -lRDKitSmilesParse -lRDKitGenericGroups -lRDKitRingDecomposerLib -lRDKitFilterCatalog -lRDKitCatalogs -lRDKitSubgraphs -lRDKitPartialCharges -lRDKitRDGeneral -lRDKitDataStructs -lboost_serialization -lboost_iostreams"

# this must always be true now
ENHANCED_LIGAND_TOOLS_CPPFLAGS="-DMAKE_ENHANCED_LIGAND_TOOLS"

# echo ENHANCED_LIGAND_TOOLS_CPPFLAGS $ENHANCED_LIGAND_TOOLS_CPPFLAGS

AC_ARG_ENABLE(layla-only, [  --enable-layla-only build layla only (not coot) ],
 layla_only="yes",
 layla_only="no")
AM_CONDITIONAL([COND_LAYLA_ONLY], [test "$layla_only" = yes])

AC_ARG_ENABLE(libcootapi-only, [  --enable-libcootapi-only build libcootapi only (not coot) ],
[case "${enableval}" in
  yes) libcootapi_only=true ;;
  no)  libcootapi_only=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
esac],[libcootapi_only=false])
AM_CONDITIONAL([COND_LIBCOOTAPI_ONLY], [test "$libcootapi_only" = true])

AC_ARG_ENABLE(relocatable-binaries, [  --enable-relocatable-binaries try to make the binaries work independent of the installation directory],
   relocatable=yes, relocatable=no)
AM_CONDITIONAL([RELOCATABLE_BINARIES], [test "$relocatable" =  yes])

# echo testingA libcootapi_only is $libcootapi_only
# echo testingA COND_LIBCOOTAPI_ONLY is $COOT_LIBCOOTAPI_ONLY
# echo testingB lidia_only is $lidia_only

AM_CONDITIONAL(MAKE_ENHANCED_LIGAND_TOOLS, test x$coot_enhanced_ligand_additions = xtrue)
AC_SUBST(ENHANCED_LIGAND_TOOLS_CPPFLAGS)

AC_SUBST(RDKIT_CXXFLAGS)
AC_SUBST(RDKIT_LIBS)

AC_SUBST(MOLTRIANGLES_CPPFLAGS)
AC_SUBST(MOLTRIANGLES_LIBS)

AC_SUBST(COOT_BUILD_INFO_STRING)

dnl check for guile modules www and gtk here.

# gemmi doesn't have a pkgconfig file. Frustrating.
AC_ARG_WITH([gemmi], AS_HELP_STRING([--with-gemmi location of the GEMMI package]), with_gemmi="$withval", with_gemmi=false)

# is this the wrong logic actually, should the checking message only be seen if with_gemmi is not false?
AC_MSG_CHECKING([for GEMMI])
if test "x$with_gemmi" = xfalse ; then
    AC_MSG_RESULT(no)
else
    AC_LANG_PUSH(C++)
    GEMMI_CXXFLAGS="-I$with_gemmi/include -DUSE_GEMMI=1"
    save_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $GEMMI_CXXFLAGS"
    AC_CHECK_HEADER([gemmi/to_json.hpp], AC_MSG_NOTICE(GEMMI headers were found), AC_MSG_FAILURE(GEMMI headers were not found),)
    CPPFLAGS="$save_CPPFLAGS"
    AC_LANG_POP
    AC_MSG_RESULT(yes)
    AC_SUBST(GEMMI_CXXFLAGS)
    AC_SUBST(GEMMI_LIBS, "-L$withval/lib -lgemmi_cpp")
fi

# 20090613 We need to be able to determine the build type, and embed
# that value in to coot binary.  We need to do that for the updating
# Coot.  THe build type will be set by COOT_SYS_BUILD_TYPE.  Which is
# for now one variable.  I may split it into os/os-ver/python/arch etc
# later.
AM_COOT_SYS_BUILD_TYPE


AC_ARG_WITH([sound], AS_HELP_STRING([--with-sound compile with ogg-vorbis and OpenAL libs]), with_sound="withval", with_sound=false)
AC_MSG_CHECKING([for sound])
if test "x$with_sound" = xfalse ; then
   :
   AC_MSG_RESULT(no)
else
   AC_MSG_RESULT(yes)
   PKG_CHECK_MODULES(VORBISFILE, [vorbisfile >= 1.2])
   PKG_CHECK_MODULES(OpenAL, [OpenAL >= 1.23])
   # uname_s=$(uname -s) # should autoconf do this for us? If so, I can't see where it is stored
   # if test "$uname_s" = Linux ; then
   #   PKG_CHECK_MODULES(ALSA, [alsa >= 1.2])
   #else
   #   # it's a mac then
   #   PKG_CHECK_MODULES(OpenAL, [OpenAL >= 1.23])
   #fi
   # Use Core Audio on macOS. Consider PortAudio instead of using alsa directly.
   AC_SUBST(SOUND_FLAGS, -DWITH_SOUND)
fi


PKG_CHECK_MODULES([RSVG], [librsvg-2.0 >= 2.40], [
                            AC_DEFINE(HAVE_RSVG, 1, [Define if librsvg is available])
], [
])

# Add the compiler and linker flags from pkg-config
AC_SUBST(RSVG_CFLAGS)
AC_SUBST(RSVG_LIBS)



AM_WITH_MYSQL_DATABASE

COOT_CURL

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
doc/Makefile
validation-graphs/Makefile
python/coot_load_modules.py
data/Makefile
data/cho-models/Makefile
data/cho-acedrg/Makefile
data/cho-links/Makefile
data/pdb-templates/Makefile
data/rama-coeffs/Makefile
data/metal/Makefile
monomers/Makefile
glTF/Makefile
pixmaps/Makefile
pixmaps/icons/Makefile
textures/Makefile
sounds/Makefile
rama-data/Makefile
man/Makefile
reference-structures/Makefile
api/Makefile
pli/Makefile
fonts/Makefile
pyrogen/Makefile
pyrogen/pyrogen
layla/Makefile
layla/layla
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
src/Makefile
src/coot
src/coot-relocatable-binaries
ui/Makefile
])
AC_OUTPUT

# removed setup/Makefile
# po/Makefile.in

# 2022-10-28
# rcrane/Makefile
# src/Makefile
# src/coot

# GETTEXT related things: When AM_GNU_GETTEXT get enabled, these get enabled.
# po/Makefile
# intl/Makefile - sigh.  Fix this later, not now.


# or maybe:
#],[sed -e "/POTFILES =/r po/POTFILES" po/Makefile.in > po/Makefile])

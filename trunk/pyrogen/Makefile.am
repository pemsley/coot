

# we need files for make dist configuresd without MAKE_ENHANCED_LIGAND_TOOLS and USE_PYTHON

# pyexec_LTLIBRARIES allows us to make libtool libraries that get the
# file names we like and installed where we like, not forced to be
# libxxx.la in $(libdir) as lib_LTLIBRARIES.
#
pyexec_LTLIBRARIES = 

# rdkit includes are in the rdkit directory, boost includes are the (normal) prefix above that
#
AM_CPPFLAGS = $(CCP4SRS_CXXFLAGS) $(BOOST_CPPFLAGS) $(RDKIT_CXXFLAGS) $(PYTHON_CFLAGS) -I$(top_srcdir)

AM_CXXFLAGS = -DPKGDATADIR='"$(pkgdatadir)"' $(COOT_WINDOWS_CFLAGS)

pyrogen_swig.py: pyrogen-swig.cc

EXTRA_DIST = pyrogen-swig.i

pyrogen-swig.cc: $(top_srcdir)/pyrogen/pyrogen-swig.i $(top_srcdir)/pyrogen/restraints.hh
	swig $(AM_CPPFLAGS) -o pyrogen-swig.cc -python -shadow -c++ $(top_srcdir)/pyrogen/pyrogen-swig.i


if COOT_USE_PYTHON
if MAKE_ENHANCED_LIGAND_TOOLS

pyexec_LTLIBRARIES += _pyrogen_swig.la pyrogen_boost.la
_pyrogen_swig_la_SOURCES = pyrogen-swig.cc restraints.cc restraints.hh

## _pyrogen_swig_la_CFLAGS = $(PYTHON_CFLAGS) -I$(top_srcdir)/include $(PYTHON_CPPFLAGS) 

_pyrogen_swig_la_LIBADD = libcoot_py_restraints.la \
                          $(top_builddir)/ideal/libcoot-ideal.la             \
			  $(top_builddir)/analysis/libcoot-analysis.la       \
			  $(top_builddir)/geometry/libcoot-geometry.la       \
                          $(RDKIT_LIBS) $(L_BOOST_PYTHON) $(CLIPPER_LIBS) $(MMDB_LIBS) \
			  $(PYTHON_LIBS)

# What does this do? If we add it, we don't get this message:
#  `_pyrogen_swig.la' is not a standard libtool library name
# pyrogen/Makefile.am:30: did you mean `lib_pyrogen_swig.la'?
# I doubt the need for the other args, so just use -module for now.
#
# _pyrogen_swig_la_LDFLAGS = -avoid-version -module -shared -export-dynamic
#
_pyrogen_swig_la_LDFLAGS = -module

pyrogen_boost_la_LDFLAGS = -module


# These get installed in libdir.
#
lib_LTLIBRARIES=libcoot_py_restraints.la

# this is a bona fide library
#
libcoot_py_restraints_la_SOURCES = mmff-restraints.cc mmff-restraints.hh \
	py-restraints.cc py-restraints.hh restraints-private.hh

# 20141014 Add $RDKIT_LIBS for Bill
libcoot_py_restraints_la_LIBADD =  \
                                $(top_builddir)/geometry/libcoot-geometry.la \
                                $(RDKIT_LIBS) $(PYTHON_LIBS)

libcoot_py_restraints_la_LDFLAGS = -no-undefined 

# The name of this starts with lib so it fits with boost::python wrapper
pyrogen_boost_la_SOURCES = restraints-boost.cc restraints.hh
pyrogen_boost_la_LIBADD = libcoot_py_restraints.la \
                                $(top_builddir)/ideal/libcoot-ideal.la             \
                                $(top_builddir)/geometry/libcoot-geometry.la       \
                                $(top_builddir)/coords/libcoot-coords.la           \
			        $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
                                $(RDKIT_LIBS) $(L_BOOST_PYTHON)

libpyrogen_boost_la_LIBDFLAGS = -no-undefined

dist_bin_SCRIPTS = pyrogen

dist_pkgpython_PYTHON = pyrogen.py atom_types.py jay_util.py tautomer.py pyrogen_swig.py \
	coot_svn_repo_revision.py

# Something here is broken, I think.
# on $ bash autogen.sh, I get:
# pyrogen/Makefile.am:79: OS_WIN32 does not appear in AM_CONDITIONAL
#
# if OS_WIN32
# # Python .pyd modules are simply DLLs, but they have to be called .pyd for
# # python to find them, and libtool only creates .dll.
# install-exec-hook:
# 	mv $(pyexecdir)/_pyrogen_swig.dll $(pyexecdir)/_pyrogen_swig.pyd
# endif

# clean needs to remove _pyrogen_swig.so
#
clean-local: 
	rm -f _pyrogen_swig.so

endif
endif

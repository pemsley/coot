

# we need files for make dist configuresd without MAKE_ENHANCED_LIGAND_TOOLS and USE_PYTHON

# pyexec_LTLIBRARIES allows us to make libtool libraries that get the
# file names we like and installed where we like, not forced to be
# libxxx.la in $(libdir) as lib_LTLIBRARIES.
#
pyexec_LTLIBRARIES = 

INCLUDES = $(RDKIT_CXXFLAGS) $(PYTHON_CFLAGS) -I$(top_srcdir)

AM_CXXFLAGS = -DPKGDATADIR='"$(pkgdatadir)"'

pyrogen_swig.py: pyrogen-swig.cc

EXTRA_DIST = pyrogen-swig.i

pyrogen-swig.cc: $(top_srcdir)/pyrogen/pyrogen-swig.i $(top_srcdir)/pyrogen/restraints.hh
	swig $(INCLUDES) -o pyrogen-swig.cc -python -shadow -c++ $(top_srcdir)/pyrogen/pyrogen-swig.i


if COOT_USE_PYTHON
if MAKE_ENHANCED_LIGAND_TOOLS

pyexec_LTLIBRARIES += _pyrogen_swig.la libpyrogen_boost.la
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


# These get installed in libdir - libpyrogen_boost.la needs to be
# installed in $(libdir) because now it is a standard libtool library
# not a pyexec lib.
lib_LTLIBRARIES=libcoot_py_restraints.la

# this is a bona fide library
#
libcoot_py_restraints_la_SOURCES = mmff-restraints.cc mmff-restraints.hh \
	py-restraints.cc py-restraints.hh restraints-private.hh
libcoot_py_restraints_la_LIBADD =  $(top_builddir)/geometry/libcoot-geometry.la

# The name of this starts with lib so it fits with boost::python wrapper
libpyrogen_boost_la_SOURCES = restraints-boost.cc restraints.hh
libpyrogen_boost_la_LIBADD = libcoot_py_restraints.la \
                                $(top_builddir)/ideal/libcoot-ideal.la             \
                                $(top_builddir)/geometry/libcoot-geometry.la       \
                                $(top_builddir)/coords/libcoot-coords.la           \
			        $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
                                $(RDKIT_LIBS) $(L_BOOST_PYTHON)


dist_bin_SCRIPTS = pyrogen

dist_pkgpython_SCRIPTS = pyrogen.py jay_util.py tautomer.py pyrogen_swig.py coot_svn_repo_revision.py

# old
# pyrogendir=$(pkgdatadir)/python

# clean needs to remove _pyrogen_swig.so
#
clean-local: 
	rm -f _pyrogen_swig.so

endif
endif

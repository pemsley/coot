

EXTRA_DIST = wrap-coot-libs.i

if COOT_USE_PYTHON
if MAKE_ENHANCED_LIGAND_TOOLS

INCLUDES = $(RDKIT_CXXFLAGS) $(PYTHON_CFLAGS) -I$(top_srcdir) -I$(MMDBSSM_CXXFLAGS)

# I'd like this to be pyrogen_PYTHON
#
pyrogen_SCRIPTS = pyrogen.py jay_util.py tautomer.py

pyrogendir=$(pkgdatadir)/pyrogen

# lib_LIBRARIES=_coot_libs.so restraints_boost.so

lib_LTLIBRARIES=lib_coot_libs.la librestraints_boost.la

lib_coot_libs_la_LIBADD = $(RDKIT_LIBS) $(CLIPPER_LIBS) $(MMDB_LIBS)

lib_coot_libs_la_SOURCES = wrap-coot-libs.cc \
	restraints.cc restraints.hh          \
	py-restraints.cc py-restraints.hh

librestraints_boost_la_SOURCES = restraints-boost.cc restraints.hh

wrap-coot-libs.cc: wrap-coot-libs.i restraints.hh
	swig $(INCLUDES) -o wrap-coot-libs.cc -python -shadow -c++ $(top_srcdir)/pyrogen/wrap-coot-libs.i

# Hmmm....
# _coot_libs.so: lib_coot_libs.so
# cp lib_coot_libs.so _coot_libs.so

endif
endif


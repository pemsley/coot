

EXTRA_DIST = wrap-coot-libs.i

if COOT_USE_PYTHON
if MAKE_ENHANCED_LIGAND_TOOLS

INCLUDES = $(RDKIT_CXXFLAGS) $(PYTHON_CFLAGS) -I$(top_srcdir) -I$(MMDBSSM_CXXFLAGS)

# I'd like this to be pyrogen_PYTHON
#
# these need to go in $(prefix)/share/coot/python
pyrogen_SCRIPTS = pyrogen.py jay_util.py tautomer.py

# pyrogendir=$(pkgdatadir)/pyrogen
pyrogendir=$(pkgdatadir)/python

# lib_LIBRARIES=_coot_libs.so restraints_boost.so

lib_LTLIBRARIES=libcoot_py_restraints.la lib_coot_libs.la librestraints_boost.la

lib_coot_libs_la_LIBADD = libcoot_py_restraints.la \
                          $(top_builddir)/ideal/libcoot-ideal.la             \
			  $(top_builddir)/analysis/libcoot-analysis.la       \
			  $(top_builddir)/geometry/libcoot-geometry.la       \
                          $(RDKIT_LIBS) -lboost_python $(CLIPPER_LIBS) $(MMDB_LIBS)


#                         $(top_builddir)/coot-utils/libcoot-coord-utils.la 


# this is a bone fide library
#
libcoot_py_restraints_la_SOURCES = py-restraints.cc py-restraints.hh

libcoot_py_restraints_la_LIBADD =  $(top_builddir)/geometry/libcoot-geometry.la

lib_coot_libs_la_SOURCES = wrap-coot-libs.cc \
	restraints.cc restraints.hh          

librestraints_boost_la_SOURCES = restraints-boost.cc restraints.hh


librestraints_boost_la_LIBADD = libcoot_py_restraints.la \
                                $(top_builddir)/ideal/libcoot-ideal.la             \
                                $(top_builddir)/geometry/libcoot-geometry.la       \
                                $(top_builddir)/coords/libcoot-coords.la           \
                                $(RDKIT_LIBS) -lboost_python

# I am not sure that you can do this.  This makes two all: targets in the Makefile.
all: all-am _coot_libs.so restraints_boost.so

_coot_libs.so: .libs/lib_coot_libs.so.0.0.0
	cp $^ $@

restraints_boost.so: .libs/librestraints_boost.so.0.0.0
	cp $^ $@

wrap-coot-libs.cc: wrap-coot-libs.i $(top_srcdir)/pyrogen/restraints.hh
	swig $(INCLUDES) -o wrap-coot-libs.cc -python -shadow -c++ $(top_srcdir)/pyrogen/wrap-coot-libs.i

# Hmmm....
# _coot_libs.so: lib_coot_libs.so
# cp lib_coot_libs.so _coot_libs.so

endif
endif


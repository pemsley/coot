
# Python is used by specs.cc/hh
AM_CPPFLAGS = \
	-DUSE_PYTHON \
	-I$(top_srcdir) \
	$(ENHANCED_LIGAND_TOOLS_CPPFLAGS) \
	$(COOT_GSL_FLAGS)		\
	$(GSL_CFLAGS)			\
	$(CLIPPER_CXXFLAGS)             \
	$(MMDB_CXXFLAGS) \
	$(PYTHON_CPPFLAGS)  \
	$(RDKIT_CXXFLAGS) $(BOOST_CPPFLAGS) \
	$(COOT_WINDOWS_CFLAGS)

lib_LTLIBRARIES=libcoot-pli.la

# put specs into a python-loving environment when we need access to coot-pli from python

# remove (python-using) specs from coot-pli
libcoot_pli_la_SOURCES = protein-ligand-interactions.cc flev-annotations.cc pi-stacking.cc \
	sdf-interface.hh         \
	sdf-internal.hh          \
	sdf-interface.cc

libcoot_pli_la_LIBADD = $(top_builddir)/coot-utils/libcoot-coord-utils.la \
                        $(top_builddir)/geometry/libcoot-geometry.la      \
                        $(top_builddir)/lidia-core/libcoot-lidia-core.la  \
                        $(CLIPPER_LIBS) $(MMDB_LIBS)                      \
                        $(RDKIT_LIBS) $(BOOST_LDFLAGS)

libcoot_pli_la_LDFLAGS = $(SHARED_LDFLAGS)

pkginclude_HEADERS = flev-annotations.hh protein-ligand-interactions.hh pi-stacking.hh specs.hh \
	sdf-interface.hh sdf-internal.hh sdf-interface-for-export.hh

pkgincludedir = $(includedir)/coot/pli
